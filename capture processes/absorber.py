import math

from eq import eq_loading
from kinetics import reaction_rate
from properties import co2_partial_pressure, mea_moles, water_mass_fraction

DEFAULT_SEGMENTS = 50
DEFAULT_SOLVENT_DENSITY = 1050.0
DEFAULT_CO2_BALANCE_TOLERANCE = 0.05


def gas_liq_transfer(loading_driving_force, mea_conc, mass_coeff, interfacial):

    # FIXED: Moved validation checks above early exit check
    if mea_conc <= 0:
        raise ValueError("MEA concentration must be positive.")

    if mass_coeff <= 0:
        raise ValueError("Mass-transfer coefficient must be positive.")

    if interfacial <= 0:
        raise ValueError("Interfacial area must be positive.")

    if loading_driving_force <= 0:
        return 0.0

    # CO2 concentration driving force:
    # Delta C_CO2 = Delta alpha * C_MEA
    co2_driving_force_mol_m3 = loading_driving_force * mea_conc

    # Mass-transfer capacity:
    # mol/s = K_L * A * Delta C
    transfer_rate = mass_coeff * interfacial * co2_driving_force_mol_m3

    return max(transfer_rate, 0.0)


def simulate_absorber(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    temperature_K,
    pressure_Pa,
    column_height,
    column_area,
    mass_transfer_coefficient,
    specific_interfacial,
    initial_loading,
    number_of_segments=DEFAULT_SEGMENTS,  # FIXED: Assigned default constant
    solvent_density_kg_m3=DEFAULT_SOLVENT_DENSITY,  # FIXED: Assigned default constant
    co2_balance_tolerance=DEFAULT_CO2_BALANCE_TOLERANCE,  # FIXED: Assigned default constant
):

    if not 0 < co2_inlet_fraction < 1:
        raise ValueError("CO2 inlet fraction must be between 0 and 1.")

    if gas_flow <= 0:
        raise ValueError("Gas flow must be positive.")

    if solvent_flow <= 0:
        raise ValueError("Solvent flow must be positive.")

    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")

    if temperature_K <= 273.15:
        raise ValueError("Temperature must be above 0°C.")

    if pressure_Pa <= 0:
        raise ValueError("Pressure must be positive.")

    if column_height <= 0:
        raise ValueError("Column height must be positive.")

    if column_area <= 0:
        raise ValueError("Column area must be positive.")

    if mass_transfer_coefficient <= 0:
        raise ValueError("Mass-transfer coefficient must be positive.")

    if specific_interfacial <= 0:
        raise ValueError("Specific interfacial area must be positive.")

    if initial_loading < 0:
        raise ValueError("Initial CO2 loading cannot be negative.")

    if number_of_segments <= 0:
        raise ValueError("Number of segments must be positive.")

    if solvent_density_kg_m3 <= 0:
        raise ValueError("Solvent density must be positive.")

    if co2_balance_tolerance <= 0:
        raise ValueError("CO2 balance tolerance must be positive.")

    dz = column_height / number_of_segments
    segment_volume = column_area * dz
    total_interfacial = specific_interfacial * segment_volume

    mea_moles_per_kg = mea_moles(mea_mass_fraction)
    water_mass_fraction_value = water_mass_fraction(mea_mass_fraction)
    mea_molar_flow = solvent_flow * mea_moles_per_kg
    solvent_volumetric_flow = solvent_flow / solvent_density_kg_m3
    mea_concentration = mea_molar_flow / max(solvent_volumetric_flow, 1e-12)

    inlet_co2_flow = co2_inlet_fraction * gas_flow
    inert_gas_flow = gas_flow - inlet_co2_flow

    def evaluate_column(bottom_loading):

        current_gas_co2 = inlet_co2_flow
        current_liquid_loading = bottom_loading
        total_absorbed = 0.0

        profile = []
        for segment in range(number_of_segments):

            total_gas = inert_gas_flow + current_gas_co2
            gas_co2_fraction = current_gas_co2 / max(total_gas, 1e-12)
            Pp_co2 = co2_partial_pressure(gas_co2_fraction, pressure_Pa)

            eq_loading_val = eq_loading(
                temperature_K,
                pressure_Pa,
                Pp_co2,
                mea_mass_fraction,
                current_liquid_loading,
            )

            loading_driving_force = max(
                eq_loading_val - current_liquid_loading, 0.0
            )

            transfer_capacity = gas_liq_transfer(
                loading_driving_force,
                mea_concentration,
                mass_transfer_coefficient,
                total_interfacial,
            )

            # -------------------------------------------------
            # Reaction kinetics
            # -------------------------------------------------

            co2_concentration = (
                current_liquid_loading
                * mea_concentration
            )

            available_mea_fraction = max(
                1.0 - 2.0 * current_liquid_loading,
                0.0
            )

            available_mea_concentration = (
                available_mea_fraction
                * mea_concentration
            )

            reaction = reaction_rate(
                co2_concentration=
                    max(co2_concentration, 0.0),
                mea_concentration=
                    max(
                        available_mea_concentration,
                        0.0
                    ),
                temperature_K=temperature_K
            )

            # Reaction capacity in mol/s
            reaction_capacity = (
                reaction * segment_volume
            )

            # -------------------------------------------------
            # Physical limits
            # -------------------------------------------------

            # Maximum CO2 that can be transferred before
            # reaching equilibrium in this segment.
            max_loading_capacity = (
                loading_driving_force
                * mea_molar_flow
            )

            # CO2 transfer is limited by:
            #
            # 1. Mass-transfer capacity
            # 2. Reaction capacity
            # 3. CO2 available in gas
            # 4. Available liquid loading capacity
            #
            co2_removed = min(
                transfer_capacity,
                reaction_capacity,
                current_gas_co2,
                max_loading_capacity
            )

            co2_removed = max(
                co2_removed,
                0.0
            )

            # -------------------------------------------------
            # Store segment information
            # -------------------------------------------------

            profile.append({
                "segment": segment + 1,

                "height_m":
                    (segment + 1) * dz,

                "gas_CO2_fraction":
                    gas_co2_fraction,

                "CO2_partial_pressure_Pa":
                    co2_partial_pressure_value,

                "equilibrium_loading":
                    equilibrium_loading_value,

                "loading":
                    current_liquid_loading,

                "loading_driving_force":
                    loading_driving_force,

                "CO2_removed_mol_s":
                    co2_removed,

                "transfer_capacity":
                    transfer_capacity,

                "reaction_capacity":
                    reaction_capacity,

                "reaction_rate":
                    reaction
            })

            # -------------------------------------------------
            # Update gas stream
            # -------------------------------------------------

            current_gas_co2 -= co2_removed

            current_gas_co2 = max(
                current_gas_co2,
                0.0
            )

            # -------------------------------------------------
            # Update liquid loading
            # -------------------------------------------------

            loading_change = (
                co2_removed
                / max(mea_molar_flow, 1e-12)
            )

            current_liquid_loading -= (
                loading_change
            )

            current_liquid_loading = max(
                current_liquid_loading,
                0.0
            )

            total_absorbed += co2_removed

        # -----------------------------------------------------
        # Shooting-method boundary error
        # -----------------------------------------------------

        error = (
            current_liquid_loading
            - initial_loading
        )

        return (
            error,
            profile,
            current_gas_co2,
            total_absorbed
        )

    # ---------------------------------------------------------
    # 6. Bisection search for rich loading
    # ---------------------------------------------------------

    low_guess = initial_loading

    high_guess = 0.50

    mid_guess = low_guess

    for _ in range(50):

        mid_guess = (
            low_guess + high_guess
        ) / 2.0

        error, _, _, _ = (
            evaluate_column_profile(
                mid_guess
            )
        )

        if abs(error) < 1e-5:
            break

        if error > 0:

            # Calculated top loading is too high.
            # Reduce bottom/rich loading.
            high_guess = mid_guess

        else:

            # Calculated top loading is too low.
            # Increase bottom/rich loading.
            low_guess = mid_guess

    # ---------------------------------------------------------
    # 7. Final absorber calculation
    # ---------------------------------------------------------

    (
        final_error,
        final_profile,
        final_gas_co2,
        total_co2_absorbed
    ) = evaluate_column_profile(
        mid_guess
    )

    # ---------------------------------------------------------
    # 8. Independent CO2 calculation from solvent loading
    #
    # At steady state:
    #
    # CO2 absorbed =
    # MEA flow × (rich loading - lean loading)
    # ---------------------------------------------------------

    loading_based_co2_absorbed = (
        mea_molar_flow
        * max(
            mid_guess - initial_loading,
            0.0
        )
    )

    # ---------------------------------------------------------
    # 9. CO2 balance check
    # ---------------------------------------------------------

    co2_balance_error = (
        abs(
            total_co2_absorbed
            - loading_based_co2_absorbed
        )
        / max(
            total_co2_absorbed,
            1e-12
        )
    )

    # ---------------------------------------------------------
    # 10. Final overall gas-side metrics
    # ---------------------------------------------------------

    capture_percentage = (
        total_co2_absorbed
        / inlet_co2_flow
        * 100.0
        if inlet_co2_flow > 0
        else 0.0
    )

    final_total_gas = (
        inert_gas_flow
        + final_gas_co2
    )

    outlet_co2_fraction = (
        final_gas_co2
        / max(final_total_gas, 1e-12)
    )

    # ---------------------------------------------------------
    # 11. Return results
    # ---------------------------------------------------------

    return {
        "capture_percentage":
            capture_percentage,

        "outlet_co2_fraction":
            outlet_co2_fraction,

        "co2_absorbed_mol_s":
            total_co2_absorbed,

        "loading_based_co2_absorbed_mol_s":
            loading_based_co2_absorbed,

        "co2_balance_error":
            co2_balance_error,

        "co2_balance_error_percentage":
            co2_balance_error * 100.0,

        "absorber_balance_passed":
            co2_balance_error
            <= co2_balance_tolerance,

        "rich_loading":
            mid_guess,

        "MEA_molar_flow":
            mea_molar_flow,

        "water_mass_fraction":
            water_mass_fraction_value,

        "profile":
            final_profile
    }
