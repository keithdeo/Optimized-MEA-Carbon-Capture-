import math

from eq import eq_loading
from kinetics import reaction_rate
from properties import co2_partial_pressure, mea_moles, water_mass_fraction

DEFAULT_SEGMENTS = float(input("Enter the number of absorber segments (between 5.0 and 500.0): "))
DEFAULT_SOLVENT_DENSITY = float(input("Enter the solvent density (kg/m³, between 950.0 and 1250.0): "))
DEFAULT_CO2_BALANCE_TOLERANCE = float(input("Enter the CO2 balance tolerance (%, between 0.1% and 20%): "))
DEFAULT_CO2_BALANCE_TOLERANCE = DEFAULT_CO2_BALANCE_TOLERANCE / 100.0 

if not 5.0 <= DEFAULT_SEGMENTS <= 500.0:
    raise ValueError("Number of absorber segments must be between 5.0 and 500.0.")

if not 950.0 <= DEFAULT_SOLVENT_DENSITY <= 1250.0:
    raise ValueError("Solvent density must be between 950.0 and 1250.0 kg/m³.")

if not 0.001 <= DEFAULT_CO2_BALANCE_TOLERANCE <= 0.20:
    raise ValueError("CO2 balance tolerance must be between 0.1% and 20%.")

if not 0 < DEFAULT_CO2_BALANCE_TOLERANCE <= 0.20:
    raise ValueError("CO2 balance tolerance must be between 0.1% and 20%.")

def gas_liq_transfer(loading_driving_force, mea_conc, mass_coeff, interfacial):

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
    temperature_k,
    pressure_Pa,
    column_height,
    column_area,
    mass_transfer_coefficient,
    specific_interfacial,
    initial_loading,
    number_of_segments=DEFAULT_SEGMENTS,
    solvent_density_kg_m3=DEFAULT_SOLVENT_DENSITY,  
    co2_balance_tolerance=DEFAULT_CO2_BALANCE_TOLERANCE,  
):

    if not 0 < co2_inlet_fraction < 1:
        raise ValueError("CO2 inlet fraction must be between 0 and 1.")

    if gas_flow <= 0:
        raise ValueError("Gas flow must be positive.")

    if solvent_flow <= 0:
        raise ValueError("Solvent flow must be positive.")

    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")

    if temperature_k <= 273.15:
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
                temperature_k,
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

            co2_concentration = current_liquid_loading * mea_concentration
            
            available_mea_fraction = max(1.0 - 2.0 * current_liquid_loading,0.0)
            available_mea_concentration = available_mea_fraction * mea_concentration

            reaction = reaction_rate(temperature_k,max(co2_concentration, 0.0),
                    max(available_mea_concentration,0.0))

            reaction_capacity = reaction * segment_volume
            max_loading_capacity = loading_driving_force * mea_molar_flow

            co2_removed = min(transfer_capacity,reaction_capacity,
                current_gas_co2,max_loading_capacity)
            co2_removed = max(co2_removed,0.0)

            profile.append({
                "segment":segment + 1,
                "height_m":(segment + 1) * dz,
                "gas_CO2_fraction":gas_co2_fraction,
                "co2_partial_pressure":Pp_co2,
                "equilibrium_loading":eq_loading_val,
                "loading":current_liquid_loading,
                "loading_driving_force":loading_driving_force,
                "co2_removed_mol_s":co2_removed,
                "transfer_capacity":transfer_capacity,
                "reaction_capacity":reaction_capacity,
                "reaction_rate":reaction})

            current_gas_co2 -= co2_removed
            current_gas_co2 = max(current_gas_co2,0.0)

            loading_change = co2_removed / max(mea_molar_flow, 1e-12)
            current_liquid_loading -= loading_change
            current_liquid_loading = max(current_liquid_loading,0.0)
            total_absorbed += co2_removed

        error = current_liquid_loading - initial_loading

        return (error,profile,current_gas_co2,total_absorbed)

    low_guess = initial_loading
    high_guess = 0.50

    for _ in range(50):

        mid_guess = (low_guess + high_guess) / 2.0
        error, _, _, _ = (evaluate_column(mid_guess))

        if abs(error) < 1e-5:
            break
        if error > 0:
            # Calculated top loading is too high.
            high_guess = mid_guess
        else:
            # Calculated top loading is too low.
            low_guess = mid_guess

    (final_error,final_profile,
    final_gas_co2,total_co2_absorbed) = evaluate_column(mid_guess)

    loading_based_co2_absorbed = (mea_molar_flow *
    max(mid_guess - initial_loading,0.0))
    
    co2_balance_error = (
        abs(total_co2_absorbed - loading_based_co2_absorbed)
        / max(total_co2_absorbed,1e-12))

    if inlet_co2_flow > 0:
        capture_percentage = (
            total_co2_absorbed / inlet_co2_flow * 100.0)
    else:
        capture_percentage = 0.0
    
    final_total_gas = inert_gas_flow + final_gas_co2
    outlet_co2_fraction = final_gas_co2 / max(final_total_gas, 1e-12)
    
    return {
        "capture_percentage":capture_percentage,
        "outlet_co2_fraction":outlet_co2_fraction,
        "co2_absorbed_mol_s":total_co2_absorbed,
        "loading_based_co2_absorbed_mol_s":loading_based_co2_absorbed,
        "co2_balance_error":co2_balance_error,
        "co2_balance_error_percentage":co2_balance_error * 100.0,
        "absorber_balance_passed":co2_balance_error <= co2_balance_tolerance,
        "rich_loading":mid_guess,
        "mea_molar_flow":mea_molar_flow,
        "water_mass_fraction":water_mass_fraction_value,
        "profile":final_profile}

if __name__ == "__main__":

    co2_inlet_fraction = float(
        input("Enter the CO2 inlet fraction (between 0 and 1): ")
    )
    if not (0 < co2_inlet_fraction < 1):
        raise ValueError("CO2 inlet fraction must be strictly between 0 and 1.")

    gas_flow = float(input("Enter the gas flow (mol/s, positive): "))
    if gas_flow <= 0:
        raise ValueError("Gas flow must be positive.")

    solvent_flow = float(input("Enter the solvent flow (mol/s, positive): "))
    if solvent_flow <= 0:
        raise ValueError("Solvent flow must be positive.")

    mea_mass_fraction = float(
    input("Enter the MEA mass fraction (between 0 and 1): ")
)
    if not (0 < mea_mass_fraction < 1):
        raise ValueError("MEA mass fraction must be strictly between 0 and 1.")

    temperature_k = float(input("Enter the temperature (K, above 273.15): "))
    if temperature_k <= 273.15:
        raise ValueError("Temperature must be above 273.15 K.")

    pressure_Pa = float(input("Enter the pressure (Pa, positive): "))
    if pressure_Pa <= 0:
        raise ValueError("Pressure must be positive.")

    column_height = float(input("Enter the column height (m, positive): "))
    if column_height <= 0:
        raise ValueError("Column height must be positive.")

    column_area = float(input("Enter the column area (m², positive): "))
    if column_area <= 0:
        raise ValueError("Column area must be positive.")

    mass_transfer_coefficient = float(
    input("Enter the mass-transfer coefficient (m/s, positive): ")
)
    if mass_transfer_coefficient <= 0:
        raise ValueError("Mass-transfer coefficient must be positive.")

    specific_interfacial = float(
    input("Enter the specific interfacial area (m²/m³, positive): ")
)
    if specific_interfacial <= 0:
        raise ValueError("Specific interfacial area must be positive.")

    initial_loading = float(
    input("Enter the initial CO2 loading (mol CO2/mol MEA, positive): ")
)
    if initial_loading < 0:
        raise ValueError("Initial CO2 loading cannot be negative.")

    results = simulate_absorber(
        co2_inlet_fraction,
        gas_flow,
        solvent_flow,
        mea_mass_fraction,
        temperature_k,
        pressure_Pa,
        column_height,
        column_area,
        mass_transfer_coefficient,
        specific_interfacial,
        initial_loading,
        DEFAULT_SEGMENTS,
        DEFAULT_SOLVENT_DENSITY,
        DEFAULT_CO2_BALANCE_TOLERANCE
    )

    print(f"\nCO2 capture: {results['capture_percentage']:.2f}%")
    print(f"Outlet CO2: {results['outlet_co2_fraction'] * 100:.2f}%")
    print(f"CO2 absorbed: {results['co2_absorbed_mol_s']:.4f} mol/s")
    print(f"Rich loading: {results['rich_loading']:.4f} mol CO2/mol MEA")
