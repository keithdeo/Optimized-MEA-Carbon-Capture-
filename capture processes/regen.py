from properties import (
    water_mass_fraction,
    mea_moles_per_kg_solution
)

# Default operating values for the regeneration model.
DEFAULT_STEAM_TEMPERATURE = 393.15
DEFAULT_REBOILER_EFFICIENCY = 0.85
LATENT_HEAT_WATER_J_KG = 2260000.0
CO2_MOLAR_MASS_G_MOL = 44.01
WATER_MOLAR_MASS_KG_MOL = 0.018015


# Checks that all regeneration inputs are physically valid.
def validate_inputs(
    rich_loading,
    lean_loading,
    solvent_flow_kg_s,
    mea_mass_fraction,
    regeneration_temperature_K,
    feed_temperature_K,
    cp_solution_kJ_kgK,
    heat_of_desorption_J_mol,
    stripping_ratio_mol_water_per_mol_co2,
    reboiler_efficiency
):

    if rich_loading < 0:
        raise ValueError(
            "Rich loading cannot be negative."
        )

    if lean_loading < 0:
        raise ValueError(
            "Lean loading cannot be negative."
        )

    if rich_loading <= lean_loading:
        raise ValueError(
            "Rich loading must be strictly greater than lean loading."
        )

    if solvent_flow_kg_s <= 0:
        raise ValueError(
            "Solvent flow must be positive."
        )

    if not 0 < mea_mass_fraction < 1:
        raise ValueError(
            "MEA mass fraction must be between 0 and 1."
        )

    if regeneration_temperature_K <= 273.15:
        raise ValueError(
            "Regeneration temperature must be above 0°C (273.15 K)."
        )

    if feed_temperature_K <= 273.15:
        raise ValueError(
            "Feed temperature must be above 0°C (273.15 K)."
        )

    if feed_temperature_K >= regeneration_temperature_K:
        raise ValueError(
            "Regeneration temperature must be greater than feed temperature."
        )

    if cp_solution_kJ_kgK <= 0:
        raise ValueError(
            "Specific heat capacity must be positive."
        )

    if heat_of_desorption_J_mol <= 0:
        raise ValueError(
            "Heat of desorption must be positive."
        )

    if stripping_ratio_mol_water_per_mol_co2 < 0:
        raise ValueError(
            "Stripping ratio cannot be negative."
        )

    if not 0 < reboiler_efficiency <= 1:
        raise ValueError(
            "Reboiler efficiency must be between 0 and 1."
        )


# Calculates the total molar flow of MEA entering the stripper.
def calculate_mea_flow(
    solvent_flow_kg_s,
    mea_mass_fraction
):

    mea_moles_per_kg = (
        mea_moles_per_kg_solution(
            mea_mass_fraction
        )
    )

    return (
        solvent_flow_kg_s
        * mea_moles_per_kg
    )


# Calculates how much CO2 is released during regeneration.
def calculate_co2_desorbed(
    mea_molar_flow_mol_s,
    rich_loading,
    lean_loading
):

    loading_difference = (
        rich_loading
        - lean_loading
    )

    co2_desorbed = (
        mea_molar_flow_mol_s
        * loading_difference
    )

    return max(
        co2_desorbed,
        0.0
    )


# Converts the CO2 molar flow into kg/s.
def calculate_co2_mass_flow(
    co2_desorbed_mol_s
):

    return (
        co2_desorbed_mol_s
        * CO2_MOLAR_MASS_G_MOL
        / 1000.0
    )


# Calculates the amount of water vaporized during stripping.
def calculate_water_evaporation(
    co2_desorbed_mol_s,
    stripping_ratio_mol_water_per_mol_co2
):

    water_evaporated_mol_s = (
        co2_desorbed_mol_s
        * stripping_ratio_mol_water_per_mol_co2
    )

    water_evaporated_kg_s = (
        water_evaporated_mol_s
        * WATER_MOLAR_MASS_KG_MOL
    )

    return (
        water_evaporated_mol_s,
        water_evaporated_kg_s
    )


# Calculates the different components of reboiler heat duty.
def calculate_reboiler_duty(
    solvent_mass_flow_kg_s,
    regeneration_temperature_K,
    feed_temperature_K,
    cp_solution_kJ_kgK,
    co2_desorbed_mol_s,
    heat_of_desorption_J_mol,
    stripping_ratio_mol_water_per_mol_co2,
    reboiler_efficiency
):

    sensible_heat_kW = (
        solvent_mass_flow_kg_s
        * cp_solution_kJ_kgK
        * (
            regeneration_temperature_K
            - feed_temperature_K
        )
    )

    desorption_heat_kW = (
        co2_desorbed_mol_s
        * heat_of_desorption_J_mol
        / 1000.0
    )

    (
        water_evaporated_mol_s,
        water_evaporated_kg_s
    ) = calculate_water_evaporation(
        co2_desorbed_mol_s,
        stripping_ratio_mol_water_per_mol_co2
    )

    evaporation_heat_kW = (
        water_evaporated_kg_s
        * LATENT_HEAT_WATER_J_KG
        / 1000.0
    )

    total_ideal_heat_kW = (
        sensible_heat_kW
        + desorption_heat_kW
        + evaporation_heat_kW
    )

    actual_heat_kW = (
        total_ideal_heat_kW
        / reboiler_efficiency
    )

    return {
        "sensible_heat_kW":
            sensible_heat_kW,

        "desorption_heat_kW":
            desorption_heat_kW,

        "evaporation_heat_kW":
            evaporation_heat_kW,

        "water_evaporated_mol_s":
            water_evaporated_mol_s,

        "water_evaporated_kg_s":
            water_evaporated_kg_s,

        "ideal_reboiler_duty_kW":
            total_ideal_heat_kW,

        "actual_reboiler_duty_kW":
            actual_heat_kW
    }


# Calculates the energy required per tonne of CO2 regenerated.
def calculate_specific_energy(
    reboiler_duty_kW,
    co2_desorbed_mol_s
):

    if co2_desorbed_mol_s <= 0:
        return 0.0

    co2_mass_flow_tonne_h = (
        co2_desorbed_mol_s
        * CO2_MOLAR_MASS_G_MOL
        * 0.0036
    )

    if co2_mass_flow_tonne_h <= 0:
        return 0.0

    energy_GJ_per_tonne = (
        reboiler_duty_kW
        * 0.0036
        / co2_mass_flow_tonne_h
    )

    return energy_GJ_per_tonne


# Runs the complete regeneration calculation.
def simulate_regeneration(
    rich_loading,
    lean_loading,
    solvent_flow_kg_s,
    mea_mass_fraction,
    regeneration_temperature_K,
    feed_temperature_K,
    cp_solution_kJ_kgK=4.0,
    heat_of_desorption_J_mol=75000.0,
    stripping_ratio_mol_water_per_mol_co2=1.2,
    reboiler_efficiency=DEFAULT_REBOILER_EFFICIENCY
):

    validate_inputs(
        rich_loading=rich_loading,
        lean_loading=lean_loading,
        solvent_flow_kg_s=solvent_flow_kg_s,
        mea_mass_fraction=mea_mass_fraction,
        regeneration_temperature_K=regeneration_temperature_K,
        feed_temperature_K=feed_temperature_K,
        cp_solution_kJ_kgK=cp_solution_kJ_kgK,
        heat_of_desorption_J_mol=heat_of_desorption_J_mol,
        stripping_ratio_mol_water_per_mol_co2=
            stripping_ratio_mol_water_per_mol_co2,
        reboiler_efficiency=reboiler_efficiency
    )

    mea_molar_flow_mol_s = (
        calculate_mea_flow(
            solvent_flow_kg_s,
            mea_mass_fraction
        )
    )

    co2_desorbed_mol_s = (
        calculate_co2_desorbed(
            mea_molar_flow_mol_s,
            rich_loading,
            lean_loading
        )
    )

    co2_mass_flow_kg_s = (
        calculate_co2_mass_flow(
            co2_desorbed_mol_s
        )
    )

    heat_results = (
        calculate_reboiler_duty(
            solvent_mass_flow_kg_s=
                solvent_flow_kg_s,

            regeneration_temperature_K=
                regeneration_temperature_K,

            feed_temperature_K=
                feed_temperature_K,

            cp_solution_kJ_kgK=
                cp_solution_kJ_kgK,

            co2_desorbed_mol_s=
                co2_desorbed_mol_s,

            heat_of_desorption_J_mol=
                heat_of_desorption_J_mol,

            stripping_ratio_mol_water_per_mol_co2=
                stripping_ratio_mol_water_per_mol_co2,

            reboiler_efficiency=
                reboiler_efficiency
        )
    )

    specific_energy = (
        calculate_specific_energy(
            heat_results[
                "actual_reboiler_duty_kW"
            ],
            co2_desorbed_mol_s
        )
    )

    water_fraction = (
        water_mass_fraction(
            mea_mass_fraction
        )
    )

    return {

        "rich_loading":
            rich_loading,

        "lean_loading":
            lean_loading,

        "loading_difference":
            rich_loading
            - lean_loading,

        "MEA_molar_flow_mol_s":
            mea_molar_flow_mol_s,

        "CO2_desorbed_mol_s":
            co2_desorbed_mol_s,

        "CO2_desorbed_kg_s":
            co2_mass_flow_kg_s,

        "water_mass_fraction":
            water_fraction,

        "regeneration_temperature_K":
            regeneration_temperature_K,

        "reboiler_duty_kW":
            heat_results[
                "actual_reboiler_duty_kW"
            ],

        "sensible_heat_kW":
            heat_results[
                "sensible_heat_kW"
            ],

        "desorption_heat_kW":
            heat_results[
                "desorption_heat_kW"
            ],

        "evaporation_heat_kW":
            heat_results[
                "evaporation_heat_kW"
            ],

        "water_evaporated_mol_s":
            heat_results[
                "water_evaporated_mol_s"
            ],

        "water_evaporated_kg_s":
            heat_results[
                "water_evaporated_kg_s"
            ],

        "specific_energy_GJ_per_tonne_CO2":
            specific_energy
    }
