from properties import mea_moles, water_mass_fraction
from eq import eq_loading

DEFAULT_STEAM_TEMPERATURE = 393.15
DEFAULT_PREHEATED_FEED_TEMP = 378.15  # Preheated rich solvent from cross-HEX (~105 °C)
DEFAULT_REBOILER_EFFICIENCY = 0.85
DEFAULT_STRIPPER_P_CO2 = 20000.0      # Typical CO2 partial pressure in stripper bottom (20 kPa)
LATENT_HEAT_WATER_J_KG = 2260000.0
CO2_MOLAR_MASS_G_MOL = 44.01
WATER_MOLAR_MASS_G_MOL = 18.015


def validate_inputs(
    rich_loading,
    lean_loading,
    solvent_flow,
    mea_mass_fraction,
    regeneration_temp,
    feed_temp,
    cp,
    heat_desorption,
    stripping_water_per_co2,
    reboiler_efficiency,
):
    if rich_loading < 0:
        raise ValueError("Rich loading cannot be negative.")
    if lean_loading < 0:
        raise ValueError("Lean loading cannot be negative.")
    if rich_loading <= lean_loading:
        raise ValueError("Rich loading must be strictly greater than lean loading.")
    if solvent_flow <= 0:
        raise ValueError("Solvent flow must be positive.")
    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")
    if regeneration_temp <= 273.15:
        raise ValueError("Regeneration temperature must be above 0°C (273.15 K).")
    if feed_temp <= 273.15:
        raise ValueError("Feed temperature must be above 0°C (273.15 K).")
    if feed_temp >= regeneration_temp:
        raise ValueError("Regeneration temperature must be greater than feed temperature.")
    if cp <= 0:
        raise ValueError("Specific heat capacity must be positive.")
    if heat_desorption <= 0:
        raise ValueError("Heat of desorption must be positive.")
    if stripping_water_per_co2 < 0:
        raise ValueError("Stripping ratio cannot be negative.")
    if not 0 < reboiler_efficiency <= 1:
        raise ValueError("Reboiler efficiency must be between 0 and 1.")


def mea_flow(mea_mass_fraction, solvent_flow):
    return mea_moles(mea_mass_fraction, solvent_flow)


def co2_desorbed(mea_molar_flow_mol_s, rich_loading, lean_loading):
    loading_difference = rich_loading - lean_loading
    co2_desorbed_val = mea_molar_flow_mol_s * loading_difference
    return max(co2_desorbed_val, 0.0)


def co2_flow_kg_s(co2_desorbed_mol_s):
    return co2_desorbed_mol_s * CO2_MOLAR_MASS_G_MOL / 1000.0


def water_evaporation(co2_desorbed_mol_s, stripping_water_per_co2):
    water_evaporated_mol_s = co2_desorbed_mol_s * stripping_water_per_co2
    water_evaporated_kg_s = water_evaporated_mol_s * WATER_MOLAR_MASS_G_MOL / 1000.0
    return (water_evaporated_mol_s, water_evaporated_kg_s)


def reboiler_duty(
    solvent_mass_flow,
    regeneration_temp,
    feed_temp,
    cp,
    co2_desorbed_mol_s,
    heat_desorption,
    stripping_water_per_co2,
    reboiler_efficiency,
):
    sensible_heat_kW = solvent_mass_flow * cp * (regeneration_temp - feed_temp)
    desorption_heat_kW = co2_desorbed_mol_s * heat_desorption / 1000.0

    water_evaporated_mol_s, water_evaporated_kg_s = water_evaporation(
        co2_desorbed_mol_s, stripping_water_per_co2
    )

    evaporation_heat_kW = water_evaporated_kg_s * LATENT_HEAT_WATER_J_KG / 1000.0
    total_heat_kW = sensible_heat_kW + desorption_heat_kW + evaporation_heat_kW
    actual_heat_kW = total_heat_kW / reboiler_efficiency

    return {
        "sensible_heat_kW": sensible_heat_kW,
        "desorption_heat_kW": desorption_heat_kW,
        "evaporation_heat_kW": evaporation_heat_kW,
        "water_evaporated_mol_s": water_evaporated_mol_s,
        "water_evaporated_kg_s": water_evaporated_kg_s,
        "ideal_reboiler_duty_kW": total_heat_kW,
        "actual_reboiler_duty_kW": actual_heat_kW,
    }


def specific_energy(reboiler_duty_kW, co2_desorbed_mol_s):
    if co2_desorbed_mol_s <= 0 or reboiler_duty_kW <= 0:
        return 0.0
    co2_mass_flow_kg_s = co2_desorbed_mol_s * CO2_MOLAR_MASS_G_MOL / 1000.0
    reboiler_duty_MW = reboiler_duty_kW / 1000.0
    return reboiler_duty_MW / co2_mass_flow_kg_s


def simulate_regeneration(
    rich_loading,
    lean_loading=None,
    solvent_flow=0.20,
    mea_mass_fraction=0.30,
    regeneration_temp=DEFAULT_STEAM_TEMPERATURE,
    feed_temp=DEFAULT_PREHEATED_FEED_TEMP,
    cp=4.0,
    heat_desorption=75000.0,
    stripping_water_per_co2=1.2,
    reboiler_efficiency=DEFAULT_REBOILER_EFFICIENCY,
    p_co2_stripper=DEFAULT_STRIPPER_P_CO2,
    pressure=101325.0,
):
    # Dynamically compute thermodynamic lean loading if not explicitly provided
    if lean_loading is None:
        lean_loading = eq_loading(
            t_user=regeneration_temp,
            pressure=pressure,
            Pp_co2=p_co2_stripper,
            mea_mass_fraction=mea_mass_fraction,
            loading=0.0,
        )

    validate_inputs(
        rich_loading=rich_loading,
        lean_loading=lean_loading,
        solvent_flow=solvent_flow,
        mea_mass_fraction=mea_mass_fraction,
        regeneration_temp=regeneration_temp,
        feed_temp=feed_temp,
        cp=cp,
        heat_desorption=heat_desorption,
        stripping_water_per_co2=stripping_water_per_co2,
        reboiler_efficiency=reboiler_efficiency,
    )

    mea_molar_flow_mol_s = mea_flow(mea_mass_fraction, solvent_flow)
    co2_desorbed_mol_s = co2_desorbed(mea_molar_flow_mol_s, rich_loading, lean_loading)
    co2_mass_flow_kg_s = co2_flow_kg_s(co2_desorbed_mol_s)

    heat_results = reboiler_duty(
        solvent_mass_flow=solvent_flow,
        regeneration_temp=regeneration_temp,
        feed_temp=feed_temp,
        cp=cp,
        co2_desorbed_mol_s=co2_desorbed_mol_s,
        heat_desorption=heat_desorption,
        stripping_water_per_co2=stripping_water_per_co2,
        reboiler_efficiency=reboiler_efficiency,
    )

    specific_energy_val = specific_energy(
        heat_results["actual_reboiler_duty_kW"], co2_desorbed_mol_s
    )
    water_fraction = water_mass_fraction(mea_mass_fraction)

    return {
        "rich_loading": rich_loading,
        "lean_loading": lean_loading,
        "loading_difference": rich_loading - lean_loading,
        "MEA_molar_flow_mol_s": mea_molar_flow_mol_s,
        "CO2_desorbed_mol_s": co2_desorbed_mol_s,
        "CO2_desorbed_kg_s": co2_mass_flow_kg_s,
        "water_mass_fraction": water_fraction,
        "regeneration_temperature_K": regeneration_temp,
        "reboiler_duty_kW": heat_results["actual_reboiler_duty_kW"],
        "sensible_heat_kW": heat_results["sensible_heat_kW"],
        "desorption_heat_kW": heat_results["desorption_heat_kW"],
        "evaporation_heat_kW": heat_results["evaporation_heat_kW"],
        "water_evaporated_mol_s": heat_results["water_evaporated_mol_s"],
        "water_evaporated_kg_s": heat_results["water_evaporated_kg_s"],
        "specific_energy_GJ_per_tonne_CO2": specific_energy_val,
    }