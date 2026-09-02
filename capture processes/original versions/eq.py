import math

# Reference states
R = 8.314462618  # Universal gas constant, J/(mol·K)
T_REF = 298.15  # Baseline reference temperature in Kelvin (25.0 °C)

K_REF_298 = {
    "K_HYDRATION": 1.70e-3,  # CO2(aq) + H2O <-> H2CO3
    "K_1": 2.51e-4,  # H2CO3 <-> H+ + HCO3-
    "K_2": 4.69e-11,  # HCO3- <-> H+ + CO3--
    "K_H2O": 1.008e-14,  # H2O <-> H+ + OH-
    "K_MEA": 1.023e-10,  # MEAH+ <-> MEA + H+
    "K_CARBAMATE": 42.5,  # CO2 + MEA <-> MEACOO-
}

# Standard reaction enthalpies (J/mol at 298.15 K)
DELTA_H = {
    "H_HYDRATION": 19000.0,  # CO2(aq) + H2O <-> H2CO3     [endothermic]
    "H_1": 9100.0,  # H2CO3 <-> H+ + HCO3-          [endothermic]
    "H_2": 14900.0,  # HCO3- <-> H+ + CO3--          [endothermic]
    "H_H2O": 55800.0,  # H2O <-> H+ + OH-              [endothermic]
    "H_MEA": 52000.0,  # MEAH+ <-> MEA + H+           [endothermic]
    "H_CARBAMATE": -48000.0,  # CO2 + MEA <-> MEACOO-         [exothermic]
}


def user_temp(temp_celsius):
    t_user = temp_celsius + 273.15
    return t_user


def calculate_k(t_user, k_298, deltah):
    exponent = -(deltah / R) * ((1.0 / t_user) - (1.0 / T_REF))
    calibrated_k = k_298 * math.exp(exponent)
    return calibrated_k


def equilibrium_constants(t_user):

    k_hydration = calculate_k(
        t_user, K_REF_298["K_HYDRATION"], DELTA_H["H_HYDRATION"]
    )
    k_1 = calculate_k(t_user, K_REF_298["K_1"], DELTA_H["H_1"])
    k_2 = calculate_k(t_user, K_REF_298["K_2"], DELTA_H["H_2"])
    k_h2o = calculate_k(t_user, K_REF_298["K_H2O"], DELTA_H["H_H2O"])
    k_mea = calculate_k(t_user, K_REF_298["K_MEA"], DELTA_H["H_MEA"])
    k_carbamate = calculate_k(
        t_user, K_REF_298["K_CARBAMATE"], DELTA_H["H_CARBAMATE"]
    )

    return {
        "k_hydration": k_hydration,
        "k_1": k_1,
        "k_2": k_2,
        "k_h2o": k_h2o,
        "k_mea": k_mea,
        "k_carbamate": k_carbamate,
    }


def carbonate_species(t_user, co2_conc, h_conc):
    if co2_conc < 0:
        raise ValueError("CO2 concentration cannot be negative.")

    if h_conc <= 0:
        raise ValueError("Hydrogen-ion concentration cannot be negative.")

    constants = equilibrium_constants(t_user)

    K_hydration = constants["k_hydration"]
    K1 = constants["k_1"]
    K2 = constants["k_2"]

    H2CO3 = K_hydration * co2_conc
    HCO3 = K1 * (H2CO3 / h_conc)
    CO3 = K2 * (HCO3 / h_conc)

    return {"CO2": co2_conc, "H2CO3": H2CO3, "HCO3-": HCO3, "CO3--": CO3}


def mea_species(t_user, total_mea, h_conc):
    # [MEAH+] = [MEA][H+]/Ka
    if total_mea < 0:
        raise ValueError("Total MEA concentration cannot be negative.")

    if h_conc <= 0:
        raise ValueError("Hydrogen-ion concentration must be positive.")

    constants = equilibrium_constants(t_user)
    Ka_MEA = constants["k_mea"]  # FIXED: Key lowercase 'k_mea'

    # [MEAH+] = [MEA][H+]/Ka
    free_mea = total_mea / (1 + h_conc / Ka_MEA)
    MEAH = total_mea - free_mea

    return {"MEA": free_mea, "MEAH+": MEAH}


def carbamate(t_user, co2_conc, free_mea):
    # CO2 + MEA <-> MEACOO-
    if co2_conc < 0:
        raise ValueError("CO2 concentration cannot be negative.")

    if free_mea < 0:
        raise ValueError("Free MEA concentration cannot be negative.")

    constants = equilibrium_constants(t_user)
    K_carbamate = constants["k_carbamate"]
    carbamate_conc = K_carbamate * co2_conc * free_mea

    return carbamate_conc


def hydroxide_conc(t_user, h_conc):
    # Kw = [H+][OH-]
    if h_conc <= 0:
        raise ValueError("Hydrogen-ion concentration must be positive.")

    constants = equilibrium_constants(t_user)
    Kw = constants["k_h2o"]

    return Kw / h_conc


def equilibrium_state(t_user, total_co2, total_mea, h_conc):
    # Calculate the chemical species present in the MEA-CO2-H2O system.
    # total_co2 = Total dissolved CO2 concentration, mol/L.
    # total_MEA = Total MEA concentration, mol/L.
    # H_conc = Hydrogen-ion concentration, mol/L.
    if total_co2 < 0:
        raise ValueError("Total CO2 concentration cannot be negative.")

    if total_mea <= 0:
        raise ValueError("Total MEA concentration must be greater than zero.")

    if h_conc <= 0:
        raise ValueError("Hydrogen-ion concentration must be positive.")

    carbonate = carbonate_species(t_user, total_co2, h_conc)
    mea = mea_species(t_user, total_mea, h_conc)
    carbamate_val = carbamate(
        t_user, total_co2, mea["MEA"]
    )  # FIXED: Renamed variable to avoid shadowing function
    OH = hydroxide_conc(t_user, h_conc)

    return {
        **carbonate,
        **mea,
        "MEACOO-": carbamate_val,
        "H+": h_conc,
        "OH-": OH,
    }


def eq_loading(t_user, pressure, Pp_co2, mea_mass_fraction, loading=0.0):
    if t_user <= 0:
        raise ValueError("Temperature must be greater than zero.")

    if pressure <= 0:
        raise ValueError("Pressure must be greater than zero.")

    if not 0 <= Pp_co2 <= pressure:
        raise ValueError(
            "CO2 partial pressure must be between zero and total pressure."
        )

    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")

    pressure_factor = Pp_co2 / (Pp_co2 + 10000.0)
    concentration_factor = mea_mass_fraction / 0.30
    temperature_factor = 313.15 / t_user

    return min(
        0.50,
        max(
            loading,
            0.50 * pressure_factor * concentration_factor * temperature_factor,
        ),
    )
