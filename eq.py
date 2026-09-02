import math

# Reference constants
R = 8.314462618  # Universal gas constant, J/(mol·K)
T_REF = 298.15   # Reference temperature (25 °C) in Kelvin

K_REF_298 = {
    "K_HYDRATION": 1.70e-3,   # CO2(aq) + H2O <-> H2CO3
    "K_1": 2.51e-4,           # H2CO3 <-> H+ + HCO3-
    "K_2": 4.69e-11,          # HCO3- <-> H+ + CO3--
    "K_H2O": 1.008e-14,       # H2O <-> H+ + OH-
    "K_MEA": 1.023e-10,       # MEAH+ <-> MEA + H+
    "K_CARBAMATE": 42.5,      # CO2 + MEA <-> MEACOO-
}

DELTA_H = {
    "H_HYDRATION": 19000.0,
    "H_1": 9100.0,
    "H_2": 14900.0,
    "H_H2O": 55800.0,
    "H_MEA": 52000.0,
    "H_CARBAMATE": -48000.0,
}

def user_temp(temp):
    """Ensures temperature is returned in Kelvin."""
    return temp + 273.15 if temp < 200.0 else temp

def calculate_k(t_user, k_298, deltah):
    t_k = user_temp(t_user)
    exponent = -(deltah / R) * ((1.0 / t_k) - (1.0 / T_REF))
    return k_298 * math.exp(exponent)

def equilibrium_constants(t_user):
    t_k = user_temp(t_user)
    return {
        "k_hydration": calculate_k(t_k, K_REF_298["K_HYDRATION"], DELTA_H["H_HYDRATION"]),
        "k_1": calculate_k(t_k, K_REF_298["K_1"], DELTA_H["H_1"]),
        "k_2": calculate_k(t_k, K_REF_298["K_2"], DELTA_H["H_2"]),
        "k_h2o": calculate_k(t_k, K_REF_298["K_H2O"], DELTA_H["H_H2O"]),
        "k_mea": calculate_k(t_k, K_REF_298["K_MEA"], DELTA_H["H_MEA"]),
        "k_carbamate": calculate_k(t_k, K_REF_298["K_CARBAMATE"], DELTA_H["H_CARBAMATE"]),
    }

def carbonate_species(t_user, co2_conc, h_conc):
    if co2_conc < 0: raise ValueError("CO2 concentration cannot be negative.")
    if h_conc <= 0: raise ValueError("Hydrogen-ion concentration must be positive.")

    constants = equilibrium_constants(t_user)
    H2CO3 = constants["k_hydration"] * co2_conc
    HCO3 = constants["k_1"] * (H2CO3 / h_conc)
    CO3 = constants["k_2"] * (HCO3 / h_conc)

    return {"CO2": co2_conc, "H2CO3": H2CO3, "HCO3-": HCO3, "CO3--": CO3}

def mea_species(t_user, total_mea, h_conc):
    if total_mea < 0: raise ValueError("Total MEA concentration cannot be negative.")
    if h_conc <= 0: raise ValueError("Hydrogen-ion concentration must be positive.")

    constants = equilibrium_constants(t_user)
    Ka_MEA = constants["k_mea"]
    free_mea = total_mea / (1.0 + h_conc / Ka_MEA)
    MEAH = total_mea - free_mea

    return {"MEA": free_mea, "MEAH+": MEAH}

def carbamate(t_user, co2_conc, free_mea):
    if co2_conc < 0: raise ValueError("CO2 concentration cannot be negative.")
    if free_mea < 0: raise ValueError("Free MEA concentration cannot be negative.")

    constants = equilibrium_constants(t_user)
    return constants["k_carbamate"] * co2_conc * free_mea

def hydroxide_conc(t_user, h_conc):
    if h_conc <= 0: raise ValueError("Hydrogen-ion concentration must be positive.")
    return equilibrium_constants(t_user)["k_h2o"] / h_conc

def eq_loading(t_user, pressure, Pp_co2, mea_mass_fraction, loading=0.0):
    """
    Thermodynamic VLE equilibrium model for CO2-MEA-H2O system.
    Returns maximum loading in mol CO2 / mol MEA.
    """
    if t_user <= 0: raise ValueError("Temperature must be positive.")
    if pressure <= 0: raise ValueError("Pressure must be positive.")
    if not 0 <= Pp_co2 <= pressure: raise ValueError("CO2 partial pressure out of range.")
    if not 0 < mea_mass_fraction < 1: raise ValueError("MEA mass fraction must be between 0 and 1.")

    if Pp_co2 <= 1e-6:
        return 0.0

    t_k = user_temp(t_user)

    # Henry's Law constant for CO2 in solution, Pa/(mol/L)
    H_co2 = 2.94e6 * math.exp(-2400.0 * (1.0 / t_k - 1.0 / 298.15))
    co2_aq = Pp_co2 / max(H_co2, 1e-12)

    # Carbamate equilibrium contribution (Kent-Eisenberg empirical VLE fit)
    ln_Kc = -11.6 + 6000.0 / t_k
    Kc = math.exp(ln_Kc)
    Y = math.sqrt(max(Kc * co2_aq, 0.0))
    alpha_carb = Y / (1.0 + 2.0 * Y)

    # Bicarbonate formation contribution at higher CO2 partial pressures
    ln_Kbic = -4.0 + 2800.0 / t_k
    Kbic = math.exp(ln_Kbic)
    alpha_bic = (Kbic * co2_aq) / (1.0 + Kbic * co2_aq) * 0.15

    # Concentration scaling normalized to 30 wt% MEA
    conc_factor = (mea_mass_fraction / 0.30) ** 0.1
    alpha_total = (alpha_carb + alpha_bic) * conc_factor

    return max(0.0, min(0.58, alpha_total))