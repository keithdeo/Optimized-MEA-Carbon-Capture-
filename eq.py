import math

from properties import R

# Reference states

# temperature
T_REF = 298.15       # K


# CO2 hydration
# CO2 + H2O <-> H2CO3
K_HYDRATION_REF = 1.7e-3


# Carbonic acid first dissociation
# H2CO3 <-> H+ + HCO3-
K1_REF = 2.5e-4


# Carbonic acid second dissociation
# HCO3- <-> H+ + CO3--
K2_REF = 4.7e-11


# Water dissociation
# H2O <-> H+ + OH-
KW_REF = 1.07e-14


# MEA protonation
# MEAH+ <-> MEA + H+
MEA_ACIDITY_REF = 1.02e-10


# MEA carbamate formation
# CO2 + MEA <-> MEACOO-
CARBAMATE_K_REF = 40.0


# TEMPERATURE CORRELATION

def temperature_equilibrium_constant(
    K_ref,
    delta_H,
    temperature_K
):
    """
    Estimate an equilibrium constant at a different
    temperature using a van't Hoff relationship.

        ln(K2/K1)
        =
        -(delta_H/R)(1/T2 - 1/T1)

    Parameters
    ----------
    K_ref : float
        Equilibrium constant at T_REF.

    delta_H : float
        Enthalpy change, J/mol.

    temperature_K : float
        Temperature in Kelvin.
    """

    if temperature_K <= 0:

        raise ValueError(
            "Temperature must be greater than zero."
        )

    exponent = (
        -delta_H / R
        * (
            1 / temperature_K
            - 1 / T_REF
        )
    )

    return (
        K_ref
        * math.exp(exponent)
    )


# EQUILIBRIUM CONSTANTS

def get_equilibrium_constants(
    temperature_K
):

    # Calculate temperature-adjusted equilibrium constants.

    K_hydration = temperature_equilibrium_constant(
        K_HYDRATION_REF,
        delta_H=15000,
        temperature_K=temperature_K
    )

    K1 = temperature_equilibrium_constant(
        K1_REF,
        delta_H=1500,
        temperature_K=temperature_K
    )

    K2 = temperature_equilibrium_constant(
        K2_REF,
        delta_H=14000,
        temperature_K=temperature_K
    )

    Kw = temperature_equilibrium_constant(
        KW_REF,
        delta_H=55700,
        temperature_K=temperature_K
    )

    Ka_MEA = temperature_equilibrium_constant(
        MEA_ACIDITY_REF,
        delta_H=41000,
        temperature_K=temperature_K
    )

    K_carbamate = (
        CARBAMATE_K_REF
    )

    return {
        "K_hydration": K_hydration,
        "K1": K1,
        "K2": K2,
        "Kw": Kw,
        "Ka_MEA": Ka_MEA,
        "K_carbamate": K_carbamate
    }


# SPECIES CALCULATION

def calculate_carbonate_species(
    co2_concentration,
    H_concentration,
    temperature_K
):
    """
    Calculate approximate concentrations of:

        CO2
        H2CO3
        HCO3-
        CO3--

    based on equilibrium relationships
    """

    if co2_concentration < 0:

        raise ValueError(
            "CO2 concentration cannot be negative."
        )

    if H_concentration <= 0:

        raise ValueError(
            "Hydrogen-ion concentration must be positive."
        )

    constants = get_equilibrium_constants(
        temperature_K
    )

    K_hydration = constants["K_hydration"]
    K1 = constants["K1"]
    K2 = constants["K2"]


    # CO2 hydration
    # CO2 + H2O <-> H2CO3
    # [H2CO3] = K_hydration [CO2]

    H2CO3 = (
        K_hydration
        * co2_concentration
    )


    # First dissociation
    # H2CO3 <-> H+ + HCO3-
    # K1 = [H+][HCO3-]/[H2CO3]
    # Therefore:
    # [HCO3-] = K1[H2CO3]/[H+]

    HCO3 = (
        K1
        * H2CO3
        / H_concentration
    )


    # Second dissociation
    # HCO3- <-> H+ + CO3--
    # K2 = [H+][CO3--]/[HCO3-]

    CO3 = (
        K2
        * HCO3
        / H_concentration
    )


    return {
        "CO2": co2_concentration,
        "H2CO3": H2CO3,
        "HCO3-": HCO3,
        "CO3--": CO3
    }


# MEA SPECIES CALCULATION

def calculate_mea_species(
    total_MEA,
    H_concentration,
    temperature_K
):
    """
    Calculate free MEA and protonated MEA.

        MEAH+ <-> MEA + H+

    Ka = [MEA][H+]/[MEAH+]

    Therefore:

        [MEAH+] = [MEA][H+]/Ka
    """

    if total_MEA < 0:

        raise ValueError(
            "Total MEA concentration cannot be negative."
        )

    if H_concentration <= 0:

        raise ValueError(
            "Hydrogen-ion concentration must be positive."
        )

    constants = get_equilibrium_constants(
        temperature_K
    )

    Ka_MEA = constants["Ka_MEA"]


    # total_MEA = [MEA] + [MEAH+]
    # [MEAH+] = [MEA][H+]/Ka

    free_MEA = (
        total_MEA
        /
        (
            1
            +
            H_concentration
            / Ka_MEA
        )
    )


    MEAH = (
        total_MEA
        - free_MEA
    )


    return {
        "MEA": free_MEA,
        "MEAH+": MEAH
    }


# MEA CARBAMATE CALCULATION

def calculate_carbamate(
    co2_concentration,
    free_MEA,
    temperature_K
):
    """
    Calculate approximate MEA carbamate concentration.

        CO2 + MEA <-> MEACOO-

    This represents the important carbamate
    formation pathway in aqueous MEA carbon capture.
    """

    if co2_concentration < 0:

        raise ValueError(
            "CO2 concentration cannot be negative."
        )

    if free_MEA < 0:

        raise ValueError(
            "Free MEA concentration cannot be negative."
        )

    constants = get_equilibrium_constants(
        temperature_K
    )

    K_carbamate = (
        constants["K_carbamate"]
    )


    carbamate = (
        K_carbamate
        * co2_concentration
        * free_MEA
    )


    return carbamate


# HYDROXIDE CONCENTRATION

def hydroxide_concentration(
    H_concentration,
    temperature_K
):
    """
    Calculate OH- concentration using:

        Kw = [H+][OH-]
    """

    if H_concentration <= 0:

        raise ValueError(
            "Hydrogen-ion concentration must be positive."
        )

    constants = get_equilibrium_constants(
        temperature_K
    )

    Kw = constants["Kw"]


    return (
        Kw
        / H_concentration
    )


# SPECIES SUMMARY

def equilibrium_state(
    total_co2,
    total_MEA,
    H_concentration,
    temperature_K
):
    """
    Calculate the chemical species present in the
    MEA-CO2-H2O system.

    Parameters
    ----------
    total_co2 : float
        Total dissolved CO2 concentration, mol/L.

    total_MEA : float
        Total MEA concentration, mol/L.

    H_concentration : float
        Hydrogen-ion concentration, mol/L.

    temperature_K : float
        Temperature, K.
    """

    if total_co2 < 0:

        raise ValueError(
            "Total CO2 concentration cannot be negative."
        )

    if total_MEA <= 0:

        raise ValueError(
            "Total MEA concentration must be greater than zero."
        )

    if H_concentration <= 0:

        raise ValueError(
            "Hydrogen-ion concentration must be positive."
        )


    carbonate = calculate_carbonate_species(
        total_co2,
        H_concentration,
        temperature_K
    )


    mea = calculate_mea_species(
        total_MEA,
        H_concentration,
        temperature_K
    )


    carbamate = calculate_carbamate(
        carbonate["CO2"],
        mea["MEA"],
        temperature_K
    )


    OH = hydroxide_concentration(
        H_concentration,
        temperature_K
    )


    return {
        **carbonate,
        **mea,
        "MEACOO-": carbamate,
        "H+": H_concentration,
        "OH-": OH
    }