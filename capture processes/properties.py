#Physical properties for an aqueous MEA-CO2-H2O system.
#This file contains:
#Physical constants
#Process conditions
#CO2 partial pressure
#Water/MEA composition
#MEA amount
#CO2 in rich MEA solution

from dataclasses import dataclass


# ============================================================
# 1. PHYSICAL CONSTANTS
# ============================================================

R = 8.314462618
# Universal gas constant
# Units: J/(mol*K)


MW_CO2 = 44.0095
# Molecular weight of CO2
# Units: g/mol


MW_MEA = 61.083
# Molecular weight of monoethanolamine (MEA)
# Units: g/mol


MW_H2O = 18.01528
# Molecular weight of water
# Units: g/mol


# ============================================================
# 2. PROCESS CONDITIONS
# ============================================================

#class set used to hold process conditions
@dataclass
class ProcessConditions:

    temperature_K: float
    # Temperature in Kelvin

    pressure_Pa: float
    # Total pressure in Pascal

    co2_mole_fraction: float
    # CO2 mole fraction in the gas

    mea_mass_fraction: float
    # MEA mass fraction in the liquid

    gas_molar_flow: float
    # Gas flow rate in mol/s

    solvent_molar_flow: float
    # Solvent flow rate in mol/s


# ============================================================
# 3. CO2 PARTIAL PRESSURE
# ============================================================

def co2_partial_pressure(
    co2_mole_fraction,
    pressure_Pa
):

    #calculates the partial pressure of CO2 in the gas phase
    if not 0 <= co2_mole_fraction <= 1:

        raise ValueError(
            "CO2 mole fraction must be between 0 and 1."
        )


    if pressure_Pa <= 0:

        raise ValueError(
            "Pressure must be greater than zero."
        )


    return (
        co2_mole_fraction
        * pressure_Pa
    )


# ============================================================
# 4. WATER MASS FRACTION
# ============================================================

def water_mass_fraction(
    mea_mass_fraction
):
    # calculates mass fraction of water in the MEA/water solution
    if not 0 < mea_mass_fraction < 1:

        raise ValueError(
            "MEA mass fraction must be between 0 and 1."
        )

    return (
        1.0 - mea_mass_fraction
    )


# ============================================================
# 5. SOLVENT MOLECULAR WEIGHT
# ============================================================

def solvent_molecular_weight(
    mea_mass_fraction
):
    # calculates the approximate molecular weight of the MEA/water solution
    water_fraction = (
        water_mass_fraction(
            mea_mass_fraction
        )
    )


    return (
        mea_mass_fraction * MW_MEA
        +
        water_fraction * MW_H2O
    )


# ============================================================
# 6. MEA MOLAR AMOUNT
# ============================================================

def mea_moles_per_kg_solution(
    mea_mass_fraction
):
    #calculates the number of moles of MEA in 1 kg of solution (approximate)

    if not 0 < mea_mass_fraction < 1:

        raise ValueError(
            "MEA mass fraction must be between 0 and 1."
        )


    mass_MEA_kg = (
        mea_mass_fraction
    )


    mass_MEA_g = (
        mass_MEA_kg * 1000
    )


    moles_MEA = (
        mass_MEA_g
        / MW_MEA
    )


    return moles_MEA


# ============================================================
# 7. CO2 LOADING
# ============================================================

def co2_loading(
    co2_moles,
    mea_moles
):
    # calculates the ratio in mol CO2/mol MEA
    #  high ratios indicate high loading of CO2 in the solution

    if co2_moles < 0:

        raise ValueError(
            "CO2 amount cannot be negative."
        )


    if mea_moles <= 0:

        raise ValueError(
            "MEA amount must be greater than zero."
        )


    return (
        co2_moles
        / mea_moles
    )
