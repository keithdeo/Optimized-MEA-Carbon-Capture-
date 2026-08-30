from dataclasses import dataclass

R = 8.314462618
# Universal gas constant

MW_CO2 = 44.0095
# Molecular weight of CO2

MW_MEA = 61.083
# Molecular weight of MEA

MW_H2O = 18.01528
# Molecular weight of water

@dataclass
class ProcessConditions:
    temperature: float #kelvin
    pressure: float #pascal
    co2_mole_fraction: float
    mea_mass_fraction: float
    gas_molar_flow: float #mol/s
    solvent_molar_flow: float #mol/s

def co2_partial_pressure(co2_mole_fraction, pressure):

    if not 0 <= co2_mole_fraction <= 1 or not pressure > 0:
        raise ValueError("CO2 mole fraction must be between 0 and 1. and pressure must be positive.")

    return (co2_mole_fraction * pressure)

def water_mass_fraction(mea_mass_fraction):
    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")

    return (1.0 - mea_mass_fraction)

def solvent_mw(mea_mass_fraction):
    water_fraction = water_mass_fraction(mea_mass_fraction)

    return 1.0 / ((mea_mass_fraction / MW_MEA) + (water_fraction / MW_H2O))

def solution_mass(mea_mass_fraction, solution_volume, solution_density):
    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")
    if solution_volume <= 0:
        raise ValueError("Solution volume must be positive.")
    if solution_density <= 0:
        raise ValueError("Solution density must be positive.")

    return (solution_volume * solution_density)

def mea_moles(mea_mass_fraction, solution_mass):
    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")
    if solution_mass <= 0:
        raise ValueError("Solution mass must be positive.")

    mea_mass = mea_mass_fraction * solution_mass

    return (mea_mass / MW_MEA)

def co2_moles(co2_absorbed_mass):
    if co2_absorbed_mass < 0:
        raise ValueError("CO2 absorbed mass cannot be negative.")

    return (co2_absorbed_mass / MW_CO2)

def co2_loading(co2_absorbed_moles, mea_mass_fraction, solution_mass):
    # calculates the ratio in mol CO2/mol MEA
    # high ratios indicate high loading of CO2 in the solution

    if co2_absorbed_moles < 0:
        raise ValueError("CO2 amount cannot be negative.")

    n_mea = mea_moles(mea_mass_fraction, solution_mass)

    if n_mea <= 0:
        raise ValueError("MEA amount must be greater than zero.")

    return (co2_absorbed_moles / n_mea)
