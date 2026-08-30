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

if __name__ == "__main__":

    print("PROPERTY TEST")

    # Example process conditions

    temperature = float(input("Enter temperature in Kelvin: "))

    pressure = float(input("Enter pressure in Pascals: "))

    co2_fraction = float(input("Enter mole fraction of CO2: "))

    mea_fraction = float(input("Enter mass fraction of MEA: "))

    volume = float(input("Enter volume of solution in mL: "))

    density = float(input("Enter density of solution in g/mL: "))

    co2_absorbed_mass = float(input("Enter absorbed CO2 mass in g: "))

    Pp_CO2 = co2_partial_pressure(co2_fraction, pressure)
    print(f"\nCO2 partial pressure: {Pp_CO2:.2f} Pa")

    water_fraction = water_mass_fraction(mea_fraction)
    print(f"Water mass fraction: {water_fraction:.2f}")

    MW_solution = solvent_mw(mea_fraction)
    print(f"Approximate solvent MW: {MW_solution:.2f} g/mol")

    calc_solution_mass = solution_mass(mea_fraction, volume, density)

    calc_mea_moles = mea_moles(mea_fraction, calc_solution_mass)

    mea_in_kg_sln = mea_fraction * 1.0

    print(f"MEA in moles: {calc_mea_moles:.2f} mol")
    print(f"MEA mass in 1 kg of solution {mea_in_kg_sln}")

    co2_amount = co2_moles(co2_absorbed_mass)

    calc_co2_loading = co2_loading(co2_amount, mea_fraction, calc_solution_mass)

    print(f"Example CO2 loading: {calc_co2_loading:.3f} mol CO2/mol MEA")
