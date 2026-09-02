import math

from eq import eq_loading
from kinetics import reaction_rate
from properties import (
    co2_partial_pressure, 
    mea_moles, 
    water_mass_fraction,
    henrys_constant,
    co2_diffusivity
)

DEFAULT_SEGMENTS = 100  # Increased for better mathematical resolution
DEFAULT_SOLVENT_DENSITY = 1050.0
DEFAULT_CO2_BALANCE_TOLERANCE = 0.05

def calculate_enhancement_factor(reaction_rate_vol, free_co2_conc, mass_coeff, diffusivity):
    if free_co2_conc <= 0:
        return 1.0
        
    k_app = reaction_rate_vol / free_co2_conc 
    hatta_sq = (k_app * diffusivity) / (mass_coeff ** 2)
    hatta = math.sqrt(max(hatta_sq, 0.0))
    return math.sqrt(1.0 + hatta_sq)

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

    if not 0 < co2_inlet_fraction < 1: raise ValueError("CO2 inlet fraction must be between 0 and 1.")
    if gas_flow <= 0: raise ValueError("Gas flow must be positive.")
    if solvent_flow <= 0: raise ValueError("Solvent flow must be positive.")
    if not 0 < mea_mass_fraction < 1: raise ValueError("MEA mass fraction must be between 0 and 1.")
    if temperature_k <= 273.15: raise ValueError("Temperature must be above 0°C.")
    if pressure_Pa <= 0: raise ValueError("Pressure must be positive.")
    if column_height <= 0: raise ValueError("Column height must be positive.")
    if column_area <= 0: raise ValueError("Column area must be positive.")
    if mass_transfer_coefficient <= 0: raise ValueError("Mass-transfer coefficient must be positive.")
    if specific_interfacial <= 0: raise ValueError("Specific interfacial area must be positive.")
    if initial_loading < 0: raise ValueError("Initial CO2 loading cannot be negative.")
    if number_of_segments <= 0: raise ValueError("Number of segments must be positive.")
    if solvent_density_kg_m3 <= 0: raise ValueError("Solvent density must be positive.")
    if co2_balance_tolerance <= 0: raise ValueError("CO2 balance tolerance must be positive.")

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
    
    H_co2 = henrys_constant(temperature_k, mea_mass_fraction)
    D_co2 = co2_diffusivity(temperature_k, mea_mass_fraction)

    def evaluate_column(bottom_loading):

        current_gas_co2 = inlet_co2_flow
        current_liquid_loading = bottom_loading
        total_absorbed = 0.0

        profile = []
        for segment in range(number_of_segments):

            total_gas = inert_gas_flow + current_gas_co2
            gas_co2_fraction = current_gas_co2 / max(total_gas, 1e-12)
            Pp_co2 = co2_partial_pressure(gas_co2_fraction, pressure_Pa)
            free_co2_concentration = Pp_co2 / H_co2

            eq_loading_val = eq_loading(
                temperature_k,
                pressure_Pa,
                Pp_co2,
                mea_mass_fraction,
                current_liquid_loading,
            )

            loading_driving_force = max(eq_loading_val - current_liquid_loading, 0.0)
            available_mea_fraction = max(1.0 - 2.0 * current_liquid_loading, 0.0)
            available_mea_concentration = available_mea_fraction * mea_concentration

            reaction = reaction_rate(
                temperature_k, 
                max(free_co2_concentration, 0.0), 
                max(available_mea_concentration, 0.0)
            )

            E = calculate_enhancement_factor(
                reaction_rate_vol=reaction, 
                free_co2_conc=free_co2_concentration, 
                mass_coeff=mass_transfer_coefficient, 
                diffusivity=D_co2
            )
            
            effective_transfer_rate = mass_transfer_coefficient * E * total_interfacial * free_co2_concentration 
            
            # --- STIFF ODE FIX APPLIED HERE ---
            # By capping the transfer to 25% of the available gap per segment, 
            # we force the solver to take smooth steps toward equilibrium instead of violently crashing.
            max_loading_capacity = (loading_driving_force * mea_molar_flow) * 0.25
            max_gas_capacity = current_gas_co2 * 0.25
            
            co2_removed = min(effective_transfer_rate, max_gas_capacity, max_loading_capacity)
            co2_removed = max(co2_removed, 0.0)
            # ----------------------------------

            profile.append({
                "segment": segment + 1,
                "height_m": (segment + 1) * dz,
                "gas_CO2_fraction": gas_co2_fraction,
                "co2_partial_pressure": Pp_co2,
                "equilibrium_loading": eq_loading_val,
                "loading": current_liquid_loading,
                "free_co2_conc": free_co2_concentration,
                "enhancement_factor": E,
                "co2_removed_mol_s": co2_removed
            })

            current_gas_co2 -= co2_removed
            current_gas_co2 = max(current_gas_co2, 0.0)
            loading_change = co2_removed / max(mea_molar_flow, 1e-12)
            current_liquid_loading -= loading_change
            current_liquid_loading = max(current_liquid_loading, 0.0)
            
            total_absorbed += co2_removed

        error = current_liquid_loading - initial_loading
        return (error, profile, current_gas_co2, current_liquid_loading)

    low_guess = initial_loading  
    high_guess = 1.0 

    for _ in range(50):
        mid_guess = (low_guess + high_guess) / 2.0
        error, _, _, _ = evaluate_column(mid_guess)
        if abs(error) < 1e-5:
            break
        if error > 0:
            high_guess = mid_guess
        else:
            low_guess = mid_guess

    (final_error, final_profile, final_gas_co2, final_top_loading) = evaluate_column(mid_guess)
    total_co2_absorbed_gas_side = inlet_co2_flow - final_gas_co2

    co2_in = inlet_co2_flow + (mea_molar_flow * initial_loading)
    co2_out = final_gas_co2 + (mea_molar_flow * mid_guess)
    co2_balance_error = abs(co2_in - co2_out) / max(co2_in, 1e-12)

    if inlet_co2_flow > 0:
        capture_percentage = (total_co2_absorbed_gas_side / inlet_co2_flow) * 100.0
    else:
        capture_percentage = 0.0
    
    final_total_gas = inert_gas_flow + final_gas_co2
    outlet_co2_fraction = final_gas_co2 / max(final_total_gas, 1e-12)
    
    return {
        "capture_percentage": capture_percentage,
        "outlet_co2_fraction": outlet_co2_fraction,
        "co2_absorbed_mol_s": total_co2_absorbed_gas_side,
        "co2_balance_error": co2_balance_error,
        "co2_balance_error_percentage": co2_balance_error * 100.0,
        "absorber_balance_passed": co2_balance_error <= co2_balance_tolerance,
        "rich_loading": mid_guess,
        "mea_molar_flow": mea_molar_flow,
        "water_mass_fraction": water_mass_fraction_value,
        "profile": final_profile
    }