from absorber import simulate_absorber
from regen import simulate_regeneration

DEFAULT_MAX_ITERATIONS = 50
DEFAULT_CONVERGENCE_TOLERANCE = 0.0001
DEFAULT_CO2_BALANCE_TOLERANCE = 0.05

def validate_inputs(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    temperature_k,
    pressure_Pa,
    column_height,
    column_area,
    mass_coeff,
    initial_loading,
    regeneration_temp,
    feed_temp
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
        raise ValueError("Absorber temperature must be above 0°C (273.15 K).")
    if pressure_Pa <= 0:
        raise ValueError("Pressure must be positive.")
    if column_height <= 0:
        raise ValueError("Column height must be positive.")
    if column_area <= 0:
        raise ValueError("Column area must be positive.")
    if mass_coeff <= 0:
        raise ValueError("Mass-transfer coefficient must be positive.")
    if initial_loading < 0:
        raise ValueError("Initial loading cannot be negative.")
    if regeneration_temp <= 273.15:
        raise ValueError("Regeneration temperature must be above 0°C (273.15 K).")
    if feed_temp <= 273.15:
        raise ValueError("Feed temperature must be above 0°C (273.15 K).")
    if feed_temp >= regeneration_temp:
        raise ValueError("Regeneration temperature must be strictly greater than feed temperature.")

def calculate_relative_error(value_1, value_2):
    denominator = max(abs(value_1), abs(value_2), 1e-12)
    return abs(value_1 - value_2) / denominator

def simulate_process(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    absorber_temp,
    pressure_Pa,
    column_height,
    column_area,
    mass_coeff,
    regeneration_temp,
    feed_temp,
    initial_loading,
    cp,
    heat_desorption,
    stripping_water_per_co2,
    reboiler_efficiency,
    regeneration_effectiveness,
    damping_factor,
    iterations=DEFAULT_MAX_ITERATIONS,
    convergence_tolerance=DEFAULT_CONVERGENCE_TOLERANCE,
    co2_balance_tolerance=DEFAULT_CO2_BALANCE_TOLERANCE
):

    validate_inputs(
        co2_inlet_fraction,
        gas_flow,
        solvent_flow,
        mea_mass_fraction,
        absorber_temp,
        pressure_Pa,
        column_height,
        column_area,
        mass_coeff,
        initial_loading,
        regeneration_temp,
        feed_temp
    )

    if iterations <= 0:
        raise ValueError("Maximum iterations must be positive.")
    if convergence_tolerance <= 0:
        raise ValueError("Convergence tolerance must be positive.")
    if co2_balance_tolerance <= 0:
        raise ValueError("CO2 balance tolerance must be positive.")
    if not 0 < regeneration_effectiveness <= 1:
        raise ValueError("Regeneration effectiveness must be between 0 and 1.")
    if not 0 < damping_factor <= 1:
        raise ValueError("Damping factor must be greater than 0 and less than or equal to 1.")

    converged = False
    iteration_history = []
    completed_iterations = 0
    loading_diff = float("inf")
    lean_loading = initial_loading

    for iteration in range(1, iterations + 1):
        completed_iterations = iteration
        ab_results = simulate_absorber(
            co2_inlet_fraction,
            gas_flow,
            solvent_flow,
            mea_mass_fraction,
            absorber_temp,
            pressure_Pa,
            column_height,
            column_area,
            mass_coeff,
            lean_loading
        )
        
        rich_loading = ab_results["rich_loading"]
        absorber_co2_absorbed = max(ab_results["co2_absorbed_mol_s"], 0.0)

        regeneration_results = simulate_regeneration(
            rich_loading,
            lean_loading,
            solvent_flow,
            mea_mass_fraction,
            regeneration_temp,
            feed_temp,
            cp,
            heat_desorption,
            stripping_water_per_co2,
            reboiler_efficiency
        )
        
        mea_molar_flow = regeneration_results["MEA_molar_flow_mol_s"]
        loading_change = max(rich_loading - lean_loading, 0.0)
        physical_co2_desorption = mea_molar_flow * loading_change
        req_co2_desorption = max(regeneration_results["CO2_desorbed_mol_s"], 0.0)
        actual_co2_desorption = min(req_co2_desorption, physical_co2_desorption)

        if mea_molar_flow > 0:
            regenerated_loading = actual_co2_desorption / mea_molar_flow
        else:
            regenerated_loading = 0.0

        effective_loading_change = regeneration_effectiveness * regenerated_loading
        calc_lean_loading = max(rich_loading - effective_loading_change, 0.0)
        next_lean_loading = (damping_factor * calc_lean_loading + 
                             (1.0 - damping_factor) * lean_loading)
        loading_diff = abs(next_lean_loading - lean_loading)
        cyclic_co2_transfer = (mea_molar_flow * max(rich_loading - next_lean_loading, 0.0))

        iteration_history.append({
            "iteration": iteration,
            "lean_loading": lean_loading,
            "rich_loading": rich_loading,
            "next_lean_loading": next_lean_loading,
            "co2_absorbed_mol_s": absorber_co2_absorbed,
            "co2_desorbed_mol_s": actual_co2_desorption,
            "cyclic_co2_transfer_mol_s": cyclic_co2_transfer,
            "capture_percentage": ab_results["capture_percentage"],
            "reboiler_duty_kW": regeneration_results["reboiler_duty_kW"],
            "specific_energy_GJ_per_tonne_CO2": regeneration_results["specific_energy_GJ_per_tonne_CO2"],
            "loading_difference": loading_diff
        })

        lean_loading = next_lean_loading
        if loading_diff <= convergence_tolerance:
            converged = True
            break

    if not converged:
        raise RuntimeError(f"Process did not converge within {iterations} iterations. Final tolerance difference: {loading_diff:.6f}")

    final_ab_results = simulate_absorber(
        co2_inlet_fraction,
        gas_flow,
        solvent_flow,
        mea_mass_fraction,
        absorber_temp,
        pressure_Pa,
        column_height,
        column_area,
        mass_coeff,
        lean_loading
    )
    
    final_rich_loading = final_ab_results["rich_loading"]
    final_regeneration_results = simulate_regeneration(
        final_rich_loading,
        lean_loading,
        solvent_flow,
        mea_mass_fraction,
        regeneration_temp,
        feed_temp,
        cp,
        heat_desorption,
        stripping_water_per_co2,
        reboiler_efficiency
    )

    mea_molar_flow = final_regeneration_results["MEA_molar_flow_mol_s"]
    cyclic_co2_transfer = (mea_molar_flow * max(final_rich_loading - lean_loading, 0.0))
    absorber_co2_absorbed = max(final_ab_results["co2_absorbed_mol_s"], 0.0)
    
    # For a converged steady-state cyclic system, these should be approximately equal.
    co2_balance_error = calculate_relative_error(absorber_co2_absorbed, cyclic_co2_transfer)
    if co2_balance_error > co2_balance_tolerance:
        raise RuntimeError(f"CO2 mass balance failed. Absorbed = {absorber_co2_absorbed:.6f} mol/s, Cyclic transfer = {cyclic_co2_transfer:.6f} mol/s")

    reboiler_duty_kW = final_regeneration_results["reboiler_duty_kW"]
    reboiler_duty_W = reboiler_duty_kW * 1000.0
    co2_mass_flow_kg_s = absorber_co2_absorbed * 0.04401 
    
    if co2_mass_flow_kg_s > 1e-12:
        specific_energy = reboiler_duty_W / co2_mass_flow_kg_s
    else:
        specific_energy = float("inf")

    return {
        "converged": converged,
        "iterations": completed_iterations,
        "lean_loading": lean_loading,
        "rich_loading": final_rich_loading,
        "capture_percentage": final_ab_results["capture_percentage"],
        "outlet_co2_fraction": final_ab_results["outlet_co2_fraction"],
        "co2_absorbed_mol_s": absorber_co2_absorbed,
        "co2_desorbed_mol_s": cyclic_co2_transfer,
        "co2_balance_error": co2_balance_error,
        "co2_balance_error_percentage": co2_balance_error * 100,
        "reboiler_duty_kW": reboiler_duty_kW,
        "specific_energy_J_per_kg": specific_energy,
        "sensible_heat_kW": final_regeneration_results["sensible_heat_kW"],
        "desorption_heat_kW": final_regeneration_results["desorption_heat_kW"],
        "evaporation_heat_kW": final_regeneration_results["evaporation_heat_kW"],
        "iteration_history": iteration_history,
        "absorber_profile": final_ab_results.get("profile", [])}
