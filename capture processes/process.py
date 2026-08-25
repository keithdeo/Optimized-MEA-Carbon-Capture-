from absorber import simulate_absorber
from regen import simulate_regeneration


# Set system-wide default convergence parameters
DEFAULT_MAX_ITERATIONS = 50
DEFAULT_CONVERGENCE_TOLERANCE = 0.0001

def validate_inputs(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    temperature_K,
    pressure_Pa,
    column_height,
    column_area,
    mass_transfer_coefficient,
    initial_loading,
    regeneration_temperature_K,
    feed_temperature_K
):
    # Verify that all values are positive
    if not 0 < co2_inlet_fraction < 1:
        raise ValueError("CO2 inlet fraction must be between 0 and 1.")

    if gas_flow <= 0:
        raise ValueError("Gas flow must be positive.")

    if solvent_flow <= 0:
        raise ValueError("Solvent flow must be positive.")

    if not 0 < mea_mass_fraction < 1:
        raise ValueError("MEA mass fraction must be between 0 and 1.")

    if temperature_K <= 273.15:
        raise ValueError("Absorber temperature must be above 0°C (273.15 K).")

    if pressure_Pa <= 0:
        raise ValueError("Pressure must be positive.")

    if column_height <= 0:
        raise ValueError("Column height must be positive.")

    if column_area <= 0:
        raise ValueError("Column area must be positive.")

    if mass_transfer_coefficient <= 0:
        raise ValueError("Mass-transfer coefficient must be positive.")

    if initial_loading < 0:
        raise ValueError("Initial loading cannot be negative.")

    if regeneration_temperature_K <= 273.15:
        raise ValueError("Regeneration temperature must be above 0°C (273.15 K).")

    if feed_temperature_K <= 273.15:
        raise ValueError("Feed temperature must be above 0°C (273.15 K).")

    if feed_temperature_K >= regeneration_temperature_K:
        raise ValueError("Regeneration temperature must be strictly greater than feed temperature.")


#process function to simulate the entire CO2 capture process
def simulate_process(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    absorber_temperature_K,
    pressure_Pa,
    column_height,
    column_area,
    mass_transfer_coefficient,
    regeneration_temperature_K,
    feed_temperature_K,
    initial_loading=0.20,
    cp_solution_kJ_kgK=4.0,
    heat_of_desorption_J_mol=84000.0,
    stripping_ratio_mol_water_per_mol_co2=1.2,
    reboiler_efficiency=0.85,
    regeneration_effectiveness=1.0,
    damping_factor=0.5,
    max_iterations=DEFAULT_MAX_ITERATIONS,
    convergence_tolerance=DEFAULT_CONVERGENCE_TOLERANCE
):
    # Validate physical bounds of all simulation input variables
    validate_inputs(
        co2_inlet_fraction=co2_inlet_fraction,
        gas_flow=gas_flow,
        solvent_flow=solvent_flow,
        mea_mass_fraction=mea_mass_fraction,
        temperature_K=absorber_temperature_K,
        pressure_Pa=pressure_Pa,
        column_height=column_height,
        column_area=column_area,
        mass_transfer_coefficient=mass_transfer_coefficient,
        initial_loading=initial_loading,
        regeneration_temperature_K=regeneration_temperature_K,
        feed_temperature_K=feed_temperature_K
    )

    # Check that loop boundary condition limits are non-zero positive numbers
    if max_iterations <= 0:
        raise ValueError("Maximum iterations must be positive.")

    # Ensure convergence tolerance threshold is non-zero positive number
    if convergence_tolerance <= 0:
        raise ValueError("Convergence tolerance must be positive.")

    if not 0 < regeneration_effectiveness <= 1:
        raise ValueError(
            "Regeneration effectiveness must be between 0 and 1."
        )

    if not 0 < damping_factor <= 1:
        raise ValueError(
            "Damping factor must be greater than 0 and less than or equal to 1."
        )

    # Initialize recycled lean solvent CO2 loading from starting estimate
    lean_loading = initial_loading

    # Set loop tracker flags for convergence state
    converged = False

    # Initialize empty list to record convergence steps
    iteration_history = []

    # Initialize counter to track executed iteration cycles
    completed_iterations = 0

    # Execute tear-stream recycling numerical loop until convergence
    for iteration in range(1, max_iterations + 1):

        # Update completed iteration counter step
        completed_iterations = iteration

        # Simulate gas absorption column using current lean loading
        absorber_results = simulate_absorber(
            co2_inlet_fraction=co2_inlet_fraction,
            gas_flow=gas_flow,
            solvent_flow=solvent_flow,
            mea_mass_fraction=mea_mass_fraction,
            temperature_K=absorber_temperature_K,
            pressure_Pa=pressure_Pa,
            column_height=column_height,
            column_area=column_area,
            mass_transfer_coefficient=mass_transfer_coefficient,
            initial_loading=lean_loading
        )

        # Extract resulting rich solvent CO2 loading leaving absorber bottom
        rich_loading = absorber_results["rich_loading"]

        # Simulate stripper reboiler column energy performance and CO2 desorption
        regeneration_results = simulate_regeneration(
            rich_loading=rich_loading,
            lean_loading=lean_loading,
            solvent_flow_kg_s=solvent_flow,
            mea_mass_fraction=mea_mass_fraction,
            regeneration_temperature_K=regeneration_temperature_K,
            feed_temperature_K=feed_temperature_K,
            cp_solution_kJ_kgK=cp_solution_kJ_kgK,
            heat_of_desorption_J_mol=heat_of_desorption_J_mol,
            stripping_ratio_mol_water_per_mol_co2=stripping_ratio_mol_water_per_mol_co2,
            reboiler_efficiency=reboiler_efficiency
        )

        # Calculate the total CO2 loading removed by regeneration
        regenerated_loading_change = (
            regeneration_results["CO2_desorbed_mol_s"]
            / regeneration_results["MEA_molar_flow_mol_s"]
        )

        # Limit regeneration to the physically available loading difference
        maximum_loading_change = max(
            rich_loading - lean_loading,
            0.0
        )

        regenerated_loading_change = min(
            regenerated_loading_change,
            maximum_loading_change
        )

        # Apply regeneration effectiveness
        effective_loading_change = (
            regeneration_effectiveness
            * regenerated_loading_change
        )

        # Calculate the new lean loading leaving the regenerator
        calculated_lean_loading = max(
            rich_loading - effective_loading_change,
            0.0
        )

        # Apply numerical relaxation damping to avoid loop instability
        next_lean_loading = (
            damping_factor * calculated_lean_loading
            + (1.0 - damping_factor) * lean_loading
        )

        # Calculate absolute change in lean solvent loading for iteration
        loading_difference = abs(
            next_lean_loading - lean_loading
        )

        # Record step output parameters in historical tracking dictionary
        iteration_history.append({
            "iteration": iteration,
            "lean_loading": lean_loading,
            "rich_loading": rich_loading,
            "next_lean_loading": next_lean_loading,
            "capture_percentage": absorber_results["capture_percentage"],
            "reboiler_duty_kW": regeneration_results["reboiler_duty_kW"],
            "specific_energy_GJ_per_tonne_CO2":
                regeneration_results["specific_energy_GJ_per_tonne_CO2"],
            "loading_difference": loading_difference
        })

        # Update recycled lean solvent loading for the next iteration
        lean_loading = next_lean_loading

        # Check if change in solvent state falls within target convergence criteria
        if loading_difference <= convergence_tolerance:
            converged = True
            break

    # Raise runtime error if recycling loop fails to stabilize within max steps
    if not converged:
        raise RuntimeError(
            f"Process did not converge within {max_iterations} iterations. "
            f"Final tolerance diff: {loading_difference:.6f}"
        )

    # Perform final verification pass of absorber column with converged solvent loop
    final_absorber_results = simulate_absorber(
        co2_inlet_fraction=co2_inlet_fraction,
        gas_flow=gas_flow,
        solvent_flow=solvent_flow,
        mea_mass_fraction=mea_mass_fraction,
        temperature_K=absorber_temperature_K,
        pressure_Pa=pressure_Pa,
        column_height=column_height,
        column_area=column_area,
        mass_transfer_coefficient=mass_transfer_coefficient,
        initial_loading=lean_loading
    )

    # Extract final rich loading under fully converged conditions
    final_rich_loading = final_absorber_results["rich_loading"]

    # Perform final verification pass of regeneration column with converged solvent loop
    final_regeneration_results = simulate_regeneration(
        rich_loading=final_rich_loading,
        lean_loading=lean_loading,
        solvent_flow_kg_s=solvent_flow,
        mea_mass_fraction=mea_mass_fraction,
        regeneration_temperature_K=regeneration_temperature_K,
        feed_temperature_K=feed_temperature_K,
        cp_solution_kJ_kgK=cp_solution_kJ_kgK,
        heat_of_desorption_J_mol=heat_of_desorption_J_mol,
        stripping_ratio_mol_water_per_mol_co2=stripping_ratio_mol_water_per_mol_co2,
        reboiler_efficiency=reboiler_efficiency
    )

    # Return key-value dictionary containing all final converged plant performance metrics
    return {
        "converged": converged,
        "iterations": completed_iterations,
        "lean_loading": lean_loading,
        "rich_loading": final_rich_loading,
        "capture_percentage": final_absorber_results["capture_percentage"],
        "outlet_co2_fraction": final_absorber_results["outlet_co2_fraction"],
        "co2_absorbed_mol_s": final_absorber_results["co2_absorbed_mol_s"],
        "co2_desorbed_mol_s": final_regeneration_results["CO2_desorbed_mol_s"],
        "reboiler_duty_kW": final_regeneration_results["reboiler_duty_kW"],
        "specific_energy_GJ_per_tonne_CO2":
            final_regeneration_results["specific_energy_GJ_per_tonne_CO2"],
        "sensible_heat_kW": final_regeneration_results["sensible_heat_kW"],
        "desorption_heat_kW": final_regeneration_results["desorption_heat_kW"],
        "evaporation_heat_kW": final_regeneration_results["evaporation_heat_kW"],
        "iteration_history": iteration_history,
        "absorber_profile": final_absorber_results.get("profile", [])
    }
