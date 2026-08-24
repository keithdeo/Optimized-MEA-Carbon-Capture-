import math

#importing the necessary functions from other modules
from properties import (
    co2_partial_pressure,
    co2_loading,
    water_mass_fraction,
    mea_moles_per_kg_solution
)

from eq import (
    equilibrium_state
)

from kinetics import (
    reaction_rate
)

#trial value for the number of segments in the absorber, more segments = more specificity
DEFAULT_SEGMENTS = 50

def gas_liquid_transfer(
    gas_co2_fraction,
    pressure_Pa,
    liquid_co2_concentration,
    equilibrium_co2_concentration,
    mass_transfer_coefficient,
    interfacial_area
):

    if not 0 <= gas_co2_fraction <= 1:

        raise ValueError(
            "Gas CO2 fraction must be between 0 and 1."
        )
    
     #calculates the tendancy of CO2 to move from the gas phase to the liquid phase
    driving_force = max(
        0.0,
        equilibrium_co2_concentration
        - liquid_co2_concentration
    )

    #calculates the rate of CO2 transfer from the gas phase to the liquid phase
    transfer_rate = (
        mass_transfer_coefficient
        * interfacial_area
        * driving_force
    )

    return transfer_rate

def simulate_absorber(
    co2_inlet_fraction,
    gas_flow,
    solvent_flow,
    mea_mass_fraction,
    temperature_K,
    pressure_Pa,
    column_height,
    column_area,
    mass_transfer_coefficient,
    initial_loading=0.20,
    number_of_segments=DEFAULT_SEGMENTS
):

    #mathemical checks to ensure that the input values are within reasonable ranges
    if not 0 < co2_inlet_fraction < 1:
        raise ValueError(
            "CO2 inlet fraction must be between 0 and 1."
        )

    if gas_flow <= 0:
        raise ValueError(
            "Gas flow must be positive."
        )

    if solvent_flow <= 0:
        raise ValueError(
            "Solvent flow must be positive."
        )

    if not 0 < mea_mass_fraction < 1:
        raise ValueError(
            "MEA mass fraction must be between 0 and 1."
        )

    if temperature_K <= 273.15:
        raise ValueError(
            "Temperature must be above 0°C."
        )

    if pressure_Pa <= 0:
        raise ValueError(
            "Pressure must be positive."
        )

    if column_height <= 0:
        raise ValueError(
            "Column height must be positive."
        )

    if column_area <= 0:
        raise ValueError(
            "Column area must be positive."
        )

    if mass_transfer_coefficient <= 0:
        raise ValueError(
            "Mass-transfer coefficient must be positive."
        )

    if initial_loading < 0:
        raise ValueError(
            "Initial CO2 loading cannot be negative."
        )

    if number_of_segments <= 0:
        raise ValueError(
            "Number of segments must be positive."
        )

    #height per segment and volume of each segment in the absorber column
    dz = column_height / number_of_segments
    segment_volume = column_area * dz

    #variable declarations
    gas_total_flow = gas_flow
    solvent_total_flow = solvent_flow

    gas_co2_fraction = co2_inlet_fraction
    loading = initial_loading

    total_co2_absorbed = 0.0

    profile = []

    # Calculate how many moles of MEA are present
    # in each kg of the MEA/water solution
    mea_moles_per_kg = (
        mea_moles_per_kg_solution(
            mea_mass_fraction
        )
    )
    # Calculate how much of the solution's mass is water.
    water_mass_fraction_value = (
        water_mass_fraction(
            mea_mass_fraction
        )
    )
    # represents how much MEA/water solution is entering the absorber.
    solvent_mass_flow = solvent_total_flow

    # Calculate the molar flow of MEA entering the absorber.
    mea_molar_flow = (
        solvent_mass_flow
        * mea_moles_per_kg
    )

    for segment in range(number_of_segments):

        # Calculate the CO2 partial pressure in the gas at this point in the absorber.

        co2_partial_pressure_value = (
            co2_partial_pressure(
                gas_co2_fraction,
                pressure_Pa
            )
        )

         # Calculate the CO2 loading that the solvent
        # would have at equilibrium under the current
        # temperature, pressure, etc
        equilibrium_loading = (
            equilibrium_state(
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                co2_partial_pressure=
                    co2_partial_pressure_value,
                mea_mass_fraction=
                    mea_mass_fraction,
                loading=loading
            )
        )

        # Calculate the difference between the
        # equilibrium loading and the current loading.
        loading_driving_force = max(
            equilibrium_loading - loading,
            0.0
        )

        # Estimate the current amount of CO2 in
        # the liquid phase within this segment.
        liquid_co2_concentration = (
            loading
            * mea_molar_flow
            / segment_volume
        )

        # Calculate the CO2 concentration the liquid
        # would have if it reached its equilibrium
        # loading.
        equilibrium_co2_concentration = (
            equilibrium_loading
            * mea_molar_flow
            / segment_volume
        )

        # Calculate how quickly CO2 transfers from
        # the gas phase into the liquid phase.
        transfer_rate = (
            gas_liquid_transfer(
                gas_co2_fraction=
                    gas_co2_fraction,

                pressure_Pa=
                    pressure_Pa,

                liquid_co2_concentration=
                    liquid_co2_concentration,

                equilibrium_co2_concentration=
                    equilibrium_co2_concentration,

                mass_transfer_coefficient=
                    mass_transfer_coefficient,

                interfacial_area=
                    segment_volume
            )
        )

        # Calculate how much MEA is still available
        # to react with incoming CO2.
        available_mea = max(
            1.0 - loading,
            0.0
        )

        # Calculate the chemical reaction rate between
        # dissolved CO2 and available MEA.
        reaction = (
            reaction_rate(
                co2_concentration=
                    max(
                        liquid_co2_concentration,
                        0.0
                    ),

                mea_concentration=
                    available_mea
                    * mea_molar_flow
                    / segment_volume,

                temperature_K=
                    temperature_K
            )
        )

        # Convert the mass-transfer rate into the
        # maximum amount of CO2 that can be transferred
        # within this segment.
        transfer_capacity = (
            transfer_rate
            * segment_volume
        )

        # Convert the reaction rate into the maximum
        # amount of CO2 that can react within this
        # segment.
        reaction_capacity = (
            reaction
            * segment_volume
        )

        # Calculate how much CO2 is entering this
        # segment with the gas.
        inlet_co2_flow = (
            gas_co2_fraction
            * gas_total_flow
        )

         # Make sure the amount of available CO2 isnt negative
        maximum_available_co2 = (
            max(
                inlet_co2_flow,
                0.0
            )
        )

        # Determine the effective amount of CO2
        # that can actually be removed.
        effective_reaction_capacity = min(
            transfer_capacity,
            reaction_capacity
        )


        co2_removed = min(
            effective_reaction_capacity,
            maximum_available_co2
        )

        co2_removed = max(
            co2_removed,
            0.0
        )

        # Calculate how much CO2 remains in the gas
        # after passing through this segment.
        outlet_co2_flow = (
            inlet_co2_flow
            - co2_removed
        )

        # Prevent the outlet CO2 flow from becoming
        # negative due to numerical calculations.
        outlet_co2_flow = max(
            outlet_co2_flow,
            0.0
        )

        # Calculate the new CO2 fraction in the gas
        # after this segment.
        gas_co2_fraction = (
            outlet_co2_flow
            / gas_total_flow
        )

        # Calculate how much the solvent's CO2 loading
        # increases because of the CO2 absorbed
        loading_change = (
            co2_removed
            / mea_molar_flow
        )

        # Add the newly absorbed CO2 to the solvent's
        loading += loading_change

        loading = max(
            loading,
            0.0
        )

        # Add the CO2 absorbed in this segment to
        # the total CO2 absorbed by the entire absorber.
        total_co2_absorbed += co2_removed

        # Save the results from this segment.
        profile.append(
            {
                "segment": segment + 1,
                "height_m": (segment + 1) * dz,
                "gas_CO2_fraction": gas_co2_fraction,
                "CO2_partial_pressure_Pa":
                    co2_partial_pressure_value,
                "equilibrium_loading":
                    equilibrium_loading,
                "loading":
                    loading,
                "loading_driving_force":
                    loading_driving_force,
                "CO2_removed_mol_s":
                    co2_removed,
                "reaction_rate":
                    reaction
            }
        )

    inlet_co2_flow = (
        co2_inlet_fraction
        * gas_total_flow
    )

    if inlet_co2_flow > 0:

        capture_percentage = (
            total_co2_absorbed
            / inlet_co2_flow
        ) * 100.0

    else:

        capture_percentage = 0.0

    return {
        "capture_percentage":
            capture_percentage,

        "outlet_co2_fraction":
            gas_co2_fraction,

        "co2_absorbed_mol_s":
            total_co2_absorbed,

        "rich_loading":
            loading,

        "profile":
            profile
    }

if __name__ == "__main__":

    print("          CO2 ABSORBER TEST")


    results = simulate_absorber(

        co2_inlet_fraction=0.12,

        gas_flow=100.0,

        solvent_flow=120.0,

        mea_mass_fraction=0.30,

        temperature_K=313.15,

        pressure_Pa=101325,

        column_height=20.0,

        column_area=10.0,

        mass_transfer_coefficient=0.01,

        initial_loading=0.20,

        number_of_segments=50
    )


    print(
        f"\nCO2 capture: "
        f"{results['capture_percentage']:.2f}%"
    )


    print(
        f"Outlet CO2: "
        f"{results['outlet_co2_fraction'] * 100:.2f}%"
    )


    print(
        f"CO2 absorbed: "
        f"{results['co2_absorbed_mol_s']:.4f} mol/s"
    )


    print(
        f"Rich loading: "
        f"{results['rich_loading']:.4f} mol CO2/mol MEA"
    )


    print("\nFirst 5 absorber segments:")

    for row in results["profile"][:5]:

        print(row)


    print("\n======================================")
