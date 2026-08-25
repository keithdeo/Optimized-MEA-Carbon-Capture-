import random
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "Math"))

from process import simulate_process

NUMBER_OF_SAMPLES = 7500
results = []

for i in range(NUMBER_OF_SAMPLES):

    co2_inlet = random.uniform(0.05, 0.20)
    gas_flow = random.uniform(80.0, 150.0)
    solvent_flow = random.uniform(80.0, 180.0)
    mea_mass_fraction = random.uniform(0.20, 0.40)
    absorber_temperature_K = random.uniform(303.15,323.15)
    pressure_Pa = random.uniform(101325,200000)
    column_height = random.uniform(10.0,30.0)
    column_area = random.uniform(5.0,15.0)
    mass_transfer_coefficient = random.uniform(0.005,0.02)
    regeneration_temperature_K = random.uniform(383.15,403.15)
    feed_temperature_K = random.uniform(303.15,323.15)
    initial_loading = random.uniform(0.15,0.25)

    try:
        process_results = simulate_process(
            co2_inlet_fraction=co2_inlet,
            gas_flow=gas_flow,
            solvent_flow=solvent_flow,
            mea_mass_fraction=mea_mass_fraction,
            absorber_temperature_K=absorber_temperature_K,
            pressure_Pa=pressure_Pa,
            column_height=column_height,
            column_area=column_area,
            mass_transfer_coefficient=mass_transfer_coefficient,
            regeneration_temperature_K=regeneration_temperature_K,
            feed_temperature_K=feed_temperature_K,
            initial_loading=initial_loading
        )

        results.append({
            "CO2_inlet":
                co2_inlet,
            "gas_flow_kg_s":
                gas_flow,
            "solvent_flow_kg_s":
                solvent_flow,
            "MEA_mass_fraction":
                mea_mass_fraction,
            "absorber_temperature_K":
                absorber_temperature_K,
            "pressure_Pa":
                pressure_Pa,
            "column_height_m":
                column_height,
            "column_area_m2":
                column_area,
            "mass_transfer_coefficient":
                mass_transfer_coefficient,
            "regeneration_temperature_K":
                regeneration_temperature_K,
            "feed_temperature_K":
                feed_temperature_K,
            "initial_loading":
                initial_loading,
            "capture_percentage":
                process_results["capture_percentage"],
            "CO2_absorbed_mol_s":
                process_results["co2_absorbed_mol_s"],
            "CO2_desorbed_mol_s":
                process_results["co2_desorbed_mol_s"],
            "rich_loading":
                process_results["rich_loading"],
            "lean_loading":
                process_results["lean_loading"],
            "reboiler_duty_kW":
                process_results["reboiler_duty_kW"],
            "specific_energy_GJ_per_tonne_CO2":
                process_results[
                    "specific_energy_GJ_per_tonne_CO2"
                ]
        })

    except (ValueError, RuntimeError):
        continue

data = pd.DataFrame(results)

data.to_csv(
    "carbon_capture_dataset.csv",
    index=False
)

print("DATASET GENERATION COMPLETE")
print(f"Requested samples: {NUMBER_OF_SAMPLES}")
print(f"Successful simulations: {len(data)}")
print(f"Failed simulations: {NUMBER_OF_SAMPLES - len(data)}")
