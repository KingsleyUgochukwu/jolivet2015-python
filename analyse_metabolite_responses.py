
import numpy as np
import pandas as pd

durations = [300, 600, 1000, 1500, 2000]

results = []

for duration in durations:

    filename = f"ischemia_reperfusion_{duration}s.csv"
    data = pd.read_csv(filename)

    # Define the important time points
    ischemia_onset = 100.0
    reperfusion_onset = ischemia_onset + duration

    # Baseline: first recorded point at ischemia onset
    baseline = data.iloc[0]

    # End of ischemia: closest recorded point to reperfusion onset
    end_index = (data["time_s"] - reperfusion_onset).abs().idxmin()
    end_ischemia = data.loc[end_index]

    # Reperfusion endpoint
    end_reperfusion = data.iloc[-1]

    # Data recorded during ischemia
    ischemia_data = data[
        (data["time_s"] >= ischemia_onset) &
        (data["time_s"] <= reperfusion_onset)
    ]

    baseline_glucose = baseline["extracellular_glucose"]
    baseline_lactate = baseline["extracellular_lactate"]

    minimum_glucose = ischemia_data["extracellular_glucose"].min()
    maximum_lactate = ischemia_data["extracellular_lactate"].max()

    min_glucose_row = ischemia_data.loc[
        ischemia_data["extracellular_glucose"].idxmin()
    ]

    max_lactate_row = ischemia_data.loc[
        ischemia_data["extracellular_lactate"].idxmax()
    ]

    results.append({
        "ischemia_duration_s": duration,

        "glucose_end_ischemia_percent_baseline":
            end_ischemia["extracellular_glucose"] /
            baseline_glucose * 100,

        "glucose_minimum_percent_baseline":
            minimum_glucose / baseline_glucose * 100,

        "glucose_end_reperfusion_percent_baseline":
            end_reperfusion["extracellular_glucose"] /
            baseline_glucose * 100,

        "lactate_end_ischemia_percent_baseline":
            end_ischemia["extracellular_lactate"] /
            baseline_lactate * 100,

        "lactate_maximum_percent_baseline":
            maximum_lactate / baseline_lactate * 100,

        "lactate_end_reperfusion_percent_baseline":
            end_reperfusion["extracellular_lactate"] /
            baseline_lactate * 100,

        "time_of_minimum_glucose_s":
            min_glucose_row["time_s"] - ischemia_onset,

        "time_of_maximum_lactate_s":
            max_lactate_row["time_s"] - ischemia_onset
    })

summary = pd.DataFrame(results)

print("\nGLUCOSE AND LACTATE RESPONSE SUMMARY")
print("=" * 90)
print(summary.round(2).to_string(index=False))

summary.to_csv(
    "metabolite_response_summary.csv",
    index=False
)

print("\nSaved: metabolite_response_summary.csv")
