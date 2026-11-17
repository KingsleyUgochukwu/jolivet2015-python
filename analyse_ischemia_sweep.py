
import numpy as np
import pandas as pd
from pathlib import Path

# Ischemia durations tested (seconds)
durations = [300, 600, 1000, 1500, 2000]

# Ischemia begins at 100 seconds.
ischemia_onset = 100.0
reperfusion_duration = 1000.0

results = []

for duration in durations:

    filename = f"ischemia_reperfusion_{duration}s.csv"

    if not Path(filename).exists():
        raise FileNotFoundError(f"Missing dataset: {filename}")

    data = pd.read_csv(filename)

    # Select the recorded row closest to a target time.
    def get_nearest_row(target_time):
        index = (data["time_s"] - target_time).abs().idxmin()
        row = data.loc[index]

        difference = abs(row["time_s"] - target_time)

        if difference > 1e-5:
            raise ValueError(
                f"Could not locate target time {target_time} "
                f"in {filename}. Closest difference: {difference}"
            )

        return row

    baseline = get_nearest_row(ischemia_onset)

    end_ischemia_time = ischemia_onset + duration
    end_ischemia = get_nearest_row(end_ischemia_time)

    end_reperfusion = data.iloc[-1]

    # ATP measurements
    n_base = baseline["neuronal_ATP"]
    a_base = baseline["astrocytic_ATP"]

    n_ischemia = end_ischemia["neuronal_ATP"]
    a_ischemia = end_ischemia["astrocytic_ATP"]

    n_reperfusion = end_reperfusion["neuronal_ATP"]
    a_reperfusion = end_reperfusion["astrocytic_ATP"]

    # Percentage ATP depletion
    n_depletion = (n_base - n_ischemia) / n_base * 100
    a_depletion = (a_base - a_ischemia) / a_base * 100

    # Percentage of initial ATP depletion recovered
    n_recovery = (
        (n_reperfusion - n_ischemia) /
        (n_base - n_ischemia) * 100
    )

    a_recovery = (
        (a_reperfusion - a_ischemia) /
        (a_base - a_ischemia) * 100
    )

    results.append({
        "ischemia_duration_s": duration,

        "neuronal_ATP_baseline": n_base,
        "neuronal_ATP_end_ischemia": n_ischemia,
        "neuronal_ATP_depletion_percent": n_depletion,
        "neuronal_ATP_end_reperfusion": n_reperfusion,
        "neuronal_ATP_recovery_percent": n_recovery,

        "astrocytic_ATP_baseline": a_base,
        "astrocytic_ATP_end_ischemia": a_ischemia,
        "astrocytic_ATP_depletion_percent": a_depletion,
        "astrocytic_ATP_end_reperfusion": a_reperfusion,
        "astrocytic_ATP_recovery_percent": a_recovery,

        "glucose_baseline": baseline["extracellular_glucose"],
        "glucose_end_ischemia": end_ischemia["extracellular_glucose"],
        "glucose_end_reperfusion": end_reperfusion["extracellular_glucose"],

        "lactate_baseline": baseline["extracellular_lactate"],
        "lactate_end_ischemia": end_ischemia["extracellular_lactate"],
        "lactate_end_reperfusion": end_reperfusion["extracellular_lactate"],

        "oxygen_baseline": baseline["capillary_oxygen"],
        "oxygen_end_ischemia": end_ischemia["capillary_oxygen"],
        "oxygen_end_reperfusion": end_reperfusion["capillary_oxygen"],
    })

summary = pd.DataFrame(results)

summary.to_csv("ischemia_duration_summary.csv", index=False)

# Display the main ATP findings.
print("\nISCHEMIA DURATION SUMMARY")
print("=" * 65)

print(
    summary[
        [
            "ischemia_duration_s",
            "neuronal_ATP_depletion_percent",
            "astrocytic_ATP_depletion_percent",
            "neuronal_ATP_recovery_percent",
            "astrocytic_ATP_recovery_percent",
        ]
    ].round(2).to_string(index=False)
)

print("\nFull summary saved as ischemia_duration_summary.csv")
