
import pandas as pd
import numpy as np

durations = [300, 600, 1000, 1500, 2000]

REN = 0.4444444444444445
REG = 0.8

results = []

for duration in durations:

    df = pd.read_csv(
        f"lactate_flux_timecourse_{duration}s.csv"
    )

    reperfusion_onset = 100 + duration

    # Select the complete reperfusion period
    df = df[
        df["time_s"] >= reperfusion_onset
    ].copy()

    time = (
        df["time_s"].to_numpy()
        - reperfusion_onset
    )

    # Calculate individual contributions
    neuronal = (
        df["vLACne"].to_numpy() / REN
    )

    astrocytic = (
        df["vLACge"].to_numpy() / REG
    )

    vascular = (
        -df["vLACec"].to_numpy()
    )

    net = (
        df["net_extracellular_lactate_rate"].to_numpy()
    )

    # Integrate contributions over reperfusion
    neuronal_total = np.trapezoid(neuronal, time)
    astrocytic_total = np.trapezoid(astrocytic, time)
    vascular_total = np.trapezoid(vascular, time)
    net_total = np.trapezoid(net, time)

    # Actual extracellular lactate concentration change
    lactate_change = (
        df["extracellular_lactate"].iloc[-1]
        - df["extracellular_lactate"].iloc[0]
    )

    results.append({
        "ischemia_duration_s": duration,
        "neuronal_cumulative": neuronal_total,
        "astrocytic_cumulative": astrocytic_total,
        "vascular_cumulative": vascular_total,
        "integrated_net_balance": net_total,
        "actual_lactate_change": lactate_change,
        "integration_difference": net_total - lactate_change
    })

summary = pd.DataFrame(results)

summary.to_csv(
    "cumulative_lactate_contributions.csv",
    index=False
)

print("\n--- CUMULATIVE LACTATE CONTRIBUTIONS ---")
print(summary.round(6).to_string(index=False))

print("\nAnalysis completed.")
