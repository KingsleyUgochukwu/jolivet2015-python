
import pandas as pd

durations = [300, 600, 1000, 1500, 2000]

REN = 0.4444444444444445
REG = 0.8

sampling_times = [0, 30, 60, 120, 300, 500, 1000]

results = []

for duration in durations:

    df = pd.read_csv(
        f"lactate_flux_timecourse_{duration}s.csv"
    )

    reperfusion_onset = 100 + duration

    df["time_after_reperfusion"] = (
        df["time_s"] - reperfusion_onset
    )

    for target_time in sampling_times:

        # Find the closest available simulation time
        idx = (
            df["time_after_reperfusion"] - target_time
        ).abs().idxmin()

        row = df.loc[idx]

        if row["time_after_reperfusion"] < -1e-8:
            continue

        neuronal = row["vLACne"] / REN
        astrocytic = row["vLACge"] / REG
        vascular = -row["vLACec"]

        results.append({
            "ischemia_duration_s": duration,
            "reperfusion_time_s": row["time_after_reperfusion"],
            "neuronal_contribution": neuronal,
            "astrocytic_contribution": astrocytic,
            "vascular_contribution": vascular,
            "net_balance": (
                neuronal + astrocytic + vascular
            )
        })

summary = pd.DataFrame(results)

summary.to_csv(
    "reperfusion_contribution_summary.csv",
    index=False
)

print("\n--- REPERFUSION CONTRIBUTION SUMMARY ---")
print(summary.round(6).to_string(index=False))

print("\nAnalysis completed.")
