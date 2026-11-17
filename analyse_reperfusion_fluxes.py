
import pandas as pd

durations = [300, 600, 1000, 1500, 2000]

results = []

for duration in durations:

    filename = f"lactate_flux_timecourse_{duration}s.csv"
    df = pd.read_csv(filename)

    reperfusion_onset = 100 + duration

    # State immediately before reperfusion
    pre = df[
        df["time_s"] <= reperfusion_onset
    ].iloc[-1]

    # Reperfusion period
    reperfusion = df[
        df["time_s"] > reperfusion_onset
    ].copy()

    if reperfusion.empty:
        raise ValueError(
            f"No reperfusion data for {duration}s"
        )

    # Time of maximum extracellular lactate
    peak_index = reperfusion[
        "extracellular_lactate"
    ].idxmax()

    peak = df.loc[peak_index]

    # Sample 500 seconds after reperfusion
    later_time = reperfusion_onset + 500

    later_index = (
        df["time_s"] - later_time
    ).abs().idxmin()

    later = df.loc[later_index]

    sampling_points = {
        "Before reperfusion": pre,
        "Lactate rebound peak": peak,
        "500 s after reperfusion": later,
    }

    for stage, row in sampling_points.items():

        results.append({
            "ischemia_duration_s": duration,
            "stage": stage,
            "time_s": row["time_s"],
            "extracellular_lactate":
                row["extracellular_lactate"],
            "vLDHn": row["vLDHn"],
            "vLDHg": row["vLDHg"],
            "vLACne": row["vLACne"],
            "vLACge": row["vLACge"],
            "vLACgc": row["vLACgc"],
            "vLACec": row["vLACec"],
            "vLACc": row["vLACc"],
            "net_extracellular_lactate_rate":
                row["net_extracellular_lactate_rate"],
        })


# --------------------------------------------------
# SAVE AND DISPLAY RESULTS
# --------------------------------------------------

summary = pd.DataFrame(results)

summary.to_csv(
    "reperfusion_lactate_flux_comparison.csv",
    index=False,
)

print("\n--- REPERFUSION LACTATE FLUX COMPARISON ---")

for duration in durations:

    subset = summary[
        summary["ischemia_duration_s"] == duration
    ]

    print(f"\nISCHEMIA DURATION: {duration} seconds")

    print(
        subset.to_string(
            index=False,
            float_format=lambda x: f"{x:.7f}"
        )
    )

print("\nResults saved successfully.")
