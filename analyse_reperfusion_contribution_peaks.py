
import pandas as pd

durations = [300, 600, 1000, 1500, 2000]

REN = 0.4444444444444445
REG = 0.8

results = []

for duration in durations:

    df = pd.read_csv(
        f"lactate_flux_timecourse_{duration}s.csv"
    )

    reperfusion_onset = 100 + duration

    df = df[df["time_s"] >= reperfusion_onset].copy()

    df["time_after_reperfusion"] = (
        df["time_s"] - reperfusion_onset
    )

    df["neuronal"] = df["vLACne"] / REN
    df["astrocytic"] = df["vLACge"] / REG
    df["vascular"] = -df["vLACec"]

    variables = [
        "neuronal",
        "astrocytic",
        "vascular",
        "net_extracellular_lactate_rate"
    ]

    for variable in variables:

        maximum = df[variable].idxmax()
        minimum = df[variable].idxmin()

        results.append({
            "ischemia_duration_s": duration,
            "contribution": variable,
            "maximum_value": df.loc[maximum, variable],
            "time_of_maximum_s": df.loc[
                maximum, "time_after_reperfusion"
            ],
            "minimum_value": df.loc[minimum, variable],
            "time_of_minimum_s": df.loc[
                minimum, "time_after_reperfusion"
            ]
        })

summary = pd.DataFrame(results)

summary.to_csv(
    "reperfusion_contribution_peaks.csv",
    index=False
)

print("\n--- REPERFUSION CONTRIBUTION PEAK ANALYSIS ---")
print(summary.round(6).to_string(index=False))

print("\nAnalysis completed.")
