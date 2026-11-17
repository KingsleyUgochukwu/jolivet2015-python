
import pandas as pd
import numpy as np

durations = [300, 600, 1000, 1500, 2000]

results = []

for duration in durations:

    df = pd.read_csv(
        "cumulative_lactate_contributions.csv"
    )

    row = df[
        df["ischemia_duration_s"] == duration
    ].iloc[0]

    neuronal = row["neuronal_cumulative"]
    astrocytic = row["astrocytic_cumulative"]
    vascular = row["vascular_cumulative"]

    # Total absolute transport magnitude
    total_absolute = (
        abs(neuronal)
        + abs(astrocytic)
        + abs(vascular)
    )

    results.append({
        "ischemia_duration_s": duration,

        "neuronal_signed": neuronal,
        "astrocytic_signed": astrocytic,
        "vascular_signed": vascular,

        "neuronal_absolute_percent":
            abs(neuronal) / total_absolute * 100,

        "astrocytic_absolute_percent":
            abs(astrocytic) / total_absolute * 100,

        "vascular_absolute_percent":
            abs(vascular) / total_absolute * 100,

        "net_balance": (
            neuronal + astrocytic + vascular
        )
    })

summary = pd.DataFrame(results)

summary.to_csv(
    "relative_lactate_contributions.csv",
    index=False
)

print("\n--- RELATIVE CUMULATIVE CONTRIBUTIONS ---")
print(summary.round(4).to_string(index=False))

print("\nAnalysis completed.")
