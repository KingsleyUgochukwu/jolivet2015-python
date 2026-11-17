
import pandas as pd
import numpy as np

durations = [300, 600, 1000, 1500, 2000]

results = []

for duration in durations:

    filename = f"lactate_flux_timecourse_{duration}s.csv"
    df = pd.read_csv(filename)

    reperfusion_onset = 100 + duration

    # Exclude the exact boundary point, which belongs
    # to the final ischemia state in our diagnostic.
    reperfusion = df[
        df["time_s"] > reperfusion_onset
    ].copy()

    if reperfusion.empty:
        raise ValueError(
            f"No reperfusion data found for {duration}s"
        )

    # Maximum astrocytic LDH rate
    astro_index = reperfusion["vLDHg"].idxmax()
    astro_peak = df.loc[astro_index]

    # Maximum neuronal LDH rate
    neuron_index = reperfusion["vLDHn"].idxmax()
    neuron_peak = df.loc[neuron_index]

    # Minimum neuronal LDH rate
    neuron_min_index = reperfusion["vLDHn"].idxmin()
    neuron_min = df.loc[neuron_min_index]

    # Maximum extracellular lactate concentration
    lactate_index = reperfusion[
        "extracellular_lactate"
    ].idxmax()
    lactate_peak = df.loc[lactate_index]

    # Minimum extracellular lactate concentration
    lactate_min_index = reperfusion[
        "extracellular_lactate"
    ].idxmin()
    lactate_min = df.loc[lactate_min_index]

    results.append({
        "ischemia_duration_s": duration,

        "astrocytic_LDH_peak": astro_peak["vLDHg"],
        "astrocytic_LDH_peak_time_after_reperfusion":
            astro_peak["time_s"] - reperfusion_onset,

        "neuronal_LDH_max": neuron_peak["vLDHn"],
        "neuronal_LDH_max_time_after_reperfusion":
            neuron_peak["time_s"] - reperfusion_onset,

        "neuronal_LDH_min": neuron_min["vLDHn"],
        "neuronal_LDH_min_time_after_reperfusion":
            neuron_min["time_s"] - reperfusion_onset,

        "extracellular_lactate_max":
            lactate_peak["extracellular_lactate"],
        "extracellular_lactate_min":
            lactate_min["extracellular_lactate"],

    })

summary = pd.DataFrame(results)

summary.to_csv(
    "reperfusion_transient_summary.csv",
    index=False,
)

print("\n--- REPERFUSION TRANSIENT SUMMARY ---")
print(summary.to_string(index=False))

print("\nResults saved to reperfusion_transient_summary.csv")
