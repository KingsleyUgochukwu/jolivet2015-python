"""Create a matched-time extracellular lactate comparison during reperfusion.

Inputs:
  lactate_flux_timecourse_{duration}s.csv for durations 300, 600, 1000, 1500, 2000.

The script interpolates model outputs at 0, 60, 300, 600 and 900 s after
reperfusion onset and writes lactate_fixed_reperfusion_comparison.csv.
"""
from pathlib import Path
import numpy as np
import pandas as pd

DURATIONS = [300, 600, 1000, 1500, 2000]
REPERFUSION_TIMES = [0, 60, 300, 600, 900]
ISCHEMIA_ONSET = 100.0
VE, VN, VG = 0.2, 0.45, 0.25
REN = VE / VN
REG = VE / VG

REQUIRED = [
    "time_s", "extracellular_lactate", "vLACne",
    "vLACge", "vLACec", "net_extracellular_lactate_rate",
]

results = []
for duration in DURATIONS:
    filename = Path(f"lactate_flux_timecourse_{duration}s.csv")
    if not filename.exists():
        raise FileNotFoundError(f"Missing input dataset: {filename}")
    data = pd.read_csv(filename)
    missing = sorted(set(REQUIRED) - set(data.columns))
    if missing:
        raise ValueError(f"{filename} is missing columns: {missing}")

    reperfusion_onset = ISCHEMIA_ONSET + duration
    post = data.loc[data["time_s"] >= reperfusion_onset].copy()
    if post.empty:
        raise ValueError(f"No reperfusion data found in {filename}")
    relative_time = post["time_s"].to_numpy() - reperfusion_onset
    if relative_time[-1] < max(REPERFUSION_TIMES):
        raise ValueError(f"{filename} does not cover 900 s of reperfusion")

    # Adaptive ODE output times do not necessarily coincide with target times.
    for t in REPERFUSION_TIMES:
        def at(column):
            return float(np.interp(t, relative_time, post[column].to_numpy()))

        neuronal = at("vLACne") / REN
        astrocytic = at("vLACge") / REG
        vascular = -at("vLACec")
        calculated = neuronal + astrocytic + vascular
        model = at("net_extracellular_lactate_rate")

        results.append({
            "ischemia_duration_s": duration,
            "time_after_reperfusion_s": t,
            "extracellular_lactate": at("extracellular_lactate"),
            "neuronal_contribution": neuronal,
            "astrocytic_contribution": astrocytic,
            "vascular_removal": vascular,
            "calculated_net_balance": calculated,
            "model_net_balance": model,
            "balance_difference": calculated - model,
        })

summary = pd.DataFrame(results)
summary.to_csv("lactate_fixed_reperfusion_comparison.csv", index=False)
print(summary.round(8).to_string(index=False))
print("\\nSaved lactate_fixed_reperfusion_comparison.csv")
