
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# --------------------------------------------------
# File paths
# --------------------------------------------------

MODEL_FILE = Path("ischemia_reperfusion_1000s.csv")
EXPERIMENT_FILE = Path(
    "data/experimental/liu_2012_figure3_tissue_pO2_digitised.csv"
)

# --------------------------------------------------
# Load datasets
# --------------------------------------------------

model = pd.read_csv(MODEL_FILE)
experiment = pd.read_csv(EXPERIMENT_FILE)

# Model protocol:
# Ischemia: 100 to 1100 seconds
# Reperfusion: 1100 to 2100 seconds

BASELINE_TIME = 100
REPERFUSION_TIME = 1100
END_TIME = 2100

baseline = model.loc[
    model["time_s"].sub(BASELINE_TIME).abs().idxmin()
]

# Percentage change relative to pre-ischemia baseline
for variable in ["neuronal_oxygen", "astrocytic_oxygen"]:
    model[variable + "_change_percent"] = (
        (model[variable] - baseline[variable])
        / baseline[variable]
    ) * 100

# --------------------------------------------------
# Normalise time within each experimental phase
# --------------------------------------------------

model["phase"] = np.where(
    model["time_s"] < REPERFUSION_TIME,
    "MCAO",
    "Reperfusion"
)

model["phase_progress_percent"] = np.where(
    model["phase"] == "MCAO",
    (model["time_s"] - BASELINE_TIME)
    / (REPERFUSION_TIME - BASELINE_TIME) * 100,
    (model["time_s"] - REPERFUSION_TIME)
    / (END_TIME - REPERFUSION_TIME) * 100
)

experiment["phase_progress_percent"] = np.where(
    experiment["phase"] == "MCAO",
    experiment["phase_time_min"] / 90 * 100,
    experiment["phase_time_min"] / 30 * 100
)

# --------------------------------------------------
# Save combined model percentage-change data
# --------------------------------------------------

model.to_csv(
    "model_oxygen_percentage_changes_1000s.csv",
    index=False
)

# --------------------------------------------------
# Plot comparison
# --------------------------------------------------

fig, axes = plt.subplots(
    1, 2,
    figsize=(15, 6),
    sharey=True
)

groups = [
    ("MCAO", "Ischaemia"),
    ("Reperfusion", "Reperfusion")
]

for ax, (phase, title) in zip(axes, groups):

    model_phase = model[model["phase"] == phase]
    exp_phase = experiment[experiment["phase"] == phase]

    # Model predictions
    ax.plot(
        model_phase["phase_progress_percent"],
        model_phase["neuronal_oxygen_change_percent"],
        label="Model: neuronal oxygen",
        linewidth=2
    )

    ax.plot(
        model_phase["phase_progress_percent"],
        model_phase["astrocytic_oxygen_change_percent"],
        label="Model: astrocytic oxygen",
        linewidth=2
    )

    # Experimental observations
    for group in ["UNICAR", "BICAR"]:
        group_data = exp_phase[exp_phase["group"] == group]

        ax.plot(
            group_data["phase_progress_percent"],
            group_data["pO2_change_percent_from_baseline"],
            marker="o",
            linestyle="--",
            label=f"Liu et al. 2012: {group}"
        )

    ax.set_title(title)
    ax.set_xlabel("Progress through phase (%)")
    ax.set_xlim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

axes[0].set_ylabel("Change from pre-ischaemia baseline (%)")

fig.suptitle(
    "Modelled Cellular Oxygen and Experimental Tissue pO2",
    fontsize=13
)

fig.text(
    0.5, 0.01,
    "Model: 1,000 s ischaemia + 1,000 s reperfusion. "
    "Liu et al. (2012): 90 min MCAO + 30 min reperfusion. "
    "Experimental values are approximate digitised estimates.",
    ha="center",
    fontsize=8
)

plt.tight_layout(rect=[0, 0.06, 1, 0.93])

plt.savefig(
    "model_vs_liu2012_oxygen_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Comparison completed successfully.")
print("Generated:")
print("- model_oxygen_percentage_changes_1000s.csv")
print("- model_vs_liu2012_oxygen_comparison.png")
