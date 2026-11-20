import pandas as pd
import matplotlib.pyplot as plt

# Load model-derived metabolite summary
model = pd.read_csv("metabolite_response_summary.csv")

# Convert ischaemia duration to minutes
duration = model["ischemia_duration_s"] / 60

# Create two-panel figure
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# -----------------------------------
# Panel A: Extracellular glucose
# -----------------------------------
axes[0].plot(
    duration,
    model["glucose_end_ischemia_percent_baseline"],
    marker="o",
    linewidth=2,
    label="End of ischaemia"
)

axes[0].plot(
    duration,
    model["glucose_end_reperfusion_percent_baseline"],
    marker="s",
    linewidth=2,
    label="End of reperfusion"
)

axes[0].set_title("Extracellular glucose")
axes[0].set_xlabel("Ischaemia duration (min)")
axes[0].set_ylabel("Glucose (% of baseline)")
axes[0].grid(alpha=0.3)
axes[0].legend()

# -----------------------------------
# Panel B: Extracellular lactate
# -----------------------------------
axes[1].plot(
    duration,
    model["lactate_end_ischemia_percent_baseline"],
    marker="o",
    linewidth=2,
    label="End of ischaemia"
)

axes[1].plot(
    duration,
    model["lactate_end_reperfusion_percent_baseline"],
    marker="s",
    linewidth=2,
    label="End of reperfusion"
)

axes[1].set_title("Extracellular lactate")
axes[1].set_xlabel("Ischaemia duration (min)")
axes[1].set_ylabel("Lactate (% of baseline)")
axes[1].grid(alpha=0.3)
axes[1].legend()

# -----------------------------------
# Figure formatting and saving
# -----------------------------------
fig.suptitle(
    "Model-predicted extracellular glucose and lactate responses",
    fontsize=13
)

plt.tight_layout()

plt.savefig(
    "figure6_model_metabolite_responses.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Model metabolite response figure saved successfully.")