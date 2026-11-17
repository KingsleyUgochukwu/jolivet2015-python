
import pandas as pd
import matplotlib.pyplot as plt

# Load the combined simulation results
data = pd.read_csv("ischemia_duration_summary.csv")

duration = data["ischemia_duration_s"]

# -----------------------------------------
# Figure 1: ATP depletion
# -----------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    duration,
    data["neuronal_ATP_depletion_percent"],
    marker="o",
    linewidth=2,
    label="Neuronal ATP"
)

plt.plot(
    duration,
    data["astrocytic_ATP_depletion_percent"],
    marker="s",
    linewidth=2,
    label="Astrocytic ATP"
)

plt.xlabel("Ischemia duration (seconds)", fontsize=12)
plt.ylabel("ATP depletion (%)", fontsize=12)

plt.title(
    "Effect of Ischemia Duration on Cellular ATP Depletion",
    fontsize=14
)

plt.ylim(0, 105)
plt.xticks(duration)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    "figure1_ATP_depletion.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 1 saved as figure1_ATP_depletion.png")

# -----------------------------------------
# Figure 2: ATP recovery after reperfusion
# -----------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    duration,
    data["neuronal_ATP_recovery_percent"],
    marker="o",
    linewidth=2,
    label="Neuronal ATP"
)

plt.plot(
    duration,
    data["astrocytic_ATP_recovery_percent"],
    marker="s",
    linewidth=2,
    label="Astrocytic ATP"
)

plt.axhline(
    y=100,
    linestyle="--",
    linewidth=1,
    label="100% recovery reference"
)

plt.xlabel("Ischemia duration (seconds)", fontsize=12)
plt.ylabel("Initial ATP depletion recovered (%)", fontsize=12)

plt.title(
    "ATP Recovery Following Reperfusion",
    fontsize=14
)

plt.xticks(duration)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    "figure2_ATP_recovery.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 2 saved as figure2_ATP_recovery.png")

# -----------------------------------------
# Figure 3: Glucose and lactate dynamics
# -----------------------------------------

stages = ["Baseline", "End of ischemia", "End of reperfusion"]

# Extract glucose measurements
glucose_baseline = data["glucose_baseline"]
glucose_ischemia = data["glucose_end_ischemia"]
glucose_reperfusion = data["glucose_end_reperfusion"]

# Extract lactate measurements
lactate_baseline = data["lactate_baseline"]
lactate_ischemia = data["lactate_end_ischemia"]
lactate_reperfusion = data["lactate_end_reperfusion"]

# Create two panels
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Glucose panel
axes[0].plot(
    duration, glucose_baseline,
    marker="o", linewidth=2, label="Baseline"
)

axes[0].plot(
    duration, glucose_ischemia,
    marker="s", linewidth=2, label="End of ischemia"
)

axes[0].plot(
    duration, glucose_reperfusion,
    marker="^", linewidth=2, label="End of reperfusion"
)

axes[0].set_title("Extracellular Glucose")
axes[0].set_xlabel("Ischemia duration (seconds)")
axes[0].set_ylabel("Glucose concentration (model units)")
axes[0].set_xticks(duration)
axes[0].grid(True, alpha=0.3)
axes[0].legend()

# Lactate panel
axes[1].plot(
    duration, lactate_baseline,
    marker="o", linewidth=2, label="Baseline"
)

axes[1].plot(
    duration, lactate_ischemia,
    marker="s", linewidth=2, label="End of ischemia"
)

axes[1].plot(
    duration, lactate_reperfusion,
    marker="^", linewidth=2, label="End of reperfusion"
)

axes[1].set_title("Extracellular Lactate")
axes[1].set_xlabel("Ischemia duration (seconds)")
axes[1].set_ylabel("Lactate concentration (model units)")
axes[1].set_xticks(duration)
axes[1].grid(True, alpha=0.3)
axes[1].legend()

fig.suptitle(
    "Extracellular Glucose and Lactate During Ischemia-Reperfusion",
    fontsize=14
)

plt.tight_layout()

plt.savefig(
    "figure3_glucose_lactate_dynamics.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 3 saved as figure3_glucose_lactate_dynamics.png")

# -----------------------------------------
# Figure 4: Capillary oxygen dynamics
# -----------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    duration,
    data["oxygen_baseline"],
    marker="o",
    linewidth=2,
    label="Baseline"
)

plt.plot(
    duration,
    data["oxygen_end_ischemia"],
    marker="s",
    linewidth=2,
    label="End of ischemia"
)

plt.plot(
    duration,
    data["oxygen_end_reperfusion"],
    marker="^",
    linewidth=2,
    label="End of reperfusion"
)

plt.xlabel("Ischemia duration (seconds)", fontsize=12)

plt.ylabel(
    "Capillary oxygen concentration (model units)",
    fontsize=12
)

plt.title(
    "Capillary Oxygen Dynamics During Ischemia-Reperfusion",
    fontsize=14
)

plt.xticks(duration)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    "figure4_capillary_oxygen.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 4 saved as figure4_capillary_oxygen.png")


