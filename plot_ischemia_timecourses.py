
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------
# Load all five simulation datasets
# -----------------------------------------

durations = [300, 600, 1000, 1500, 2000]

datasets = {}

for duration in durations:
    filename = f"ischemia_reperfusion_{duration}s.csv"

    data = pd.read_csv(filename)

    # Set time zero to ischemia onset (100 seconds).
    data["elapsed_time"] = data["time_s"] - 100.0

    datasets[duration] = data

# -----------------------------------------
# Define variables for plotting
# -----------------------------------------

variables = [
    ("neuronal_ATP", "Neuronal ATP", "ATP concentration (model units)"),
    ("astrocytic_ATP", "Astrocytic ATP", "ATP concentration (model units)"),
    ("extracellular_glucose", "Extracellular Glucose", "Glucose concentration (model units)"),
    ("extracellular_lactate", "Extracellular Lactate", "Lactate concentration (model units)"),
    ("capillary_oxygen", "Capillary Oxygen", "Oxygen concentration (model units)")
]

# -----------------------------------------
# Create five time-course panels
# -----------------------------------------

fig, axes = plt.subplots(3, 2, figsize=(15, 13))

axes = axes.flatten()

for ax, (column, title, ylabel) in zip(axes, variables):

    for duration in durations:

        data = datasets[duration]

        line, = ax.plot(
            data["elapsed_time"],
            data[column],
            linewidth=1.5,
            label=f"{duration} s"
        )

        # Mark the beginning of reperfusion for this simulation.
        ax.axvline(
            x=duration,
            color=line.get_color(),
            linestyle="--",
            alpha=0.35,
            linewidth=1
        )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Time after ischemia onset (seconds)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(title="Ischemia duration")

# Remove the unused sixth panel.
axes[5].axis("off")

fig.suptitle(
    "Dynamic Metabolic Responses During Ischemia-Reperfusion",
    fontsize=16
)

plt.tight_layout()

plt.savefig(
    "figure5_metabolic_timecourses.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 5 saved as figure5_metabolic_timecourses.png")
