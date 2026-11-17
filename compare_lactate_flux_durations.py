
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# --------------------------------------------------
# 1. LOAD THE FIVE SIMULATION DATASETS
# --------------------------------------------------

durations = [300, 600, 1000, 1500, 2000]

datasets = {}

for duration in durations:

    filename = Path(
        f"lactate_flux_timecourse_{duration}s.csv"
    )

    if not filename.exists():
        raise FileNotFoundError(filename)

    datasets[duration] = pd.read_csv(filename)

print("All five lactate flux datasets loaded successfully.")


# --------------------------------------------------
# 2. DEFINE THE VARIABLES TO COMPARE
# --------------------------------------------------

variables = {
    "extracellular_lactate": (
        "Extracellular lactate concentration",
        "Concentration (model units)",
    ),
    "vLDHn": (
        "Neuronal LDH reaction rate",
        "Rate (model units)",
    ),
    "vLDHg": (
        "Astrocytic LDH reaction rate",
        "Rate (model units)",
    ),
    "vLACne": (
        "Neuron–extracellular lactate transport",
        "Rate (model units)",
    ),
    "vLACge": (
        "Astrocyte–extracellular lactate transport",
        "Rate (model units)",
    ),
   "vLACc": (
    "Capillary–arterial lactate exchange",
    "Rate (model units)",
),
    "net_extracellular_lactate_rate": (
        "Net extracellular lactate balance",
        "Rate (model units)",
    ),
}


# --------------------------------------------------
# 3. PLOT EACH VARIABLE ACROSS ALL FIVE DURATIONS
# --------------------------------------------------

for variable, labels in variables.items():

    title, ylabel = labels

    fig, ax = plt.subplots(figsize=(11, 6))

    for duration, df in datasets.items():

        # Align time to ischemia onset
        relative_time = df["time_s"] - 100.0

        ax.plot(
            relative_time,
            df[variable],
            label=f"{duration} s ischemia",
        )

    ax.axvline(
        0,
        linestyle="--",
        label="Ischemia onset",
    )

    ax.set_xlabel("Time after ischemia onset (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title + " Across Ischemia Durations")

    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()

    output = f"comparison_{variable}.png"

    plt.savefig(output, dpi=300)
    plt.show()

    print("Saved:", output)


# --------------------------------------------------
# 4. CREATE END-OF-ISCHEMIA SUMMARY
# --------------------------------------------------

summary_records = []

for duration, df in datasets.items():

    end_time = 100.0 + duration

    row = df.iloc[
        (df["time_s"] - end_time).abs().argmin()
    ]

    summary_records.append({
        "ischemia_duration_s": duration,
        "extracellular_lactate": row["extracellular_lactate"],
        "neuronal_LDH": row["vLDHn"],
        "astrocytic_LDH": row["vLDHg"],
        "neuron_extracellular_transport": row["vLACne"],
        "astrocyte_extracellular_transport": row["vLACge"],
        "extracellular_capillary_transport": row["vLACec"],
        "capillary_arterial_exchange": row["vLACc"],
        "net_extracellular_lactate_rate":
            row["net_extracellular_lactate_rate"],
    })

summary = pd.DataFrame(summary_records)

summary.to_csv(
    "cross_duration_lactate_summary.csv",
    index=False,
)

print("\n--- CROSS-DURATION END-OF-ISCHEMIA SUMMARY ---")
print(summary.to_string(index=False))

print("\nCross-duration analysis completed.")
