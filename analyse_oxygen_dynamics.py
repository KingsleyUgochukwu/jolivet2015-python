
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# Load simulation data
# --------------------------------------------------

data = pd.read_csv("ischemia_reperfusion_1000s.csv")

# Time relative to the beginning of ischemia
data["elapsed_time_s"] = data["time_s"] - 100.0

# --------------------------------------------------
# Calculate oxygen-dependent mitochondrial saturation
# --------------------------------------------------

KO2_MITO = 0.001  # mM

data["neuronal_O2_saturation_percent"] = (
    data["neuronal_oxygen"]
    / (data["neuronal_oxygen"] + KO2_MITO)
) * 100

data["astrocytic_O2_saturation_percent"] = (
    data["astrocytic_oxygen"]
    / (data["astrocytic_oxygen"] + KO2_MITO)
) * 100

# --------------------------------------------------
# Save derived dataset
# --------------------------------------------------

data.to_csv(
    "oxygen_dynamics_analysis_1000s.csv",
    index=False
)

# --------------------------------------------------
# Plot 1: Cellular oxygen concentrations
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(11, 6))

ax.plot(
    data["elapsed_time_s"],
    data["neuronal_oxygen"],
    label="Neuronal oxygen"
)

ax.plot(
    data["elapsed_time_s"],
    data["astrocytic_oxygen"],
    label="Astrocytic oxygen"
)
# Shade the experimental phases
ax.axvspan(
    0, 1000,
    color="lightcoral",
    alpha=0.15,
    label="Ischaemia"
)

ax.axvspan(
    1000, 2000,
    color="lightblue",
    alpha=0.15,
    label="Reperfusion"
)

ax.axvline(
    1000,
    color="black",
    linestyle="--",
    linewidth=1.2
)

ax.axvline(
    1000,
    linestyle="--",
    label="Reperfusion begins"
)

ax.set_xlabel("Time from ischaemia onset (s)")
ax.set_ylabel("Oxygen concentration (mM)")
ax.set_title(
    "Neuronal and Astrocytic Oxygen Concentrations "
    "During Ischaemia–Reperfusion"
)
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("cellular_oxygen_dynamics_1000s.png", dpi=300)
plt.close()

# --------------------------------------------------
# Plot 2: Oxygen-dependent mitochondrial saturation
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(11, 6))

ax.plot(
    data["elapsed_time_s"],
    data["neuronal_O2_saturation_percent"],
    label="Neuronal"
)

ax.plot(
    data["elapsed_time_s"],
    data["astrocytic_O2_saturation_percent"],
    label="Astrocytic"
)
# Shade the experimental phases
ax.axvspan(
    0, 1000,
    color="lightcoral",
    alpha=0.15,
    label="Ischaemia"
)

ax.axvspan(
    1000, 2000,
    color="lightblue",
    alpha=0.15,
    label="Reperfusion"
)

ax.axvline(
    1000,
    color="black",
    linestyle="--",
    linewidth=1.2
)

ax.axvline(
    1000,
    linestyle="--",
    label="Reperfusion begins"
)

ax.set_xlabel("Time from ischaemia onset (s)")
ax.set_ylabel("Oxygen-dependent saturation (%)")
ax.set_title(
    "Oxygen-Dependent Mitochondrial Saturation "
    "During Ischaemia–Reperfusion"
)
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("mitochondrial_oxygen_saturation_1000s.png", dpi=300)
plt.close()

# --------------------------------------------------
# Print summary
# --------------------------------------------------

for variable in [
    "neuronal_oxygen",
    "astrocytic_oxygen",
    "capillary_oxygen",
    "neuronal_O2_saturation_percent",
    "astrocytic_O2_saturation_percent",
]:

    print(f"\n{variable}")

    print(
        "Minimum:",
        data[variable].min()
    )

    print(
        "Maximum:",
        data[variable].max()
    )

    print(
        "Initial:",
        data[variable].iloc[0]
    )

    print(
        "Final:",
        data[variable].iloc[-1]
    )

print("\nOxygen analysis completed successfully.")

# --------------------------------------------------
# Identify minimum oxygen points and phase values
# --------------------------------------------------

print("\n--- OXYGEN MINIMUM ANALYSIS ---")

for variable in [
    "neuronal_oxygen",
    "astrocytic_oxygen",
    "capillary_oxygen"
]:
    idx = data[variable].idxmin()
    row = data.loc[idx]

    print(f"\n{variable}")
    print(f"Minimum value: {row[variable]:.8f} mM")
    print(f"Absolute simulation time: {row['time_s']:.2f} s")
    print(f"Time from ischemia onset: {row['elapsed_time_s']:.2f} s")

# --------------------------------------------------
# Oxygen values at important time points
# --------------------------------------------------

print("\n--- OXYGEN AT KEY TIME POINTS ---")

key_times = {
    "Pre-ischemia": 100,
    "End of ischemia / start of reperfusion": 1100,
    "End of reperfusion": 2100
}

for label, target_time in key_times.items():
    idx = (data["time_s"] - target_time).abs().idxmin()
    row = data.loc[idx]

    print(f"\n{label} ({row['time_s']:.2f} s)")
    print(f"Neuronal oxygen: {row['neuronal_oxygen']:.8f} mM")
    print(f"Astrocytic oxygen: {row['astrocytic_oxygen']:.8f} mM")
    print(f"Capillary oxygen: {row['capillary_oxygen']:.8f} mM")
