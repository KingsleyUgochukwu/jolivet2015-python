
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the extended 1000-second simulation dataset
df = pd.read_csv("ischemia_reperfusion_1000s_metabolic.csv")

# Extract time and oxygen concentrations
time = df["time_s"].to_numpy()

O2n = df["neuronal_oxygen"].to_numpy()
O2g = df["astrocytic_oxygen"].to_numpy()

# Extract mitochondrial oxygen consumption rates
Mito_n = df["neuronal_mitochondrial_oxygen_consumption"].to_numpy()
Mito_g = df["astrocytic_mitochondrial_oxygen_consumption"].to_numpy()

# Numerically estimate the rate of change of cellular oxygen
dO2n_dt = np.gradient(O2n, time)
dO2g_dt = np.gradient(O2g, time)

# Estimate oxygen transfer into each cellular compartment
transfer_n = dO2n_dt + 0.6 * Mito_n
transfer_g = dO2g_dt + 0.6 * Mito_g

# Calculate net oxygen accumulation rates
net_n = dO2n_dt
net_g = dO2g_dt

# Add calculated values to the dataset
df["neuronal_dO2_dt"] = dO2n_dt
df["astrocytic_dO2_dt"] = dO2g_dt
df["neuronal_oxygen_transfer"] = transfer_n
df["astrocytic_oxygen_transfer"] = transfer_g
df["neuronal_net_oxygen_change"] = net_n
df["astrocytic_net_oxygen_change"] = net_g

# Save the extended analysis
df.to_csv("oxygen_transfer_analysis_1000s.csv", index=False)

# Print selected checkpoints
targets = [100, 1100, 1105, 1110, 1200, 1400, 1500, 1600, 1800, 2100]

print("\nOXYGEN TRANSFER ANALYSIS")
print("Estimated transfer rates are in model flux units.\n")

for target in targets:
    idx = np.abs(time - target).argmin()
    print(
        f"Time: {time[idx]:8.3f} s | "
        f"Neuronal transfer: {transfer_n[idx]:.6f} | "
        f"Neuronal consumption: {Mito_n[idx]:.6f} | "
        f"Neuronal dO2/dt: {dO2n_dt[idx]:.6f} | "
        f"Astrocytic transfer: {transfer_g[idx]:.6f} | "
        f"Astrocytic consumption: {Mito_g[idx]:.6f} | "
        f"Astrocytic dO2/dt: {dO2g_dt[idx]:.6f}"
    )

# Plot oxygen transfer rates
plt.figure(figsize=(12, 6))

plt.plot(time, transfer_n, label="Neuronal oxygen transfer")
plt.plot(time, transfer_g, label="Astrocytic oxygen transfer")

plt.axvspan(100, 1100, alpha=0.12, color="grey", label="Ischaemia")
plt.axvspan(1100, 2100, alpha=0.10, color="green", label="Reperfusion")
plt.axvline(1100, linestyle="--", color="black", linewidth=1)

plt.xlabel("Time (s)")
plt.ylabel("Estimated oxygen transfer rate")
plt.title("Estimated Cellular Oxygen Transfer During Ischaemia–Reperfusion")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

plt.savefig("oxygen_transfer_dynamics_1000s.png", dpi=300)
plt.show()

print("\nSaved:")
print("1. oxygen_transfer_analysis_1000s.csv")
print("2. oxygen_transfer_dynamics_1000s.png")