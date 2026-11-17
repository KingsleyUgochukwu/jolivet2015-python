
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the diagnostic dataset
df = pd.read_csv("lactate_flux_timecourse_1000s.csv")

time = df["time_s"].to_numpy()

ischemia_onset = 100
reperfusion_onset = 1100

# --------------------------------------------------
# FIGURE 1: Cellular LDH reaction rates
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(time, df["vLDHn"], label="Neuronal LDH")
ax.plot(time, df["vLDHg"], label="Astrocytic LDH")

ax.axvline(ischemia_onset, linestyle="--", label="Ischemia onset")
ax.axvline(reperfusion_onset, linestyle=":", label="Reperfusion onset")

ax.set_xlabel("Time (s)")
ax.set_ylabel("LDH reaction rate (model units)")
ax.set_title("Neuronal and Astrocytic LDH Reaction Rates")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figure_7_cellular_LDH_fluxes.png", dpi=300)
plt.show()


# --------------------------------------------------
# FIGURE 2: Lactate transport rates
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 5))

for column, label in [
    ("vLACne", "Neuron–extracellular"),
    ("vLACge", "Astrocyte–extracellular"),
    ("vLACgc", "Astrocyte–capillary"),
    ("vLACec", "Extracellular–capillary"),
]:

    ax.plot(time, df[column], label=label)

ax.axvline(ischemia_onset, linestyle="--", label="Ischemia onset")
ax.axvline(reperfusion_onset, linestyle=":", label="Reperfusion onset")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Transport rate (model units)")
ax.set_title("Lactate Transport Between Compartments")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figure_8_lactate_transport_fluxes.png", dpi=300)
plt.show()


# --------------------------------------------------
# FIGURE 3: Extracellular lactate balance
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(
    time,
    df["net_extracellular_lactate_rate"],
    label="Flux-derived net rate",
)

ax.axhline(0, linestyle="--")
ax.axvline(ischemia_onset, linestyle="--", label="Ischemia onset")
ax.axvline(reperfusion_onset, linestyle=":", label="Reperfusion onset")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Net extracellular lactate rate")
ax.set_title("Net Extracellular Lactate Balance")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figure_9_extracellular_lactate_balance.png", dpi=300)
plt.show()


# --------------------------------------------------
# FIGURE 4: Numerical consistency check
# --------------------------------------------------

numerical_derivative = np.gradient(
    df["extracellular_lactate"].to_numpy(),
    time,
)

flux_derivative = df[
    "net_extracellular_lactate_rate"
].to_numpy()

difference = numerical_derivative - flux_derivative

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(
    time,
    numerical_derivative,
    label="Numerical derivative",
)

ax.plot(
    time,
    flux_derivative,
    label="Flux-derived derivative",
    alpha=0.8,
)

ax.axvline(ischemia_onset, linestyle="--", label="Ischemia onset")
ax.axvline(reperfusion_onset, linestyle=":", label="Reperfusion onset")

ax.set_xlabel("Time (s)")
ax.set_ylabel("d[Extracellular lactate]/dt")
ax.set_title("Consistency of Extracellular Lactate Balance")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figure_10_lactate_balance_validation.png", dpi=300)
plt.show()


# --------------------------------------------------
# PRINT NUMERICAL CHECK
# --------------------------------------------------

print("\n--- EXTRACELLULAR LACTATE BALANCE CHECK ---")

print(
    "Mean absolute difference:",
    np.mean(np.abs(difference)),
)

print(
    "Maximum absolute difference:",
    np.max(np.abs(difference)),
)

print(
    "Correlation:",
    np.corrcoef(numerical_derivative, flux_derivative)[0, 1],
)

print("\nFour figures saved successfully.")
