
import numpy as np
import matplotlib.pyplot as plt

from conductance import excitatory_conductance

# Load the validated full simulation trajectory
data = np.load("recovery_Ne_240_3600s.npz")

time = data["time"]
states = data["states"]

# Extract the variables used in MATLAB Figure 2
voltage = states[:, 27]
sodium_n = states[:, 0]
sodium_g = states[:, 1]

# Reconstruct the excitatory conductance
conductance = np.array([
    excitatory_conductance(
        t=float(t),
        t1=0.0,
        tstim=20.0,
        Ne=240,
        LONG="off"
    )
    for t in time
])

# Plot the three main panels
fig, axes = plt.subplots(
    3, 1, figsize=(8, 8), sharex=True
)

# Panel A: Excitatory conductance
axes[0].plot(time, conductance, color="black", linewidth=1)
axes[0].set_ylabel("Conductance")
axes[0].set_ylim(0, 0.07)
axes[0].set_title("A. Excitatory synaptic conductance")

# Panel B: Neuronal membrane potential
axes[1].plot(time, voltage, color="blue", linewidth=0.7)
axes[1].set_ylabel("Membrane potential (mV)")
axes[1].set_ylim(-100, 100)
axes[1].set_title("B. Neuronal membrane potential")

# Panel C: Intracellular sodium
axes[2].plot(time, sodium_n, color="blue", label="Neuron")
axes[2].plot(time, sodium_g, color="red", label="Glia")
axes[2].set_ylabel("Na+ concentration (mM)")
axes[2].set_ylim(5, 25)
axes[2].set_title("C. Intracellular sodium")
axes[2].legend()

# Match the MATLAB Figure 2 time range
axes[2].set_xlim(-10, 60)
axes[2].set_xlabel("Time (seconds)")

for ax in axes:
    ax.axvspan(0, 20, color="lightgray", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.tight_layout()

# --------------------------------------------------
# Figure 2: Instantaneous firing-frequency inset
# --------------------------------------------------

# Detect upward crossings of the MATLAB threshold
threshold = -20.0

crossings = (
    (time[:-1] < 20.0)
    & (voltage[:-1] < threshold)
    & (voltage[1:] >= threshold)
)

spike_times = time[1:][crossings]

# Calculate interspike intervals
isi = np.diff(spike_times)

# MATLAB uses the earlier spike time for each frequency
frequency_times = spike_times[:-1]
instantaneous_frequency = 1.0 / isi

# Reproduce MATLAB's appended zero-frequency points
frequency_times = np.concatenate([
    frequency_times,
    [frequency_times[-1], 20.0]
])

instantaneous_frequency = np.concatenate([
    instantaneous_frequency,
    [0.0, 0.0]
])

# Add the inset to Panel A
inset = axes[0].inset_axes([0.66, 0.35, 0.30, 0.55])

inset.plot(
    frequency_times,
    instantaneous_frequency,
    color="blue",
    linewidth=1
)

inset.set_xlim(-2, 20)
inset.set_ylim(0, 140)
inset.set_xticks([0, 10, 20])
inset.set_yticks([0, 50, 100, 140])
inset.set_xlabel("Time (s)", fontsize=8)
inset.set_ylabel("Hz", fontsize=8)
inset.tick_params(labelsize=8)
inset.set_title("Firing frequency", fontsize=9)


plt.savefig(
    "Figure_02_Python_initial.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
