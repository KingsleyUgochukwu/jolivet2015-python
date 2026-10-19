
import numpy as np
import matplotlib.pyplot as plt

# Load the existing 3600-second simulation
data = np.load("recovery_Ne_240_3600s.npz")

time = data["time"]
states = data["states"]

# ATP and PCr concentrations
atp_n = states[:, 14]
atp_g = states[:, 15]
pcr_n = states[:, 16]
pcr_g = states[:, 17]

# Create Figure 3 with three panels
fig = plt.figure(figsize=(11, 7))

gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1, 1],
    hspace=0.38,
    wspace=0.28
)

ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])

# Panel A: ATP and PCr over 1500 seconds
ax_a.plot(time, pcr_n, color="blue", linewidth=1)
ax_a.plot(time, pcr_g, color="red", linewidth=1)
ax_a.plot(time, atp_n, color="blue", linewidth=1)
ax_a.plot(time, atp_g, color="red", linewidth=1)

ax_a.set_xlim(-100, 1500)
ax_a.set_ylim(2.0, 5.5)
ax_a.set_ylabel("ATP / PCr (mM)")
ax_a.set_xlabel("Time (s)")
ax_a.set_title("A. ATP and PCr dynamics")

# Panel B: PCr during stimulation and early recovery
ax_b.plot(time, pcr_n, color="blue", label="Neuron")
ax_b.plot(time, pcr_g, color="red", label="Astrocyte")

ax_b.set_xlim(-10, 60)
ax_b.set_ylim(4.6, 5.0)
ax_b.set_xlabel("Time (s)")
ax_b.set_ylabel("PCr (mM)")
ax_b.set_title("B. Phosphocreatine dynamics")
ax_b.legend()

# Panel C: ATP during stimulation and early recovery
ax_c.plot(time, atp_n, color="blue", label="Neuron")
ax_c.plot(time, atp_g, color="red", label="Astrocyte")

ax_c.set_xlim(-10, 60)
ax_c.set_ylim(2.0, 2.3)
ax_c.set_xlabel("Time (s)")
ax_c.set_ylabel("ATP (mM)")
ax_c.set_title("C. ATP dynamics")
ax_c.legend()

# Highlight the 20-second stimulation
for ax in [ax_a, ax_b, ax_c]:
    ax.axvspan(0, 20, color="lightgray", alpha=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.savefig(
    "Figure_03_Python_initial.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
