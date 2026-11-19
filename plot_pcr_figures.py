"""Generate the three core PCr figures from regenerated NPZ datasets.

Run after:
    python generate_pcr_datasets.py

Outputs:
    Phosphocreatine_time_course.png
    PCr_recovery_comparison.png
    Normalized_PCr_recovery.png

The first figure shows the 0-20 s stimulation segment. The latter two use
the full 620 s simulation, including the post-stimulation PCr minimum.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

STIMULATIONS = [120, 180, 240, 300, 360]

# Figure 1: PCr concentration during the 20-second stimulation interval.
fig, ax = plt.subplots(figsize=(9, 5))
for ne in STIMULATIONS:
    data = np.load(f"stimulation_Ne_{ne}.npz")
    ax.plot(data["time"], data["states"][:, 16], label=f"Ne = {ne}")
ax.set_xlabel("Time (seconds)")
ax.set_ylabel("Neuronal phosphocreatine concentration")
ax.set_title("Phosphocreatine Dynamics During Neuronal Stimulation")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("Phosphocreatine_time_course.png", dpi=300)
plt.close(fig)

# Figure 2: Absolute PCr level expressed relative to the pre-stimulation baseline.
fig, ax = plt.subplots(figsize=(9, 5))
for ne in STIMULATIONS:
    data = np.load(f"recovery_Ne_{ne}_600s.npz")
    time = data["time"]
    pcr = data["states"][:, 16]
    baseline = pcr[0]
    ax.plot(time, 100.0 * pcr / baseline, label=f"Ne = {ne}")
ax.axvline(20.0, color="black", linestyle="--", label="End of stimulation")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Neuronal PCr (% of baseline)")
ax.set_title("Effect of Stimulation Intensity on Neuronal PCr Recovery")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("PCr_recovery_comparison.png", dpi=300)
plt.close(fig)

# Figure 3: Fraction of the depleted PCr pool replenished after its minimum.
fig, ax = plt.subplots(figsize=(9, 5))
for ne in STIMULATIONS:
    data = np.load(f"recovery_Ne_{ne}_600s.npz")
    time = data["time"]
    pcr = data["states"][:, 16]
    min_index = int(np.argmin(pcr))
    initial = pcr[0]
    minimum = pcr[min_index]
    if np.isclose(initial, minimum):
        raise ValueError(f"No measurable PCr depletion for Ne={ne}")
    replenished = 100.0 * (pcr[min_index:] - minimum) / (initial - minimum)
    relative_time = time[min_index:] - time[min_index]
    ax.plot(relative_time, replenished, label=f"Ne = {ne}")
ax.set_xlabel("Time Since PCr Minimum (s)")
ax.set_ylabel("Depleted PCr Replenished (%)")
ax.set_title("Normalized Neuronal Phosphocreatine Recovery")
ax.legend()
ax.grid(alpha=0.3)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig("Normalized_PCr_recovery.png", dpi=300)
plt.close(fig)

print("Generated:")
print("  Phosphocreatine_time_course.png")
print("  PCr_recovery_comparison.png")
print("  Normalized_PCr_recovery.png")
