
import numpy as np
import matplotlib.pyplot as plt

stimulations = [120, 180, 240, 300, 360]

plt.figure(figsize=(10, 6))

for Ne in stimulations:
    data = np.load(f"recovery_Ne_{Ne}_600s.npz")

    time = data["time"]
    pcr = data["states"][:, 16]

    # Express PCr as a percentage of its initial value
    pcr_percent = 100 * pcr / pcr[0]

    plt.plot(time, pcr_percent, label=f"Ne = {Ne}")

plt.axvline(
    x=20,
    color="black",
    linestyle="--",
    label="End of stimulation"
)

plt.xlabel("Time (s)")
plt.ylabel("Neuronal PCr (% of baseline)")
plt.title("Effect of Stimulation Intensity on Neuronal PCr Recovery")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig("PCr_recovery_comparison.png", dpi=300)
plt.show()
