
import numpy as np
import matplotlib.pyplot as plt

stimulations = [120, 180, 240, 300, 360]

plt.figure(figsize=(10, 6))

for Ne in stimulations:

    data = np.load(f"recovery_Ne_{Ne}_600s.npz")

    time = data["time"]
    pcr = data["states"][:, 16]

    # Initial and minimum phosphocreatine concentrations
    pcr_initial = pcr[0]
    min_index = np.argmin(pcr)
    pcr_min = pcr[min_index]

    # Calculate percentage of depleted PCr replenished
    recovery = (
        (pcr[min_index:] - pcr_min)
        / (pcr_initial - pcr_min)
    ) * 100

    # Time elapsed since each PCr minimum
    recovery_time = time[min_index:] - time[min_index]

    plt.plot(
        recovery_time,
        recovery,
        label=f"Ne = {Ne}"
    )

plt.xlabel("Time Since PCr Minimum (s)")
plt.ylabel("Depleted PCr Replenished (%)")
plt.title("Normalized Neuronal Phosphocreatine Recovery")
plt.legend()
plt.grid(alpha=0.3)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.tight_layout()

plt.savefig("Normalized_PCr_recovery.png", dpi=300)
print("\nTIME TO REPLENISH 40% OF DEPLETED PCr")
print("-" * 45)

for Ne in stimulations:
    data = np.load(f"recovery_Ne_{Ne}_600s.npz")

    time = data["time"]
    pcr = data["states"][:, 16]

    idx = np.argmin(pcr)

    recovery = 100 * (
        pcr[idx:] - pcr[idx]
    ) / (pcr[0] - pcr[idx])

    recovery_time = time[idx:] - time[idx]

    crossing = np.where(recovery >= 40)[0]

    if len(crossing):
        print(
            f"Ne={Ne}: "
            f"40% replenishment at "
            f"{recovery_time[crossing[0]]:.2f} s "
            f"after PCr minimum"
        )
    else:
        print(f"Ne={Ne}: 40% replenishment not reached")
plt.show()
