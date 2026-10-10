
import numpy as np

stimulations = [120, 180, 240, 300, 360]
recovery_curves = {}

# Common recovery time grid (seconds after PCr minimum)
common_time = np.arange(0, 541, 1)

for Ne in stimulations:
    data = np.load(f"recovery_Ne_{Ne}_600s.npz")

    time = data["time"]
    pcr = data["states"][:, 16]

    # Locate the minimum PCr concentration
    min_index = np.argmin(pcr)
    pcr_initial = pcr[0]
    pcr_min = pcr[min_index]

    # Time elapsed since PCr minimum
    recovery_time = time[min_index:] - time[min_index]

    # Percentage of depleted PCr replenished
    normalized_recovery = (
        100 * (pcr[min_index:] - pcr_min)
        / (pcr_initial - pcr_min)
    )

    # Interpolate onto a common time grid
    recovery_curves[Ne] = np.interp(
        common_time,
        recovery_time,
        normalized_recovery
    )

# Reference recovery curve
reference = recovery_curves[240]

print("\nNORMALIZED PCr RECOVERY RMSD")
print("Reference stimulation: Ne=240")
print("-" * 45)

for Ne in stimulations:
    difference = recovery_curves[Ne] - reference

    rmsd = np.sqrt(np.mean(difference ** 2))
    max_difference = np.max(np.abs(difference))

    print(
        f"Ne={Ne}: "
        f"RMSD={rmsd:.4f} percentage points | "
        f"Maximum difference={max_difference:.4f} percentage points"
    )
