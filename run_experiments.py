
import numpy as np

from simulation_runner import run_simulation


# Stimulation conditions
Ne_values = [120, 180, 240, 300, 360]

# Store results for each condition
results = {}

for Ne in Ne_values:

    print()
    print("=" * 50)
    print(f"RUNNING EXPERIMENT: Ne = {Ne}")
    print("=" * 50)

    time, states = run_simulation(
        Ne=Ne,
        tend=20,
        verbose=False,
    )

    results[Ne] = {
        "time": time,
        "states": states,
    }
    # Save complete simulation trajectory for this condition
    np.savez_compressed(
        f"stimulation_Ne_{Ne}.npz",
        time=time,
        states=states,
        Ne=Ne,
)
    
    print(f"Saved stimulation_Ne_{Ne}.npz")
    print(f"ATP at 20 s: {states[-1, 14]:.6f}")
    print(f"Na+ at 20 s: {states[-1, 0]:.6f}")
    print(f"Minimum ATP: {np.min(states[:, 14]):.6f}")
    print(f"Maximum Na+: {np.max(states[:, 0]):.6f}")

    # Save the experiment summary for later analysis
summary = []

for Ne, result in results.items():
    states = result["states"]

    summary.append([
        Ne,
        states[-1, 14],
        np.min(states[:, 14]),
        states[-1, 0],
        np.max(states[:, 0]),
    ])

np.savetxt(
    "stimulation_summary.csv",
    np.array(summary),
    delimiter=",",
    header="Ne,ATP_20s,ATP_min,Na_20s,Na_max",
    comments="",
    fmt="%.8f",
)

print("Summary saved to stimulation_summary.csv")
print()
print("All five stimulation experiments completed.")
