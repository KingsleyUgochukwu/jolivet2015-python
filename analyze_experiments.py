
import numpy as np
import matplotlib.pyplot as plt

# Results from the five completed simulations
Ne = np.array([120, 180, 240, 300, 360])

ATP_min = np.array([
    2.174377,
    2.139964,
    2.090334,
    2.029075,
    1.962651,
])

Na_max = np.array([
    9.982607,
    10.996363,
    11.939393,
    12.841264,
    13.713314,
])

# Plot 1: Minimum neuronal ATP
plt.figure(figsize=(8, 5))
plt.plot(Ne, ATP_min, "o-", linewidth=2)

plt.xlabel("Neuronal stimulation parameter (Ne)")
plt.ylabel("Minimum neuronal ATP")
plt.title("Effect of Ne on Neuronal ATP Depletion")
plt.xticks(Ne)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("ATP_vs_stimulation.png", dpi=300)
plt.show()

# Plot 2: Maximum neuronal sodium
plt.figure(figsize=(8, 5))
plt.plot(Ne, Na_max, "o-", linewidth=2)

plt.xlabel("Neuronal stimulation parameter (Ne)")
plt.ylabel("Maximum neuronal Na+")
plt.title("Effect of Ne on Neuronal Sodium Accumulation")
plt.xticks(Ne)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("Na_vs_stimulation.png", dpi=300)
plt.show()

print("Analysis completed successfully.")
print("Saved ATP_vs_stimulation.png")
print("Saved Na_vs_stimulation.png")

from parameters import get_parameters

p = get_parameters()

Ne_values = np.array([120, 180, 240, 300, 360])

ATP_20s = np.array([
    2.190386,
    2.182569,
    2.166067,
    2.123681,
    2.046068,
])

Na_20s = np.array([
    8.829554,
    9.346411,
    9.864685,
    10.519612,
    11.423969,
])

# Same equation used in rates.py
vPumpn_20s = (
    p["SmVn"]
    * p["kPumpn"]
    * ATP_20s
    * Na_20s
    / (1 + ATP_20s / p["KmPump"])
)

print("\nNeuronal Na/K-ATPase activity at 20 seconds")

for Ne, pump in zip(Ne_values, vPumpn_20s):
    print(f"Ne = {Ne}: Pump rate = {pump:.6f}")

plt.figure(figsize=(8, 5))
plt.plot(Ne_values, vPumpn_20s, "o-", linewidth=2)
plt.xlabel("Neuronal stimulation parameter (Ne)")
plt.ylabel("Neuronal Na/K-ATPase rate")
plt.title("Effect of Ne on Sodium Pump Activity at 20 s")
plt.xticks(Ne_values)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("Pump_vs_stimulation.png", dpi=300)
plt.show()

# ============================================================
# TIME-DEPENDENT SODIUM PUMP ANALYSIS
# ============================================================

print("\nTIME-DEPENDENT SODIUM PUMP ANALYSIS")
print("-" * 60)

pump_results = []

plt.figure(figsize=(9, 5))

for Ne_value in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne_value}.npz")

    time = data["time"]
    states = data["states"]

    ATP = states[:, 14]
    sodium = states[:, 0]

    # Neuronal Na+/K+-ATPase equation from rates.py
    pump_rate = (
        p["SmVn"]
        * p["kPumpn"]
        * ATP
        * sodium
        / (1 + ATP / p["KmPump"])
    )

    baseline_pump = pump_rate[0]
    maximum_pump = np.max(pump_rate)
    peak_time = time[np.argmax(pump_rate)]

    # Trapezoidal integration
    cumulative_pump = np.trapezoid(pump_rate, time)
    # Integrated pump activity above baseline
    baseline_auc = baseline_pump * (time[-1] - time[0])

    excess_pump_auc = cumulative_pump - baseline_auc

    print(
        f"Ne={Ne_value}: "
        f"Excess pump activity = {excess_pump_auc:.6f}"
)

    pump_results.append([
        Ne_value,
        baseline_pump,
        maximum_pump,
        peak_time,
        cumulative_pump,
    ])

    print(
        f"Ne={Ne_value}: "
        f"Baseline={baseline_pump:.6f}, "
        f"Maximum={maximum_pump:.6f}, "
        f"Peak time={peak_time:.4f}s, "
        f"AUC={cumulative_pump:.6f}"
    )

    plt.plot(time, pump_rate, label=f"Ne={Ne_value}")

plt.xlabel("Time (s)")
plt.ylabel("Neuronal Na+/K+-ATPase rate")
plt.title("Time-Dependent Neuronal Sodium Pump Activity")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("Pump_time_course.png", dpi=300)
plt.show()

# Save numerical results
np.savetxt(
    "pump_activity_summary.csv",
    np.array(pump_results),
    delimiter=",",
    header="Ne,Baseline_pump,Maximum_pump,Peak_time_s,Pump_AUC",
    comments="",
    fmt="%.8f",
)

print("\nSaved Pump_time_course.png")
print("Saved pump_activity_summary.csv")
