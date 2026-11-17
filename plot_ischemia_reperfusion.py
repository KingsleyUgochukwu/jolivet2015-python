
import pandas as pd
import matplotlib.pyplot as plt

# Load simulation results
data = pd.read_csv("ischemia_reperfusion_timeseries.csv")

# Create figure
plt.figure(figsize=(11, 6))

# Plot ATP concentrations
plt.plot(
    data["time_s"],
    data["neuronal_ATP"],
    label="Neuronal ATP",
    linewidth=2
)

plt.plot(
    data["time_s"],
    data["astrocytic_ATP"],
    label="Astrocytic ATP",
    linewidth=2
)

# Indicate reperfusion onset
plt.axvline(
    x=1100,
    color="black",
    linestyle="--",
    label="Reperfusion onset"
)

# Highlight ischemia period
plt.axvspan(
    100,
    1100,
    color="gray",
    alpha=0.15,
    label="Ischemia"
)

plt.xlabel("Simulation time (s)")
plt.ylabel("ATP concentration (model units)")
plt.title("Neuronal and Astrocytic ATP During Ischemia–Reperfusion")

plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

# Save high-resolution figure
plt.savefig(
    "ATP_ischemia_reperfusion.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# --------------------------------------------------
# Calculate ATP recovery times
# --------------------------------------------------

reperfusion = data[data["time_s"] >= 1100].copy()

baseline_neuron = 2.197280421124165
baseline_astrocyte = 2.150625724452968

print("\n--- ATP RECOVERY TIMES ---")

for column, baseline, label in [
    ("neuronal_ATP", baseline_neuron, "Neuronal ATP"),
    ("astrocytic_ATP", baseline_astrocyte, "Astrocytic ATP"),
]:
    print(f"\n{label}")

    for percentage in [50, 90]:
        threshold = baseline * percentage / 100

        recovered = reperfusion[
            reperfusion[column] >= threshold
        ]

        if recovered.empty:
            print(f"{percentage}% recovery: Not reached")
        else:
            recovery_time = (
                recovered["time_s"].iloc[0] - 1100
            )

            print(
                f"{percentage}% recovery: "
                f"{recovery_time:.2f} seconds after reperfusion"
            )

# --------------------------------------------------
# Calculate maximum ATP recovery rates
# --------------------------------------------------

import numpy as np

print("\n--- MAXIMUM ATP RECOVERY RATES ---")

reperfusion_time = reperfusion["time_s"].to_numpy()

for column, label in [
    ("neuronal_ATP", "Neuronal ATP"),
    ("astrocytic_ATP", "Astrocytic ATP"),
]:
    atp = reperfusion[column].to_numpy()

    # Numerical derivative of ATP concentration
    recovery_rate = np.gradient(atp, reperfusion_time)

    maximum_index = np.argmax(recovery_rate)

    maximum_rate = recovery_rate[maximum_index]
    maximum_time = reperfusion_time[maximum_index]

    print(f"\n{label}")
    print(f"Maximum recovery rate: {maximum_rate:.6f} model units/s")
    print(f"Time of maximum recovery: {maximum_time:.2f} s")
    print(
        "Time after reperfusion onset: "
        f"{maximum_time - 1100:.2f} s"
    )

# --------------------------------------------------
# Uniform-grid derivative check
# --------------------------------------------------

print("\n--- UNIFORM-GRID RECOVERY RATE CHECK ---")

uniform_time = np.arange(1100.0, 2100.0 + 0.001, 0.5)

for column, label in [
    ("neuronal_ATP", "Neuronal ATP"),
    ("astrocytic_ATP", "Astrocytic ATP"),
]:
    original_time = reperfusion["time_s"].to_numpy()
    original_atp = reperfusion[column].to_numpy()

    uniform_atp = np.interp(
        uniform_time,
        original_time,
        original_atp
    )

    uniform_rate = np.gradient(
        uniform_atp,
        uniform_time
    )

    peak_index = np.argmax(uniform_rate)

    print(f"\n{label}")
    print(f"Maximum rate: {uniform_rate[peak_index]:.6f} model units/s")
    print(f"Peak time: {uniform_time[peak_index]:.2f} s")
    print(
        f"Time after reperfusion: "
        f"{uniform_time[peak_index] - 1100:.2f} s"
    )

# --------------------------------------------------
# Step 60: Solver tolerance comparison
# --------------------------------------------------

original = pd.read_csv("ischemia_reperfusion_original.csv")
tighter = pd.read_csv("ischemia_reperfusion_timeseries.csv")

comparison_time = np.arange(1100.0, 2100.01, 0.5)

print("\n--- SOLVER TOLERANCE COMPARISON ---")

for column, label in [
    ("neuronal_ATP", "Neuronal ATP"),
    ("astrocytic_ATP", "Astrocytic ATP"),
]:
    original_values = np.interp(
        comparison_time,
        original["time_s"],
        original[column]
    )

    tighter_values = np.interp(
        comparison_time,
        tighter["time_s"],
        tighter[column]
    )

    absolute_difference = np.abs(
        original_values - tighter_values
    )

    maximum_difference = np.max(absolute_difference)
    mean_difference = np.mean(absolute_difference)

    print(f"\n{label}")
    print(f"Maximum absolute difference: {maximum_difference:.8f}")
    print(f"Mean absolute difference: {mean_difference:.8f}")

# --------------------------------------------------
# Maximum step-size comparison
# --------------------------------------------------

step_01 = pd.read_csv(
    "ischemia_reperfusion_rtol_1e-5.csv"
)

step_005 = pd.read_csv(
    "ischemia_reperfusion_timeseries.csv"
)

comparison_time = np.arange(1100.0, 2100.01, 0.5)

print("\n--- MAXIMUM STEP-SIZE COMPARISON ---")

for column, label in [
    ("neuronal_ATP", "Neuronal ATP"),
    ("astrocytic_ATP", "Astrocytic ATP"),
]:
    values_01 = np.interp(
        comparison_time,
        step_01["time_s"],
        step_01[column]
    )

    values_005 = np.interp(
        comparison_time,
        step_005["time_s"],
        step_005[column]
    )

    difference = np.abs(values_01 - values_005)

    print(f"\n{label}")
    print(
        f"Maximum absolute difference: "
        f"{np.max(difference):.8f}"
    )
    print(
        f"Mean absolute difference: "
        f"{np.mean(difference):.8f}"
    )

# --------------------------------------------------
# Step 66: Ischemia solver tolerance comparison
# --------------------------------------------------

previous_ischemia = pd.read_csv(
    "ischemia_reperfusion_validated_step005.csv"
)

new_ischemia = pd.read_csv(
    "ischemia_reperfusion_timeseries.csv"
)

ischemia_comparison_time = np.arange(
    100.0, 1100.01, 0.5
)

print("\n--- ISCHEMIA SOLVER TOLERANCE COMPARISON ---")

for column, label in [
    ("neuronal_ATP", "Neuronal ATP"),
    ("astrocytic_ATP", "Astrocytic ATP"),
]:
    previous_values = np.interp(
        ischemia_comparison_time,
        previous_ischemia["time_s"],
        previous_ischemia[column]
    )

    new_values = np.interp(
        ischemia_comparison_time,
        new_ischemia["time_s"],
        new_ischemia[column]
    )

    difference = np.abs(previous_values - new_values)

    print(f"\n{label}")
    print(
        f"Maximum absolute difference: "
        f"{np.max(difference):.8f}"
    )
    print(
        f"Mean absolute difference: "
        f"{np.mean(difference):.8f}"
    )
