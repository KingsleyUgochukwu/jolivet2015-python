
import numpy as np
import warnings
warnings.filterwarnings("error", category=RuntimeWarning)
from scipy.integrate import solve_ivp
from parameters import get_parameters
from initial_state import initial_state
from balance_equations import balance_equations


# --------------------------------------------------
# Ischemia-reperfusion experiment settings
# --------------------------------------------------

params = get_parameters()
y0 = np.array(initial_state, dtype=float).copy()

F0 = 0.012

ischemia_onset = 100.0
ischemia_duration = 1000.0
reperfusion_onset = ischemia_onset + ischemia_duration
ischemia_test_end = reperfusion_onset - 1.0
reperfusion_duration = 1000.0
simulation_end = reperfusion_onset + reperfusion_duration

print("Ischemia-reperfusion simulation setup")
print("Number of state variables:", len(y0))
print("Baseline cerebral blood flow:", F0)
print("Ischemia begins at:", ischemia_onset, "seconds")
print("Reperfusion begins at:", reperfusion_onset, "seconds")
print("Simulation ends at:", simulation_end, "seconds")

# --------------------------------------------------
# Differential-equation model
# --------------------------------------------------

def ischemia_model(t, y):
    oxygen = y[20]

    if not (np.isfinite(oxygen) and 0 < oxygen < params["HbOP"]):
        print("\nINVALID OXYGEN STATE DETECTED")
        print(f"Simulation time: {t:.6f} seconds")
        print(f"Capillary oxygen: {oxygen}")
        print(f"Valid range: 0 < oxygen < {params['HbOP']}")

    return balance_equations(
        t=min(t, reperfusion_onset - 1e-9),
        y=y,
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="on",
        F0=F0,
        species="rat",
        DCBF=0.4,
        Ne=240,
        LONG="off",
        ischemia=True,
        ischemia_duration=ischemia_duration,
    )


print("Ischemia differential-equation function defined.")

# --------------------------------------------------
# Phase 1: Initial stimulation (0–20 seconds)
# --------------------------------------------------

current_time = 0.0
current_state = y0.copy()

while current_time < 20.0:

    block_end = min(current_time + 0.1, 20.0)

    solution = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=1e-4,
    )

    if not solution.success:
        raise RuntimeError(
            f"Solver failed at {current_time:.4f}s: "
            f"{solution.message}"
        )

    current_state = solution.y[:, -1].copy()
    current_time = block_end

print("Initial stimulation phase completed.")
print("Simulation time:", current_time)
print("Neuronal ATP:", current_state[14])
print("Astrocytic ATP:", current_state[15])

# --------------------------------------------------
# Phase 2: Baseline recovery (20–100 seconds)
# --------------------------------------------------

while current_time < ischemia_onset:

    block_end = min(
        current_time + 10.0,
        ischemia_onset,
    )

    solution = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=1.0,
    )

    if not solution.success:
        raise RuntimeError(
            f"Baseline solver failed at {current_time:.2f}s: "
            f"{solution.message}"
        )

    current_state = solution.y[:, -1].copy()
    current_time = block_end

print("\nBaseline phase completed.")
print("Simulation time:", current_time)
print("Neuronal ATP before ischemia:", current_state[14])
print("Astrocytic ATP before ischemia:", current_state[15])
print("Extracellular glucose:", current_state[25])
print("Extracellular lactate:", current_state[26])

# --------------------------------------------------
# Phase 3: Cerebral ischemia (duration set by ischemia_duration)
# --------------------------------------------------

# Save pre-ischemia reference values
baseline_state = current_state.copy()

# Collect results throughout ischemia
ischemia_times = [current_time]
ischemia_states = [current_state.copy()]

while current_time < reperfusion_onset:

    block_end = min(current_time + 10.0, reperfusion_onset)

    solution = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-5,
        atol=1e-6,
        max_step=1.0,
    )

    if not solution.success:
        raise RuntimeError(
            f"Ischemia solver failed at {current_time:.2f}s: "
            f"{solution.message}"
        )

    ischemia_times.extend(solution.t[1:])
    ischemia_states.extend(solution.y[:, 1:].T)

    current_state = solution.y[:, -1].copy()
    current_time = block_end

print("\nIschemia phase completed.")
print("Simulation time:", current_time)

print("Neuronal ATP:", current_state[14])
print("Astrocytic ATP:", current_state[15])
print("Extracellular glucose:", current_state[25])
print("Extracellular lactate:", current_state[26])

print("\nATP changes relative to pre-ischemia baseline:")

neuron_change = (
    (current_state[14] - baseline_state[14])
    / baseline_state[14] * 100
)

astrocyte_change = (
    (current_state[15] - baseline_state[15])
    / baseline_state[15] * 100
)

print(f"Neuronal ATP change: {neuron_change:.2f}%")
print(f"Astrocytic ATP change: {astrocyte_change:.2f}%")
end_ischemia_state = current_state.copy()

# --------------------------------------------------
# Matched control: normal blood flow during the corresponding period
# --------------------------------------------------

def control_model(t, y):
    return balance_equations(
        t=t,
        y=y,
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="on",
        F0=F0,
        species="rat",
        DCBF=0.4,
        Ne=240,
        LONG="off",
        ischemia=False,
    )


control_time = ischemia_onset
control_state = baseline_state.copy()

while control_time < reperfusion_onset:

    block_end = min(control_time + 10.0, reperfusion_onset)

    solution = solve_ivp(
        control_model,
        (control_time, block_end),
        control_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=1.0,
    )

    if not solution.success:
        raise RuntimeError(
            f"Control solver failed at {control_time:.2f}s: "
            f"{solution.message}"
        )

    control_state = solution.y[:, -1].copy()
    control_time = block_end


print(f"\n--- CONTROL VS ISCHEMIA AT {reperfusion_onset:.0f} s ---")

labels = {
    14: "Neuronal ATP",
    15: "Astrocytic ATP",
    25: "Extracellular glucose",
    26: "Extracellular lactate",
}

for index, label in labels.items():
    control_value = control_state[index]
    ischemia_value = end_ischemia_state[index]

    difference = (
        (ischemia_value - control_value)
        / control_value * 100
    )

    print(f"\n{label}")
    print(f"Control: {control_value:.6f}")
    print(f"Ischemia: {ischemia_value:.6f}")
    print(f"Difference: {difference:+.2f}%")

# --------------------------------------------------
# Diagnostic: capillary oxygen during ischemia
# --------------------------------------------------

oxygen_values = np.array(ischemia_states)[:, 20]
oxygen_times = np.array(ischemia_times)

print("\n--- CAPILLARY OXYGEN DIAGNOSTIC ---")
print("Minimum oxygen:", np.min(oxygen_values))
print("Maximum oxygen:", np.max(oxygen_values))
print("HbOP upper limit:", params["HbOP"])

invalid = (
    (oxygen_values <= 0)
    | (oxygen_values >= params["HbOP"])
    | (~np.isfinite(oxygen_values))
)

print("Number of invalid recorded oxygen states:", np.sum(invalid))

if np.any(invalid):
    first_invalid = np.where(invalid)[0][0]
    print("First invalid recorded time:", oxygen_times[first_invalid])
    print("Oxygen at that time:", oxygen_values[first_invalid])
else:
    print("All recorded oxygen states are within the valid range.")

# --------------------------------------------------
# Phase 4: Reperfusion model
# --------------------------------------------------

def reperfusion_model(t, y):
    return balance_equations(
        t=t,
        y=y,
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="on",
        F0=F0,
        species="rat",
        DCBF=0.4,
        Ne=240,
        LONG="off",
        ischemia=False,
    )

print("\nReperfusion differential-equation function defined.")

# --------------------------------------------------
# Phase 5: Reperfusion (duration set by reperfusion_duration)
# --------------------------------------------------

reperfusion_start = reperfusion_onset
reperfusion_test_end = simulation_end

reperfusion_initial_state = end_ischemia_state.copy()

print("\n--- INITIAL REPERFUSION TEST ---")
print("Start time:", reperfusion_start)
print("Initial neuronal ATP:", reperfusion_initial_state[14])
print("Initial astrocytic ATP:", reperfusion_initial_state[15])
print("Initial capillary oxygen:", reperfusion_initial_state[20])

solution = solve_ivp(
    reperfusion_model,
    (reperfusion_start, reperfusion_test_end),
    reperfusion_initial_state,
    method="BDF",
    rtol=1e-5,
    atol=1e-6,
    max_step=0.05,
)

if not solution.success:
    raise RuntimeError(
        f"Reperfusion solver failed: {solution.message}"
    )

reperfusion_state = solution.y[:, -1].copy()

print("\nReperfusion test completed.")
print("Simulation time:", solution.t[-1])
print("Neuronal ATP:", reperfusion_state[14])
print("Astrocytic ATP:", reperfusion_state[15])
print("Extracellular glucose:", reperfusion_state[25])
print("Extracellular lactate:", reperfusion_state[26])
print("Capillary oxygen:", reperfusion_state[20])

oxygen_trace = solution.y[20, :]
oxygen_valid = np.all(
    np.isfinite(oxygen_trace)
    & (oxygen_trace > 0)
    & (oxygen_trace < params["HbOP"])
)

print("All recorded reperfusion oxygen states valid:", oxygen_valid)

# --------------------------------------------------
# ATP recovery relative to pre-ischemia baseline
# --------------------------------------------------

print("\n--- ATP RECOVERY ANALYSIS ---")

for index, label in [(14, "Neuronal ATP"), (15, "Astrocytic ATP")]:

    baseline = baseline_state[index]
    ischemia_end = end_ischemia_state[index]
    reperfusion_end = reperfusion_state[index]

    baseline_recovery = (reperfusion_end / baseline) * 100

    depletion_recovered = (
        (reperfusion_end - ischemia_end)
        / (baseline - ischemia_end)
        * 100
    )

    print(f"\n{label}")
    print(f"Pre-ischemia baseline: {baseline:.6f}")
    print(f"End of ischemia: {ischemia_end:.6f}")
    print(f"End of reperfusion: {reperfusion_end:.6f}")
    print(f"ATP relative to baseline: {baseline_recovery:.2f}%")
    print(f"Initial ATP depletion recovered: {depletion_recovered:.2f}%")

# --------------------------------------------------
# Save ischemia and reperfusion time-course data
# --------------------------------------------------

ischemia_array = np.array(ischemia_states)

ischemia_data = np.column_stack([
    np.array(ischemia_times),
    ischemia_array[:, 14],
    ischemia_array[:, 15],
    ischemia_array[:, 25],
    ischemia_array[:, 26],
    ischemia_array[:, 20],
])

reperfusion_data = np.column_stack([
    solution.t,
    solution.y[14, :],
    solution.y[15, :],
    solution.y[25, :],
    solution.y[26, :],
    solution.y[20, :],
])

# Avoid duplicating the shared reperfusion boundary
combined_data = np.vstack([
    ischemia_data,
    reperfusion_data[1:, :],
])

np.savetxt(
    "ischemia_reperfusion_timeseries.csv",
    combined_data,
    delimiter=",",
    header="time_s,neuronal_ATP,astrocytic_ATP,extracellular_glucose,extracellular_lactate,capillary_oxygen",
    comments="",
)

print("\nTime-course dataset saved successfully.")
print("Number of recorded time points:", len(combined_data))
print("First recorded time:", combined_data[0, 0])
print("Last recorded time:", combined_data[-1, 0])



