
import numpy as np
import pandas as pd
import sys

from parameters import get_parameters
from rates import calculate_rates
from initial_state import initial_state
from balance_equations import balance_equations
from scipy.integrate import solve_ivp
from ischemia_protocol import ischemia_blood_flow


# --------------------------------------------------
# 1. SIMULATION SETTINGS
# --------------------------------------------------

params = get_parameters()
y0 = np.array(initial_state, dtype=float).copy()

F0 = 0.012
ischemia_onset = 100.0
ischemia_duration = (
    float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0
)
reperfusion_onset = ischemia_onset + ischemia_duration
reperfusion_duration = 1000.0
simulation_end = reperfusion_onset + reperfusion_duration

print("LACTATE FLUX DIAGNOSTIC")
print("Ischemia duration:", ischemia_duration, "seconds")


# --------------------------------------------------
# 2. MODEL FUNCTIONS
# --------------------------------------------------

def ischemia_model(t, y):
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
        ischemia=True,
        ischemia_duration=ischemia_duration,
    )


def normal_flow_model(t, y):
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


# --------------------------------------------------
# 3. REPRODUCE PRE-ISCHEMIA STATE
# --------------------------------------------------

current_time = 0.0
current_state = y0.copy()

# Initial stimulation: 0–20 seconds
while current_time < 20.0:
    block_end = min(current_time + 0.1, 20.0)

    sol = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=1e-4,
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    current_state = sol.y[:, -1].copy()
    current_time = block_end

# Baseline recovery: 20–100 seconds
while current_time < ischemia_onset:
    block_end = min(current_time + 10.0, ischemia_onset)

    sol = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=1.0,
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    current_state = sol.y[:, -1].copy()
    current_time = block_end

baseline_state = current_state.copy()

print("Pre-ischemia baseline reproduced.")


# --------------------------------------------------
# 4. SIMULATE ISCHEMIA
# --------------------------------------------------

ischemia_times = [current_time]
ischemia_states = [current_state.copy()]

while current_time < reperfusion_onset:
    block_end = min(current_time + 10.0, reperfusion_onset)

    sol = solve_ivp(
        ischemia_model,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-5,
        atol=1e-6,
        max_step=1.0,
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    ischemia_times.extend(sol.t[1:])
    ischemia_states.extend(sol.y[:, 1:].T)

    current_state = sol.y[:, -1].copy()
    current_time = block_end

end_ischemia_state = current_state.copy()

print("Ischemia simulation completed.")


# --------------------------------------------------
# 5. SIMULATE REPERFUSION
# --------------------------------------------------

reperfusion_solution = solve_ivp(
    normal_flow_model,
    (reperfusion_onset, simulation_end),
    end_ischemia_state,
    method="BDF",
    rtol=1e-5,
    atol=1e-6,
    max_step=0.05,
)

if not reperfusion_solution.success:
    raise RuntimeError(reperfusion_solution.message)

print("Reperfusion simulation completed.")


# --------------------------------------------------
# 6. COMBINE FULL STATE TRAJECTORIES
# --------------------------------------------------

all_times = np.concatenate([
    np.array(ischemia_times),
    reperfusion_solution.t[1:],
])

all_states = np.vstack([
    np.array(ischemia_states),
    reperfusion_solution.y[:, 1:].T,
])


# --------------------------------------------------
# 7. CALCULATE LACTATE FLUXES
# --------------------------------------------------

def get_lactate_fluxes(t, y, ischemia_active):

    rates = calculate_rates(
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
        ischemia=ischemia_active,
        ischemia_duration=ischemia_duration,
    )

    net_extracellular = (
        rates["vLACne"] / params["ren"]
        + rates["vLACge"] / params["reg"]
        - rates["vLACec"]
    )

    return {
        "time_s": t,
        "neuronal_lactate": y[10],
        "astrocytic_lactate": y[11],
        "extracellular_lactate": y[26],
        "capillary_lactate": y[22],
        "vLDHn": rates["vLDHn"],
        "vLDHg": rates["vLDHg"],
        "vLACne": rates["vLACne"],
        "vLACge": rates["vLACge"],
        "vLACgc": rates["vLACgc"],
        "vLACec": rates["vLACec"],
        "vLACc": rates["vLACc"],
        "net_extracellular_lactate_rate": net_extracellular,
    }


records = []

for t, y in zip(all_times, all_states):

    ischemia_active = (
        ischemia_onset <= t <= reperfusion_onset
    )

    # Evaluate the end-ischemia flux immediately before
    # the blood-flow restoration boundary.
    rate_time = (
        min(t, reperfusion_onset - 1e-9)
        if ischemia_active
        else t
    )

    records.append(
        get_lactate_fluxes(rate_time, y, ischemia_active)
    )

flux_data = pd.DataFrame(records)

flux_data.to_csv(
    f"lactate_flux_timecourse_{int(ischemia_duration)}s.csv",
    index=False,
)

print("Full lactate flux dataset saved.")


# --------------------------------------------------
# 8. IDENTIFY IMPORTANT SAMPLING POINTS
# --------------------------------------------------

lactate_peak_index = np.argmax(
    flux_data["extracellular_lactate"].values[
        (flux_data["time_s"] >= ischemia_onset)
        & (flux_data["time_s"] <= reperfusion_onset)
    ]
)

ischemia_mask = (
    (flux_data["time_s"] >= ischemia_onset)
    & (flux_data["time_s"] <= reperfusion_onset)
)

ischemia_data = flux_data.loc[ischemia_mask].reset_index(drop=True)

peak_row = ischemia_data.iloc[lactate_peak_index]

sampling_points = {
    "Baseline": flux_data.iloc[
        np.argmin(np.abs(flux_data["time_s"] - ischemia_onset))
    ],
    "Early lactate peak": peak_row,
    "End of ischemia": flux_data.iloc[
        np.argmin(np.abs(flux_data["time_s"] - reperfusion_onset))
    ],
    "End of reperfusion": flux_data.iloc[
        np.argmin(np.abs(flux_data["time_s"] - simulation_end))
    ],
}


# --------------------------------------------------
# 9. PRINT COMPARATIVE RESULTS
# --------------------------------------------------

selected_columns = [
    "vLDHn",
    "vLDHg",
    "vLACne",
    "vLACge",
    "vLACgc",
    "vLACec",
    "vLACc",
    "net_extracellular_lactate_rate",
]

print("\n--- LACTATE FLUX COMPARISON ---")

for stage, row in sampling_points.items():

    print(f"\n{stage}")
    print(f"Time: {row['time_s']:.2f} seconds")
    print(
        f"Extracellular lactate: "
        f"{row['extracellular_lactate']:.6f}"
    )

    for column in selected_columns:
        print(f"{column}: {row[column]:.8f}")


# --------------------------------------------------
# 10. SAVE SELECTED SAMPLING POINTS
# --------------------------------------------------

summary = pd.DataFrame.from_dict(
    sampling_points,
    orient="index",
)

summary.index.name = "stage"
summary.to_csv(
    f"lactate_flux_stage_summary_{int(ischemia_duration)}s.csv"
)

print("\nDiagnostic completed successfully.")
