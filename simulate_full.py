import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from parameters import get_parameters
from initial_state import initial_state
from balance_equations import balance_equations


# ============================================================
# FULL JOLIVET ET AL. (2015) SIMULATION
# MATLAB-style staged integration
# ============================================================

# Load parameters and initial conditions
params = get_parameters()
y0 = np.array(initial_state, dtype=float)


# ------------------------------------------------------------
# Experimental protocol
# ------------------------------------------------------------

t0 = 0.0
t1 = 0.0
tstim = 20.0
tend = 3600.0

species = "rat"
CBF = "off"
LONG = "off"

Ne = 240
DCBF = 0.4
F0 = 0.012


# ------------------------------------------------------------
# MATLAB numerical settings
# ------------------------------------------------------------

# During stimulation
short_dt = 1e-4
short_tstep = 0.1

# After stimulation
long_dt = 1.0
long_tstep = 10.0


print("Full Jolivet 2015 simulation")
print("----------------------------")
print("Number of state variables:", len(y0))
print("Start time:", t0, "s")
print("Stimulation duration:", tstim, "s")
print("End time:", tend, "s")
print("Initial neuronal ATP:", y0[14])
print("Initial membrane potential:", y0[27])
# ============================================================
# PHASE 1: STIMULATION
# 0-20 s, integrated in 0.1-s blocks
# ============================================================

print()
print("PHASE 1: STIMULATION")
print("--------------------")

current_time = t0
current_state = y0.copy()

# Store results
time_short = [current_time]
states_short = [current_state.copy()]

block_number = 0


while current_time < (t1 + tstim):

    block_number += 1

    block_end = min(
        current_time + short_tstep,
        t1 + tstim,
    )

    # ODE function for this block
    def model_short(t, y):
        return balance_equations(
            t=t,
            y=y,
            params=params,
            TIME=0.0,
            t1=t1,
            tstim=tstim,
            CBF=CBF,
            F0=F0,
            species=species,
            DCBF=DCBF,
            Ne=Ne,
            LONG=LONG,
        )

    solution = solve_ivp(
        model_short,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=short_dt,
    )

    # Stop immediately if a block fails
    if not solution.success:
        raise RuntimeError(
            f"Solver failed at t={current_time:.4f} s: "
            f"{solution.message}"
        )

    # Save all points except the first one,
    # because it duplicates the previous block endpoint
    time_short.extend(solution.t[1:])
    states_short.extend(solution.y[:, 1:].T)

    # Final state becomes initial state of next block
    current_state = solution.y[:, -1].copy()
    current_time = block_end

    # Progress report every 1 second
    if block_number % 10 == 0:
        print(
            f"t = {current_time:5.1f} s | "
            f"Na = {current_state[0]:.6f} | "
            f"ATP = {current_state[14]:.6f} | "
            f"Vm = {current_state[27]:.3f}"
        )


# Convert stored results to NumPy arrays
time_short = np.array(time_short)
states_short = np.array(states_short)


print()
print("Stimulation phase completed successfully.")
print("Stored time points:", len(time_short))
print("Final stimulation time:", current_time)

print()
print("--- STATE AT END OF STIMULATION ---")
print("Neuronal Na+:", current_state[0])
print("Neuronal ATP:", current_state[14])
print("Membrane potential:", current_state[27])
print("Intracellular calcium:", current_state[30])
print(
    "Neuronal mitochondrial NADH:",
    current_state[31],
)
# ============================================================
# PHASE 2: POST-STIMULATION RECOVERY
# 20-3600 s, integrated in 10-s blocks
# ============================================================

print()
print("PHASE 2: POST-STIMULATION RECOVERY")
print("----------------------------------")

# Carry the final Phase-1 state into Phase 2
time_long = [current_time]
states_long = [current_state.copy()]

recovery_block = 0


def model_long(t, y):
    return balance_equations(
        t=t,
        y=y,
        params=params,
        TIME=0.0,
        t1=t1,
        tstim=tstim,
        CBF=CBF,
        F0=F0,
        species=species,
        DCBF=DCBF,
        Ne=Ne,
        LONG=LONG,
    )


while current_time < tend:

    recovery_block += 1

    block_end = min(
        current_time + long_tstep,
        tend,
    )

    solution = solve_ivp(
        model_long,
        (current_time, block_end),
        current_state,
        method="BDF",
        rtol=1e-3,
        atol=1e-6,
        max_step=long_dt,
    )

    if not solution.success:
        raise RuntimeError(
            f"Recovery solver failed at t={current_time:.1f} s: "
            f"{solution.message}"
        )

    time_long.extend(solution.t[1:])
    states_long.extend(solution.y[:, 1:].T)

    current_state = solution.y[:, -1].copy()
    current_time = block_end

    # Report every 100 seconds
    if recovery_block % 10 == 0 or current_time >= tend:
        print(
            f"t = {current_time:7.1f} s | "
            f"Na = {current_state[0]:.6f} | "
            f"ATP = {current_state[14]:.6f} | "
            f"Vm = {current_state[27]:.3f}"
        )


time_long = np.array(time_long)
states_long = np.array(states_long)


print()
print("Recovery phase completed successfully.")
print("Final simulation time:", current_time)

print()
print("--- FINAL STATE AT 3600 s ---")
print("Neuronal Na+:", current_state[0])
print("Neuronal ATP:", current_state[14])
print("Membrane potential:", current_state[27])
print("Intracellular calcium:", current_state[30])
print(
    "Neuronal mitochondrial NADH:",
    current_state[31],
)
# ============================================================
# COMBINE PHASE 1 AND PHASE 2 RESULTS
# ============================================================

# Avoid duplicating the 20-s boundary
time_full = np.concatenate(
    [time_short, time_long[1:]]
)

states_full = np.vstack(
    [states_short, states_long[1:]]
)

print()
print("Total stored time points:", len(time_full))
print("Combined state array shape:", states_full.shape)


# ============================================================
# FULL-SIMULATION PLOTS
# ============================================================

# ------------------------------------------------------------
# 1. Neuronal ATP
# ------------------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(time_full, states_full[:, 14])

plt.axvline(
    tstim,
    linestyle="--",
    label="End of stimulation",
)

plt.xlabel("Time (s)")
plt.ylabel("Neuronal ATP")
plt.title("Neuronal ATP: Stimulation and Recovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ------------------------------------------------------------
# 2. Neuronal Na+
# ------------------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(time_full, states_full[:, 0])

plt.axvline(
    tstim,
    linestyle="--",
    label="End of stimulation",
)

plt.xlabel("Time (s)")
plt.ylabel("Neuronal Na+")
plt.title("Neuronal Sodium: Stimulation and Recovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ------------------------------------------------------------
# 3. Neuronal mitochondrial NADH
# ------------------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(time_full, states_full[:, 31])

plt.axvline(
    tstim,
    linestyle="--",
    label="End of stimulation",
)

plt.xlabel("Time (s)")
plt.ylabel("Neuronal mitochondrial NADH")
plt.title("Mitochondrial NADH: Stimulation and Recovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ------------------------------------------------------------
# 4. Intracellular Ca2+
# ------------------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(time_full, states_full[:, 30])

plt.axvline(
    tstim,
    linestyle="--",
    label="End of stimulation",
)

plt.xlabel("Time (s)")
plt.ylabel("Intracellular Ca2+")
plt.title("Intracellular Calcium: Stimulation and Recovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()