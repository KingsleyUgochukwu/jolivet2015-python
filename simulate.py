import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from parameters import get_parameters
from initial_state import initial_state
from balance_equations import balance_equations


# -------------------------------------------------
# Load model
# -------------------------------------------------

params = get_parameters()

y0 = np.array(initial_state, dtype=float)


# -------------------------------------------------
# Simulation settings
# -------------------------------------------------

# First stimulated simulation
# Stimulation: 0-20 s
# Recovery: 20-25 s
t_start = 0.0
t_end = 20.0

t_eval = np.linspace(t_start, t_end, 2001)


# -------------------------------------------------
# ODE wrapper
# -------------------------------------------------

def model(t, y):
    return balance_equations(
        t=t,
        y=y,
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="off",
        F0=0.012,
        species="rat",
        DCBF=0.4,
        Ne=240,
        LONG="off",
    )


# -------------------------------------------------
# Solve the 33-state ODE system
# -------------------------------------------------

solution = solve_ivp(
    model,
    (t_start, t_end),
    y0,
    t_eval=t_eval,
    method="BDF",
    rtol=1e-3,
    atol=1e-6,
    max_step=1e-4,
)

# -------------------------------------------------
# Basic diagnostic
# -------------------------------------------------

print("Simulation successful:", solution.success)
print("Solver message:", solution.message)
print("Number of time points:", len(solution.t))
print("Number of state variables:", solution.y.shape[0])

print()
print("--- INITIAL vs FINAL VALUES ---")

print(
    "Neuronal Na+:",
    y0[0],
    "->",
    solution.y[0, -1],
)

print(
    "Neuronal ATP:",
    y0[14],
    "->",
    solution.y[14, -1],
)

print(
    "Membrane potential:",
    y0[27],
    "->",
    solution.y[27, -1],
)

print(
    "Intracellular calcium:",
    y0[30],
    "->",
    solution.y[30, -1],
)

print(
    "Neuronal mitochondrial NADH:",
    y0[31],
    "->",
    solution.y[31, -1],
)
# -------------------------------------------------
# Plot baseline simulation
# -------------------------------------------------

# 1. Membrane potential
plt.figure()
plt.plot(solution.t, solution.y[27])
plt.xlabel("Time")
plt.ylabel("Membrane potential (mV)")
plt.title("Neuronal Membrane Potential")
plt.grid(True)
plt.tight_layout()
plt.show()


# 2. Neuronal ATP
plt.figure()
plt.plot(solution.t, solution.y[14])
plt.xlabel("Time")
plt.ylabel("Neuronal ATP")
plt.title("Neuronal ATP")
plt.grid(True)
plt.tight_layout()
plt.show()


# 3. Neuronal sodium
plt.figure()
plt.plot(solution.t, solution.y[0])
plt.xlabel("Time")
plt.ylabel("Neuronal Na+")
plt.title("Neuronal Sodium")
plt.grid(True)
plt.tight_layout()
plt.show()


# 4. Intracellular calcium
plt.figure()
plt.plot(solution.t, solution.y[30])
plt.xlabel("Time")
plt.ylabel("Intracellular Ca2+")
plt.title("Neuronal Intracellular Calcium")
plt.grid(True)
plt.tight_layout()
plt.show()