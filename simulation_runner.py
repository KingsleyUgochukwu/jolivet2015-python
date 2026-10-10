
import numpy as np
from scipy.integrate import solve_ivp

from parameters import get_parameters
from initial_state import initial_state
from balance_equations import balance_equations


def run_simulation(Ne=240, tend=3600.0, verbose=True):
    """
    Run the Jolivet et al. (2015) model using the
    validated two-phase numerical integration protocol.

    Parameters
    ----------
    Ne : int
        Neuronal stimulation parameter.
    tend : float
        Final simulation time in seconds (must be >= 20).
    verbose : bool
        Whether to print progress messages.

    Returns
    -------
    time_full : numpy.ndarray
        Simulation time points.
    states_full : numpy.ndarray
        State variables, shape (time points, 33).
    """

    if tend < 20.0:
        raise ValueError("tend must be at least 20 seconds.")

    params = get_parameters()
    y0 = np.array(initial_state, dtype=float).copy()

    t0 = 0.0
    t1 = 0.0
    tstim = 20.0

    species = "rat"
    CBF = "off"
    LONG = "off"

    DCBF = 0.4
    F0 = 0.012

    short_dt = 1e-4
    short_tstep = 0.1

    long_dt = 1.0
    long_tstep = 10.0

    if verbose:
        print(f"Running simulation: Ne={Ne}, tend={tend}s")

    def model(t, y):
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

    # PHASE 1: STIMULATION (0-20 seconds)
    current_time = t0
    current_state = y0.copy()

    time_short = [current_time]
    states_short = [current_state.copy()]

    while current_time < tstim:

        block_end = min(
            current_time + short_tstep,
            tstim,
        )

        solution = solve_ivp(
            model,
            (current_time, block_end),
            current_state,
            method="BDF",
            rtol=1e-3,
            atol=1e-6,
            max_step=short_dt,
        )

        if not solution.success:
            raise RuntimeError(
                f"Stimulation solver failed at "
                f"t={current_time:.4f}s: {solution.message}"
            )

        time_short.extend(solution.t[1:])
        states_short.extend(solution.y[:, 1:].T)

        current_state = solution.y[:, -1].copy()
        current_time = block_end

    # PHASE 2: RECOVERY (20 seconds to tend)
    time_long = [current_time]
    states_long = [current_state.copy()]

    while current_time < tend:

        block_end = min(
            current_time + long_tstep,
            tend,
        )

        solution = solve_ivp(
            model,
            (current_time, block_end),
            current_state,
            method="BDF",
            rtol=1e-3,
            atol=1e-6,
            max_step=long_dt,
        )

        if not solution.success:
            raise RuntimeError(
                f"Recovery solver failed at "
                f"t={current_time:.4f}s: {solution.message}"
            )

        time_long.extend(solution.t[1:])
        states_long.extend(solution.y[:, 1:].T)

        current_state = solution.y[:, -1].copy()
        current_time = block_end

    # ====================================================
    # COMBINE STIMULATION AND RECOVERY RESULTS
    # ====================================================

    if len(time_long) > 1:

        time_full = np.concatenate([
            np.asarray(time_short),
            np.asarray(time_long[1:]),
        ])

        states_full = np.vstack([
            np.asarray(states_short),
            np.asarray(states_long[1:]),
        ])

    else:
        # No recovery phase when tend = 20 seconds
        time_full = np.asarray(time_short)
        states_full = np.asarray(states_short)

    # ====================================================
    # DISPLAY RESULTS
    # ====================================================

    if verbose:
        print("Simulation completed.")
        print("Final time:", time_full[-1])
        print("State array shape:", states_full.shape)
        print("Final ATP:", states_full[-1, 14])

    return time_full, states_full