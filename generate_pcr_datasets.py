"""Generate the NPZ inputs required by the PCr analysis scripts.

Each stimulation condition is simulated once from 0 to 620 s. The first
20 seconds are saved as the stimulation-only dataset; the complete trajectory
is saved as the 600-second post-stimulation recovery dataset.

Run from the repository root:
    python generate_pcr_datasets.py

No model equations or parameter values are changed by this script.
"""
from pathlib import Path
import numpy as np
from simulation_runner import run_simulation

STIMULATIONS = [120, 180, 240, 300, 360]
STIMULATION_DURATION_S = 20.0
RECOVERY_FOLLOWUP_S = 600.0
TOTAL_END_S = STIMULATION_DURATION_S + RECOVERY_FOLLOWUP_S


def save_dataset(path, time, states, Ne, protocol):
    np.savez_compressed(
        Path(path),
        time=np.asarray(time, dtype=float),
        states=np.asarray(states, dtype=float),
        Ne=int(Ne),
        stimulation_duration_s=STIMULATION_DURATION_S,
        recovery_followup_s=RECOVERY_FOLLOWUP_S,
        protocol=protocol,
    )
    print(
        f"Saved {path} | time points={len(time)} | "
        f"states shape={np.asarray(states).shape}"
    )


def main():
    for Ne in STIMULATIONS:
        print(f"\\nGenerating PCr datasets for Ne={Ne}")
        time, states = run_simulation(
            Ne=Ne,
            tend=TOTAL_END_S,
            verbose=True,
        )

        if states.ndim != 2 or states.shape[1] != 33:
            raise RuntimeError(f"Unexpected state shape: {states.shape}")
        if len(time) != len(states):
            raise RuntimeError("Time and state arrays have different lengths.")
        if not np.all(np.isfinite(states)):
            raise RuntimeError(f"Non-finite state values for Ne={Ne}.")

        # The stimulation-only analysis expects the 0-20 s segment.
        stim_mask = time <= STIMULATION_DURATION_S + 1e-9
        save_dataset(
            f"stimulation_Ne_{Ne}.npz",
            time[stim_mask],
            states[stim_mask],
            Ne,
            "0-20 s stimulation segment",
        )

        # Recovery analysis uses time since the original simulation start;
        # its scripts locate the PCr minimum and shift to minimum-relative time.
        save_dataset(
            f"recovery_Ne_{Ne}_600s.npz",
            time,
            states,
            Ne,
            "0-20 s stimulation plus 600 s post-stimulation follow-up",
        )

        pcr = states[:, 16]
        i_min = int(np.argmin(pcr))
        print(
            f"Ne={Ne}: minimum neuronal PCr={pcr[i_min]:.6f} "
            f"at t={time[i_min]:.2f} s"
        )

    print("\\nAll ten PCr input files generated.")
    print("Run analyze_phosphocreatine.py, compare_normalized_pcr_recovery.py,")
    print("and calculate_pcr_rmsd.py to reproduce the PCr analyses.")


if __name__ == "__main__":
    main()
