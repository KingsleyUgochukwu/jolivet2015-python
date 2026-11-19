"""Generate the NPZ inputs required by the PCr analysis scripts.

This uses the existing simulation_runner.run_simulation() implementation.
It does not modify model equations or parameter values.

Outputs (for each Ne in 120, 180, 240, 300, 360):
  stimulation_Ne_{Ne}.npz       : simulation from 0 to 20 s
  recovery_Ne_{Ne}_600s.npz     : simulation from 0 to 620 s
                                  (20 s stimulation + 600 s follow-up)

Run from the repository root:
  python generate_pcr_datasets.py

The simulations may take some time because the model is stiff and the
runner uses small integration steps during stimulation.
"""
from pathlib import Path
import numpy as np
from simulation_runner import run_simulation

STIMULATIONS = [120, 180, 240, 300, 360]
STIMULATION_DURATION_S = 20.0
RECOVERY_FOLLOWUP_S = 600.0
OUTPUT_DIR = Path(".")

def save_dataset(path, time, states, Ne, protocol):
    path = Path(path)
    np.savez_compressed(
        path,
        time=np.asarray(time, dtype=float),
        states=np.asarray(states, dtype=float),
        Ne=int(Ne),
        stimulation_duration_s=STIMULATION_DURATION_S,
        recovery_followup_s=RECOVERY_FOLLOWUP_S,
        protocol=protocol,
    )
    print(
        f"Saved {path} | "
        f"time points={len(time)} | "
        f"states shape={np.asarray(states).shape}"
    )

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for Ne in STIMULATIONS:
        print(f"\nGenerating PCr datasets for Ne={Ne}")

        # Stimulation-only output used by analyze_phosphocreatine.py
        time_stim, states_stim = run_simulation(
            Ne=Ne,
            tend=STIMULATION_DURATION_S,
            verbose=True,
        )
        save_dataset(
            OUTPUT_DIR / f"stimulation_Ne_{Ne}.npz",
            time_stim,
            states_stim,
            Ne,
            "20-second stimulation",
        )

        # Full 20 s stimulation plus 600 s post-stimulation follow-up.
        time_recovery, states_recovery = run_simulation(
            Ne=Ne,
            tend=STIMULATION_DURATION_S + RECOVERY_FOLLOWUP_S,
            verbose=True,
        )
        save_dataset(
            OUTPUT_DIR / f"recovery_Ne_{Ne}_600s.npz",
            time_recovery,
            states_recovery,
            Ne,
            "20-second stimulation plus 600-second follow-up",
        )

    print("\nAll ten PCr analysis input files have been generated.")
    print("Next run:")
    print("  python analyze_phosphocreatine.py")
    print("  python compare_normalized_pcr_recovery.py")
    print("  python calculate_pcr_rmsd.py")

if __name__ == "__main__":
    main()
