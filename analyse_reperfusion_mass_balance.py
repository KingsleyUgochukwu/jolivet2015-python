
import pandas as pd
import matplotlib.pyplot as plt

durations = [300, 600, 1000, 1500, 2000]

REN = 0.4444444444444445
REG = 0.8

fig, axes = plt.subplots(5, 1, figsize=(12, 18))

for ax, duration in zip(axes, durations):

    df = pd.read_csv(
        f"lactate_flux_timecourse_{duration}s.csv"
    )

    reperfusion_onset = 100 + duration

    # Keep only reperfusion observations
    df = df[df["time_s"] >= reperfusion_onset].copy()

    # Time relative to reperfusion onset
    df["time_after_reperfusion"] = (
        df["time_s"] - reperfusion_onset
    )

    # Individual contributions to extracellular lactate balance
    df["neuronal_contribution"] = df["vLACne"] / REN
    df["astrocytic_contribution"] = df["vLACge"] / REG
    df["vascular_contribution"] = -df["vLACec"]

    ax.plot(
        df["time_after_reperfusion"],
        df["neuronal_contribution"],
        label="Neuronal contribution"
    )

    ax.plot(
        df["time_after_reperfusion"],
        df["astrocytic_contribution"],
        label="Astrocytic contribution"
    )

    ax.plot(
        df["time_after_reperfusion"],
        df["vascular_contribution"],
        label="Vascular contribution"
    )

    ax.axhline(0, linewidth=0.8)
    ax.set_title(f"{duration}s Ischemia")
    ax.set_ylabel("Contribution")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)

axes[-1].set_xlabel("Time after reperfusion (seconds)")

plt.tight_layout()
plt.savefig(
    "reperfusion_lactate_mass_balance.png",
    dpi=300
)

plt.show()

print("Reperfusion mass-balance plots generated.")
