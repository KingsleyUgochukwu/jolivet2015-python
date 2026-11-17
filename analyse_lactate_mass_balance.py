
import pandas as pd

durations = [300, 600, 1000, 1500, 2000]

results = []

for duration in durations:

    filename = f"lactate_flux_timecourse_{duration}s.csv"
    df = pd.read_csv(filename)

    reperfusion_onset = 100 + duration

    # Select the reperfusion period
    reperfusion = df[
        df["time_s"] > reperfusion_onset
    ].copy()

    if reperfusion.empty:
        raise ValueError(f"No reperfusion data for {duration}s")

    # Identify the maximum extracellular lactate concentration
    peak_index = reperfusion[
        "extracellular_lactate"
    ].idxmax()

    peak = df.loc[peak_index]

    # Calculate individual contributions to extracellular balance
    neuronal_contribution = peak["vLACne"] / 0.4444444444444445

    astrocytic_contribution = peak["vLACge"] / 0.8

    vascular_removal = -peak["vLACec"]

    net_balance = (
        neuronal_contribution
        + astrocytic_contribution
        + vascular_removal
    )

    results.append({

        "ischemia_duration_s": duration,

        "time_after_reperfusion_s":
            peak["time_s"] - reperfusion_onset,

        "extracellular_lactate":
            peak["extracellular_lactate"],

        "neuronal_contribution":
            neuronal_contribution,

        "astrocytic_contribution":
            astrocytic_contribution,

        "vascular_removal":
            vascular_removal,

        "calculated_net_balance":
            net_balance,

        "model_net_balance":
            peak["net_extracellular_lactate_rate"],

    })


# --------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------

summary = pd.DataFrame(results)

summary.to_csv(
    "lactate_mass_balance_reperfusion.csv",
    index=False,
)

print("\n--- EXTRACELLULAR LACTATE MASS-BALANCE ANALYSIS ---")

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.8f}"
    )
)

print("\nMass-balance analysis completed.")
