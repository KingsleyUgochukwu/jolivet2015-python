import numpy as np

from parameters import get_parameters
from initial_state import initial_state
from rates import calculate_rates
from balance_equations import balance_equations


# ============================================================
# SETTINGS
# ============================================================

CBF = "off"
LONG = "off"

t0 = 0.0
tstim = 20.0
Ne = 240
DCBF = 0.4


# ============================================================
# LOAD MODEL
# ============================================================

params = get_parameters()
y0 = np.array(initial_state, dtype=float)


print()
print("=" * 70)
print("JOLIVET 2015 PYTHON MODEL VALIDATION")
print("=" * 70)

print()
print("--- EXPERIMENT SETTINGS ---")
print("CBF:", CBF)
print("LONG:", LONG)
print("Stimulation duration:", tstim, "s")
print("Ne:", Ne)
print("DCBF:", DCBF)

print()
print("--- MODEL DIMENSIONS ---")
print("Number of state variables:", len(y0))
print("Number of parameters:", len(params))


# ============================================================
# BASIC INITIAL-STATE CHECKS
# ============================================================

print()
print("=" * 70)
print("INITIAL STATE CHECKS")
print("=" * 70)

print("Contains NaN:", np.any(np.isnan(y0)))
print("Contains Inf:", np.any(np.isinf(y0)))
print("All states finite:", np.all(np.isfinite(y0)))

if len(y0) != 33:
    print("WARNING: Model should contain 33 state variables.")
else:
    print("PASS: 33 state variables detected.")


# ============================================================
# STATE NAMES
# ============================================================

state_names = [
    "Neuronal Na+",
    "Glial Na+",
    "Neuronal glucose",
    "Glial glucose",
    "Neuronal GAP",
    "Glial GAP",
    "Neuronal PEP",
    "Glial PEP",
    "Neuronal pyruvate",
    "Glial pyruvate",
    "Neuronal lactate",
    "Glial lactate",
    "Neuronal cytosolic NADH",
    "Glial cytosolic NADH",
    "Neuronal ATP",
    "Glial ATP",
    "Neuronal phosphocreatine",
    "Glial phosphocreatine",
    "Neuronal oxygen",
    "Glial oxygen",
    "Capillary oxygen",
    "Capillary glucose",
    "Capillary lactate",
    "Venous blood volume",
    "Deoxyhemoglobin",
    "Extracellular glucose",
    "Extracellular lactate",
    "Membrane potential",
    "h gating variable",
    "n gating variable",
    "Intracellular calcium",
    "Neuronal mitochondrial NADH",
    "Glial mitochondrial NADH",
]


print()
print("--- ALL 33 INITIAL STATES ---")

for i, (name, value) in enumerate(
    zip(state_names, y0),
    start=1,
):
    print(
        f"Y{i:02d}  "
        f"{name:<32} "
        f"{value:.15g}"
    )


# ============================================================
# CALCULATE INITIAL RATES
# ============================================================

print()
print("=" * 70)
print("INITIAL RATE CALCULATION")
print("=" * 70)

rates = calculate_rates(
    t0,
    y0,
    params,
    t1=0.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)

print("Number of calculated rates:", len(rates))

rate_values = np.array(
    [
        value
        for value in rates.values()
        if np.isscalar(value)
    ],
    dtype=float,
)

print(
    "Rates contain NaN:",
    np.any(np.isnan(rate_values)),
)

print(
    "Rates contain Inf:",
    np.any(np.isinf(rate_values)),
)

print(
    "All scalar rates finite:",
    np.all(np.isfinite(rate_values)),
)


# ============================================================
# PRINT ALL RATES
# ============================================================

print()
print("--- ALL CALCULATED RATES ---")

for name, value in rates.items():

    if np.isscalar(value):
        print(
            f"{name:<20} "
            f"{float(value): .15e}"
        )
    else:
        print(
            f"{name:<20} "
            f"{value}"
        )


# ============================================================
# CALCULATE INITIAL DERIVATIVES
# ============================================================

print()
print("=" * 70)
print("INITIAL DERIVATIVE CHECK")
print("=" * 70)

dydt = balance_equations(
    t=t0,
    y=y0,
    params=params,
    t1=0.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)

print("Number of derivatives:", len(dydt))

print(
    "Derivatives contain NaN:",
    np.any(np.isnan(dydt)),
)

print(
    "Derivatives contain Inf:",
    np.any(np.isinf(dydt)),
)

print(
    "All derivatives finite:",
    np.all(np.isfinite(dydt)),
)


print()
print("--- ALL 33 INITIAL DERIVATIVES ---")

for i, (name, value) in enumerate(
    zip(state_names, dydt),
    start=1,
):
    print(
        f"dY{i:02d} "
        f"{name:<32} "
        f"{value: .15e}"
    )


# ============================================================
# EQUILIBRIUM / RESIDUAL ANALYSIS
# ============================================================

print()
print("=" * 70)
print("INITIAL EQUILIBRIUM RESIDUAL ANALYSIS")
print("=" * 70)

abs_dydt = np.abs(dydt)

max_index = np.argmax(abs_dydt)

print(
    "Largest absolute derivative:",
    abs_dydt[max_index],
)

print(
    "State with largest derivative:",
    f"Y{max_index + 1}",
    state_names[max_index],
)

print(
    "Signed derivative:",
    dydt[max_index],
)


# ============================================================
# CHECK METABOLIC STATES SEPARATELY
# ============================================================

# Exclude membrane voltage because electrophysiology can
# have a different numerical scale.
metabolic_indices = [
    i
    for i in range(33)
    if i != 27
]

metabolic_residuals = abs_dydt[
    metabolic_indices
]

metabolic_max_local = np.argmax(
    metabolic_residuals
)

metabolic_index = metabolic_indices[
    metabolic_max_local
]

print()
print(
    "Largest non-Vm derivative:",
    abs_dydt[metabolic_index],
)

print(
    "Corresponding state:",
    f"Y{metabolic_index + 1}",
    state_names[metabolic_index],
)


# ============================================================
# CBF-OFF CHECK
# ============================================================

print()
print("=" * 70)
print("CBF-OFF VALIDATION")
print("=" * 70)

print(
    "dY21 capillary oxygen:",
    dydt[20],
)

print(
    "dY22 capillary glucose:",
    dydt[21],
)

print(
    "dY23 capillary lactate:",
    dydt[22],
)

print(
    "dY25 deoxyhemoglobin:",
    dydt[24],
)

cbf_off_values = np.array(
    [
        dydt[20],
        dydt[21],
        dydt[22],
        dydt[24],
    ]
)

if np.allclose(
    cbf_off_values,
    0.0,
    atol=1e-12,
):
    print(
        "PASS: CBF-off capillary derivatives "
        "are zero."
    )
else:
    print(
        "WARNING: CBF-off capillary derivatives "
        "are not all zero."
    )


# ============================================================
# IMPORTANT ELECTROPHYSIOLOGY VALUES
# ============================================================

print()
print("=" * 70)
print("ELECTROPHYSIOLOGY CHECK")
print("=" * 70)

important_currents = [
    "IL",
    "INa",
    "IK",
    "ICa",
    "ImAHP",
    "dIPump",
    "Isyne",
    "Isyni",
]

for current in important_currents:

    if current in rates:
        print(
            f"{current:<10}",
            rates[current],
        )


print()
print(
    "Initial membrane potential:",
    y0[27],
)

print(
    "Initial dVm/dt:",
    dydt[27],
)


# ============================================================
# STIMULATION ONSET CHECK
# ============================================================

print()
print("=" * 70)
print("STIMULATION ONSET CHECK")
print("=" * 70)

test_times = [
    0.0,
    0.001,
    0.005,
    0.01,
    0.1,
    1.0,
    20.0,
    20.001,
]

for test_t in test_times:

    test_rates = calculate_rates(
        test_t,
        y0,
        params,
        t1=0.0,
        tstim=tstim,
        CBF=CBF,
        LONG=LONG,
        Ne=Ne,
        DCBF=DCBF,
    )

    print(
        f"t = {test_t:7.3f} s | "
        f"ge = {test_rates['ge']: .8e} | "
        f"Isyne = {test_rates['Isyne']: .8e} | "
        f"vnstim = {test_rates['vnstim']: .8e}"
    )
# ============================================================
# RESTING STATE vs STIMULUS-ONSET VALIDATION
# ============================================================

print()
print("=" * 70)
print("RESTING STATE vs STIMULUS-ONSET VALIDATION")
print("=" * 70)

# ------------------------------------------------------------
# 1. RESTING STATE
# ------------------------------------------------------------
#
# We evaluate the same initial state just BEFORE stimulation.
#
# Setting t1 = 1.0 means stimulation is scheduled to begin
# at 1 second.
#
# At t = 0, therefore, the neuron should still be resting.
# ------------------------------------------------------------

dydt_rest = balance_equations(
    t=0.0,
    y=y0,
    params=params,
    t1=1.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)

rates_rest = calculate_rates(
    0.0,
    y0,
    params,
    t1=1.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)


# ------------------------------------------------------------
# 2. STIMULUS ONSET
# ------------------------------------------------------------
#
# Here stimulation begins at t = 0.
# ------------------------------------------------------------

dydt_stim = balance_equations(
    t=0.0,
    y=y0,
    params=params,
    t1=0.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)

rates_stim = calculate_rates(
    0.0,
    y0,
    params,
    t1=0.0,
    tstim=tstim,
    CBF=CBF,
    LONG=LONG,
    Ne=Ne,
    DCBF=DCBF,
)


# ------------------------------------------------------------
# PRINT COMPARISON
# ------------------------------------------------------------

print()
print("--- RESTING STATE ---")

print("ge:", rates_rest["ge"])
print("Isyne:", rates_rest["Isyne"])
print("vnstim:", rates_rest["vnstim"])

print("dNa_n/dt:", dydt_rest[0])
print("dNa_g/dt:", dydt_rest[1])
print("dATP_n/dt:", dydt_rest[14])
print("dVm/dt:", dydt_rest[27])
print("dCa/dt:", dydt_rest[30])
print("dNADHmito_n/dt:", dydt_rest[31])


print()
print("--- STIMULUS ONSET ---")

print("ge:", rates_stim["ge"])
print("Isyne:", rates_stim["Isyne"])
print("vnstim:", rates_stim["vnstim"])

print("dNa_n/dt:", dydt_stim[0])
print("dNa_g/dt:", dydt_stim[1])
print("dATP_n/dt:", dydt_stim[14])
print("dVm/dt:", dydt_stim[27])
print("dCa/dt:", dydt_stim[30])
print("dNADHmito_n/dt:", dydt_stim[31])


# ------------------------------------------------------------
# RESTING RESIDUAL ANALYSIS
# ------------------------------------------------------------

rest_abs = np.abs(dydt_rest)

rest_max_index = np.argmax(rest_abs)

print()
print("--- RESTING RESIDUAL ANALYSIS ---")

print(
    "Largest resting absolute derivative:",
    rest_abs[rest_max_index],
)

print(
    "Corresponding state:",
    f"Y{rest_max_index + 1}",
    state_names[rest_max_index],
)

print(
    "Signed derivative:",
    dydt_rest[rest_max_index],
)


# ------------------------------------------------------------
# METABOLIC RESTING RESIDUAL
# ------------------------------------------------------------

rest_metabolic_indices = [
    i
    for i in range(33)
    if i != 27
]

rest_metabolic_residuals = rest_abs[
    rest_metabolic_indices
]

rest_metabolic_local = np.argmax(
    rest_metabolic_residuals
)

rest_metabolic_index = rest_metabolic_indices[
    rest_metabolic_local
]

print(
    "Largest resting non-Vm derivative:",
    rest_abs[rest_metabolic_index],
)

print(
    "Corresponding state:",
    f"Y{rest_metabolic_index + 1}",
    state_names[rest_metabolic_index],
)


# ------------------------------------------------------------
# SIMPLE RESTING-STATE TEST
# ------------------------------------------------------------

print()
print("--- RESTING STATE TEST ---")

if (
    abs(rates_rest["ge"]) < 1e-12
    and abs(rates_rest["Isyne"]) < 1e-12
):
    print(
        "PASS: Excitatory stimulation is OFF "
        "before stimulus onset."
    )
else:
    print(
        "WARNING: Excitatory stimulation is not "
        "fully OFF before stimulus onset."
    )


if abs(dydt_rest[27]) < 1e-2:
    print(
        "PASS: Membrane potential is approximately "
        "stationary at rest."
    )
else:
    print(
        "CHECK: Resting membrane potential derivative =",
        dydt_rest[27],
    )


print()
print(
    "Stimulus-induced change in dVm/dt:",
    dydt_stim[27] - dydt_rest[27],
)

print(
    "Stimulus-induced change in neuronal Na derivative:",
    dydt_stim[0] - dydt_rest[0],
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)

if (
    len(y0) == 33
    and len(dydt) == 33
    and np.all(np.isfinite(y0))
    and np.all(np.isfinite(rate_values))
    and np.all(np.isfinite(dydt))
):
    print(
        "PASS: Core numerical integrity checks passed."
    )
else:
    print(
        "FAIL: One or more core numerical checks failed."
    )

print()
print(
    "Validation script completed successfully."
)