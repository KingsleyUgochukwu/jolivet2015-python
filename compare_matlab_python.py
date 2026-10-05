import numpy as np

from parameters import get_parameters
from initial_state import initial_state
from rates import calculate_rates
from balance_equations import balance_equations


# ============================================================
# JOLIVET 2015 MATLAB <-> PYTHON TRANSLATION AUDIT
# ============================================================

params = get_parameters()
y0 = np.asarray(initial_state, dtype=float)

CBF = "off"
LONG = "off"

TSTIM = 20.0
NE = 240
DCBF = 0.4


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def print_heading(title):
    print()
    print("=" * 75)
    print(title)
    print("=" * 75)


def relative_error(reference, test):
    """
    Relative percentage error.

    Returns NaN when the reference value is effectively zero,
    because percentage error is not meaningful there.
    """
    if abs(reference) < 1e-14:
        return np.nan

    return (
        abs(test - reference)
        / abs(reference)
        * 100.0
    )


# ============================================================
# 1. BASIC MODEL INFORMATION
# ============================================================

print_heading(
    "JOLIVET 2015 MATLAB <-> PYTHON TRANSLATION AUDIT"
)

print("State variables:", len(y0))
print("Parameters:", len(params))
print("CBF:", CBF)
print("LONG:", LONG)
print("Stimulation duration:", TSTIM, "s")
print("Ne:", NE)
print("DCBF:", DCBF)


# ============================================================
# 2. INITIAL STATE REFERENCE VALUES
# ============================================================

state_names = [
    "Na_n",
    "Na_g",
    "GLC_n",
    "GLC_g",
    "GAP_n",
    "GAP_g",
    "PEP_n",
    "PEP_g",
    "PYR_n",
    "PYR_g",
    "LAC_n",
    "LAC_g",
    "NADHcyt_n",
    "NADHcyt_g",
    "ATP_n",
    "ATP_g",
    "PCr_n",
    "PCr_g",
    "O2_n",
    "O2_g",
    "O2_cap",
    "GLC_cap",
    "LAC_cap",
    "Vv",
    "dHb",
    "GLC_e",
    "LAC_e",
    "Vm",
    "h",
    "n",
    "Ca",
    "NADHmito_n",
    "NADHmito_g",
]


print_heading("INITIAL STATE VECTOR")

for index, (name, value) in enumerate(
    zip(state_names, y0),
    start=1,
):
    print(
        f"Y{index:02d} "
        f"{name:<14} "
        f"{value: .16e}"
    )


# ============================================================
# 3. RESTING CONDITION
# ============================================================
#
# Stimulation is scheduled to begin at t = 1 s.
# We evaluate the model at t = 0 s.
#
# Therefore:
#
#       t < t1
#
# and excitatory stimulation should be OFF.
# ============================================================

rest_rates = calculate_rates(
    0.0,
    y0,
    params,
    t1=1.0,
    tstim=TSTIM,
    CBF=CBF,
    LONG=LONG,
    Ne=NE,
    DCBF=DCBF,
)

rest_dydt = balance_equations(
    t=0.0,
    y=y0,
    params=params,
    t1=1.0,
    tstim=TSTIM,
    CBF=CBF,
    LONG=LONG,
    Ne=NE,
    DCBF=DCBF,
)


# ============================================================
# 4. STIMULUS-ONSET CONDITION
# ============================================================
#
# Stimulation begins at t = 0.
# ============================================================

stim_rates = calculate_rates(
    0.0,
    y0,
    params,
    t1=0.0,
    tstim=TSTIM,
    CBF=CBF,
    LONG=LONG,
    Ne=NE,
    DCBF=DCBF,
)

stim_dydt = balance_equations(
    t=0.0,
    y=y0,
    params=params,
    t1=0.0,
    tstim=TSTIM,
    CBF=CBF,
    LONG=LONG,
    Ne=NE,
    DCBF=DCBF,
)


# ============================================================
# 5. IMPORTANT RATE FAMILIES
# ============================================================

important_rates = [

    # Na/K transport
    "vLeakNan",
    "vLeakNag",
    "vPumpn",
    "vPumpg",

    # Glucose / glycolysis
    "vGLCen",
    "vGLCeg",
    "vGLCcg",
    "vHKPFKn",
    "vHKPFKg",
    "vPGKn",
    "vPGKg",
    "vPKn",
    "vPKg",

    # Lactate metabolism
    "vLDHn",
    "vLDHg",
    "vLACne",
    "vLACge",
    "vLACgc",

    # Mitochondria
    "vMitoinn",
    "vMitoing",
    "vMitooutn",
    "vMitooutg",

    # NADH shuttle
    "vShuttlen",
    "vShuttleg",

    # Creatine kinase
    "vCKn",
    "vCKg",

    # Oxygen
    "vO2mcn",
    "vO2mcg",

    # Blood / transport
    "Fin",
    "O2c",
    "vO2cc",
    "vGLCcc",
    "vLACcc",

    # Electrophysiology
    "IL",
    "INa",
    "IK",
    "ICa",
    "ImAHP",
    "dIPump",

    # Synaptic stimulation
    "ge",
    "Isyne",
    "Isyni",
    "vnstim",
    "vgstim",
]


print_heading("RATE AUDIT: REST vs STIMULUS ONSET")

print(
    f"{'Rate':<18}"
    f"{'Rest':>22}"
    f"{'Stimulus':>22}"
    f"{'Difference':>22}"
)

print("-" * 84)

for name in important_rates:

    if (
        name in rest_rates
        and name in stim_rates
    ):

        rest_value = float(rest_rates[name])
        stim_value = float(stim_rates[name])

        difference = (
            stim_value - rest_value
        )

        print(
            f"{name:<18}"
            f"{rest_value:>22.12e}"
            f"{stim_value:>22.12e}"
            f"{difference:>22.12e}"
        )

    else:
        print(
            f"{name:<18}"
            f"{'NOT FOUND':>22}"
        )


# ============================================================
# 6. ALL 33 DERIVATIVES
# ============================================================

print_heading(
    "33-STATE DERIVATIVE AUDIT: REST vs STIMULUS"
)

print(
    f"{'State':<18}"
    f"{'Rest dY/dt':>22}"
    f"{'Stim dY/dt':>22}"
    f"{'Difference':>22}"
)

print("-" * 84)

for i in range(33):

    difference = (
        stim_dydt[i]
        - rest_dydt[i]
    )

    print(
        f"Y{i + 1:02d} {state_names[i]:<13}"
        f"{rest_dydt[i]:>22.12e}"
        f"{stim_dydt[i]:>22.12e}"
        f"{difference:>22.12e}"
    )


# ============================================================
# 7. CHECK WHICH STATES RESPOND INSTANTLY
# ============================================================

print_heading(
    "STATES CHANGED IMMEDIATELY BY STIMULATION"
)

instant_change = (
    stim_dydt - rest_dydt
)

threshold = 1e-12

changed_indices = np.where(
    np.abs(instant_change) > threshold
)[0]

if len(changed_indices) == 0:

    print(
        "No derivative changed above threshold."
    )

else:

    for i in changed_indices:

        print(
            f"Y{i + 1:02d} "
            f"{state_names[i]:<14} "
            f"Delta(dY/dt) = "
            f"{instant_change[i]: .12e}"
        )


# ============================================================
# 8. RESTING RESIDUAL RANKING
# ============================================================

print_heading(
    "LARGEST RESTING-STATE RESIDUALS"
)

ranking = np.argsort(
    np.abs(rest_dydt)
)[::-1]

for rank, index in enumerate(
    ranking[:10],
    start=1,
):

    print(
        f"{rank:02d}. "
        f"Y{index + 1:02d} "
        f"{state_names[index]:<14} "
        f"dY/dt = "
        f"{rest_dydt[index]: .12e}"
    )


# ============================================================
# 9. IMPORTANT PARAMETER SNAPSHOT
# ============================================================

important_parameters = [
    "Cm",
    "F",
    "MVF",
    "A",
    "qAK",
    "vATPasesn",
    "vATPasesg",
    "VMaxMitooutn",
    "VMaxMitooutg",
    "kCKnps",
    "kCKgps",
    "TnNADH",
    "TgNADH",
    "TMaxLACne",
]


print_heading(
    "IMPORTANT PARAMETER SNAPSHOT"
)

for name in important_parameters:

    if name in params:
        print(
            f"{name:<18} "
            f"{params[name]: .16e}"
        )
    else:
        print(
            f"{name:<18} NOT FOUND"
        )


# ============================================================
# 10. CONSISTENCY CHECKS
# ============================================================

print_heading(
    "AUTOMATIC CONSISTENCY CHECKS"
)

checks = {}


checks["33 states"] = (
    len(y0) == 33
)

checks["33 resting derivatives"] = (
    len(rest_dydt) == 33
)

checks["33 stimulated derivatives"] = (
    len(stim_dydt) == 33
)

checks["Initial states finite"] = (
    np.all(np.isfinite(y0))
)

checks["Rest derivatives finite"] = (
    np.all(np.isfinite(rest_dydt))
)

checks["Stim derivatives finite"] = (
    np.all(np.isfinite(stim_dydt))
)

checks["Rest ge = 0"] = (
    abs(rest_rates["ge"]) < 1e-12
)

checks["Rest Isyne = 0"] = (
    abs(rest_rates["Isyne"]) < 1e-12
)

checks["Stim ge > 0"] = (
    stim_rates["ge"] > 0
)

checks["Stim Isyne != 0"] = (
    abs(stim_rates["Isyne"]) > 1e-12
)


for name, passed in checks.items():

    status = (
        "PASS"
        if passed
        else "FAIL"
    )

    print(
        f"{status:<6} {name}"
    )


# ============================================================
# 11. KEY REFERENCE VALUES
# ============================================================

print_heading(
    "KEY PYTHON REFERENCE VALUES"
)

print(
    "Resting Vm:",
    y0[27],
)

print(
    "Resting dVm/dt:",
    rest_dydt[27],
)

print(
    "Stimulus-onset dVm/dt:",
    stim_dydt[27],
)

print(
    "Resting neuronal Na derivative:",
    rest_dydt[0],
)

print(
    "Stimulated neuronal Na derivative:",
    stim_dydt[0],
)

print(
    "Resting neuronal ATP derivative:",
    rest_dydt[14],
)

print(
    "Stimulated neuronal ATP derivative:",
    stim_dydt[14],
)

print(
    "Stimulus ge:",
    stim_rates["ge"],
)

print(
    "Stimulus Isyne:",
    stim_rates["Isyne"],
)


# ============================================================
# FINAL RESULT
# ============================================================

print_heading(
    "AUDIT RESULT"
)

if all(checks.values()):

    print(
        "PASS: Python model passed all "
        "automatic translation-audit checks."
    )

else:

    print(
        "CHECK REQUIRED: At least one "
        "automatic audit test failed."
    )

print()
print(
    "NOTE: Passing this test establishes internal "
    "Python consistency."
)

print(
    "It does NOT by itself prove numerical identity "
    "with MATLAB."
)

print(
    "Direct MATLAB reference values are required "
    "for quantitative cross-language validation."
)