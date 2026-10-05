import numpy as np


# ============================================================
# MATLAB vs PYTHON VALIDATION
# Jolivet et al. model
# ============================================================

# ------------------------------------------------------------
# Load Python simulation
# ------------------------------------------------------------

data = np.load("jolivet_full_results.npz")

time = data["time"]
states = data["states"]

print()
print("=" * 68)
print("JOLIVET MODEL: MATLAB vs PYTHON VALIDATION")
print("=" * 68)


# ------------------------------------------------------------
# State indices
# ------------------------------------------------------------

NA = 0
ATP = 14
VM = 27
CA = 30
NADH_MITO = 31


# ------------------------------------------------------------
# Extract Python variables
# ------------------------------------------------------------

na = states[:, NA]
atp = states[:, ATP]
vm = states[:, VM]
ca = states[:, CA]
nadh = states[:, NADH_MITO]


# ------------------------------------------------------------
# Find stimulation endpoint
# ------------------------------------------------------------

tstim = 20.0

stim_end_index = np.argmin(
    np.abs(time - tstim)
)


# ------------------------------------------------------------
# Calculate Python response statistics
# ------------------------------------------------------------

python_results = {
    "ATP minimum": np.min(atp),
    "ATP minimum time": time[np.argmin(atp)],

    "Na maximum": np.max(na),
    "Na maximum time": time[np.argmax(na)],

    "Ca maximum": np.max(ca),
    "Ca maximum time": time[np.argmax(ca)],

    "NADH minimum": np.min(nadh),
    "NADH minimum time": time[np.argmin(nadh)],

    "Vm maximum": np.max(vm),
    "Vm maximum time": time[np.argmax(vm)],

    "Na at 20 s": na[stim_end_index],
    "ATP at 20 s": atp[stim_end_index],
    "Vm at 20 s": vm[stim_end_index],
    "Ca at 20 s": ca[stim_end_index],
    "NADH at 20 s": nadh[stim_end_index],

    "Final Na": na[-1],
    "Final ATP": atp[-1],
    "Final Vm": vm[-1],
    "Final Ca": ca[-1],
    "Final NADH": nadh[-1],
}


# ------------------------------------------------------------
# Reference MATLAB results
# Obtained from original MATLAB implementation
# ------------------------------------------------------------

matlab_results = {
    "ATP minimum": 2.09031,
    "ATP minimum time": 8.3577,

    "Na maximum": 11.9409,
    "Na maximum time": 5.5562,

    "Ca maximum": 0.000938277,
    "Ca maximum time": 0.1118,

    "NADH minimum": 0.111791,
    "NADH minimum time": 6.8528,

    "Vm maximum": 70.5308,
    "Vm maximum time": 0.0135,

    "Na at 20 s": 9.86454,
    "ATP at 20 s": 2.16607,
    "Vm at 20 s": -64.8771,
    "Ca at 20 s": 5.66369e-05,
    "NADH at 20 s": 0.115788,

    "Final Na": 7.97372,
    "Final ATP": 2.19736,
    "Final Vm": -73.5864,
    "Final Ca": 5.10083e-05,
    "Final NADH": 0.123447,
}


# ------------------------------------------------------------
# Percentage difference function
# ------------------------------------------------------------

def percent_difference(reference, test):

    if reference == 0:
        return np.nan

    return abs(
        (test - reference) / reference
    ) * 100


# ------------------------------------------------------------
# Print comparison table
# ------------------------------------------------------------

print()
print(
    f"{'MEASUREMENT':<24}"
    f"{'MATLAB':>14}"
    f"{'PYTHON':>14}"
    f"{'ERROR (%)':>14}"
)

print("-" * 68)


errors = {}


for key in matlab_results:

    matlab_value = matlab_results[key]
    python_value = python_results[key]

    error = percent_difference(
        matlab_value,
        python_value,
    )

    errors[key] = error

    print(
        f"{key:<24}"
        f"{matlab_value:>14.6g}"
        f"{python_value:>14.6g}"
        f"{error:>14.4f}"
    )


# ------------------------------------------------------------
# Validation summary
# ------------------------------------------------------------

print()
print("=" * 68)
print("VALIDATION SUMMARY")
print("=" * 68)


# Main metabolic / ionic quantities used for validation.
# Fast electrophysiological peak timing is reported separately
# because MATLAB ode15s and SciPy BDF need not locate the exact
# spike maximum at identical internal solver times.

primary_metrics = [
    "ATP minimum",
    "Na maximum",
    "Ca maximum",
    "NADH minimum",
    "Na at 20 s",
    "ATP at 20 s",
    "Vm at 20 s",
    "Ca at 20 s",
    "NADH at 20 s",
    "Final Na",
    "Final ATP",
    "Final Vm",
    "Final Ca",
    "Final NADH",
]


primary_errors = [
    errors[key]
    for key in primary_metrics
]


mean_error = np.mean(primary_errors)
max_error = np.max(primary_errors)


print(
    f"Mean error across primary metrics: "
    f"{mean_error:.4f} %"
)

print(
    f"Maximum error across primary metrics: "
    f"{max_error:.4f} %"
)


# ------------------------------------------------------------
# Simple validation criterion
# ------------------------------------------------------------

tolerance = 1.0

print()
print(
    f"Validation tolerance: {tolerance:.1f}% "
    f"for primary response magnitudes"
)


if max_error <= tolerance:

    print()
    print("VALIDATION PASSED")
    print(
        "Python reproduces the MATLAB reference "
        "within the selected tolerance."
    )

else:

    print()
    print("VALIDATION REQUIRES REVIEW")
    print(
        "At least one primary response magnitude "
        "exceeds the selected tolerance."
    )


# ------------------------------------------------------------
# Timing diagnostics
# ------------------------------------------------------------

print()
print("=" * 68)
print("PEAK TIMING DIAGNOSTICS")
print("=" * 68)

timing_metrics = [
    "ATP minimum time",
    "Na maximum time",
    "Ca maximum time",
    "NADH minimum time",
    "Vm maximum time",
]

for key in timing_metrics:

    matlab_value = matlab_results[key]
    python_value = python_results[key]

    absolute_difference = abs(
        python_value - matlab_value
    )

    print(
        f"{key:<24}"
        f"MATLAB = {matlab_value:>10.6f} s   "
        f"Python = {python_value:>10.6f} s   "
        f"|Δt| = {absolute_difference:.6f} s"
    )


print()
print("=" * 68)
print("VALIDATION COMPLETE")
print("=" * 68)