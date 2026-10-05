import numpy as np


# ============================================================
# LOAD RESULTS
# ============================================================

data = np.load("jolivet_full_results.npz")

time = data["time"]
states = data["states"]


print()
print("========================================")
print("JOLIVET 2015 PYTHON MODEL ANALYSIS")
print("========================================")

print()
print("--- SAVED EXPERIMENT SETTINGS ---")

if "CBF" in data:
    print("CBF:", str(data["CBF"]))
else:
    print("CBF: not stored")

if "tstim" in data:
    tstim = float(data["tstim"])
    print("Stimulation duration:", tstim, "s")
else:
    tstim = 20.0
    print("Stimulation duration: not stored")

if "Ne" in data:
    print("Ne:", int(data["Ne"]))
else:
    print("Ne: not stored")


print()
print("--- DATASET ---")
print("Time points:", len(time))
print("State array shape:", states.shape)
print("Start time:", time[0], "s")
print("End time:", time[-1], "s")


# ============================================================
# EXTRACT VARIABLES
# ============================================================

na = states[:, 0]
atp = states[:, 14]
vm = states[:, 27]
ca = states[:, 30]
nadh = states[:, 31]


# ============================================================
# FIND EXTREMA
# ============================================================

na_max_idx = np.argmax(na)
atp_min_idx = np.argmin(atp)
ca_max_idx = np.argmax(ca)
nadh_min_idx = np.argmin(nadh)
vm_max_idx = np.argmax(vm)


# ============================================================
# INITIAL VALUES
# ============================================================

na0 = na[0]
atp0 = atp[0]
vm0 = vm[0]
ca0 = ca[0]
nadh0 = nadh[0]


# ============================================================
# FINAL VALUES
# ============================================================

na_final = na[-1]
atp_final = atp[-1]
vm_final = vm[-1]
ca_final = ca[-1]
nadh_final = nadh[-1]


# ============================================================
# PERCENT CHANGES
# ============================================================

na_percent = (
    (na[na_max_idx] - na0)
    / na0
    * 100
)

atp_percent = (
    (atp0 - atp[atp_min_idx])
    / atp0
    * 100
)

ca_percent = (
    (ca[ca_max_idx] - ca0)
    / ca0
    * 100
)

nadh_percent = (
    (nadh0 - nadh[nadh_min_idx])
    / nadh0
    * 100
)


# ============================================================
# END-OF-STIMULATION VALUES
# ============================================================

stim_idx = np.argmin(
    np.abs(time - tstim)
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("========================================")
print("NEURONAL ATP")
print("========================================")

print("Initial ATP:", atp0)
print("Minimum ATP:", atp[atp_min_idx])
print(
    "Time of minimum ATP:",
    time[atp_min_idx],
    "s",
)
print(
    "ATP decrease:",
    atp_percent,
    "%"
)
print(
    "ATP at end of stimulation:",
    atp[stim_idx],
)
print("Final ATP:", atp_final)


print()
print("========================================")
print("NEURONAL SODIUM")
print("========================================")

print("Initial Na+:", na0)
print("Maximum Na+:", na[na_max_idx])
print(
    "Time of maximum Na+:",
    time[na_max_idx],
    "s",
)
print(
    "Na+ increase:",
    na_percent,
    "%"
)
print(
    "Na+ at end of stimulation:",
    na[stim_idx],
)
print("Final Na+:", na_final)


print()
print("========================================")
print("INTRACELLULAR CALCIUM")
print("========================================")

print("Initial Ca2+:", ca0)
print("Maximum Ca2+:", ca[ca_max_idx])
print(
    "Time of maximum Ca2+:",
    time[ca_max_idx],
    "s",
)
print(
    "Ca2+ increase:",
    ca_percent,
    "%"
)
print(
    "Ca2+ at end of stimulation:",
    ca[stim_idx],
)
print("Final Ca2+:", ca_final)


print()
print("========================================")
print("MITOCHONDRIAL NADH")
print("========================================")

print("Initial NADH:", nadh0)
print("Minimum NADH:", nadh[nadh_min_idx])
print(
    "Time of minimum NADH:",
    time[nadh_min_idx],
    "s",
)
print(
    "NADH decrease:",
    nadh_percent,
    "%"
)
print(
    "NADH at end of stimulation:",
    nadh[stim_idx],
)
print("Final NADH:", nadh_final)


print()
print("========================================")
print("MEMBRANE POTENTIAL")
print("========================================")

print("Initial Vm:", vm0)
print("Maximum Vm:", vm[vm_max_idx])
print(
    "Time of maximum Vm:",
    time[vm_max_idx],
    "s",
)
print(
    "Vm at end of stimulation:",
    vm[stim_idx],
)
print("Final Vm:", vm_final)


# ============================================================
# RECOVERY ERROR
# ============================================================

print()
print("========================================")
print("RECOVERY TO BASELINE")
print("========================================")

print(
    "Na+ final difference:",
    na_final - na0,
)

print(
    "ATP final difference:",
    atp_final - atp0,
)

print(
    "Vm final difference:",
    vm_final - vm0,
)

print(
    "Ca2+ final difference:",
    ca_final - ca0,
)

print(
    "NADH final difference:",
    nadh_final - nadh0,
)