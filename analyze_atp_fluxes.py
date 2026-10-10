
import numpy as np

from parameters import get_parameters
from rates import calculate_rates

# Start with the validated reference condition
Ne = 240

data = np.load(f"stimulation_Ne_{Ne}.npz")

time = data["time"]
states = data["states"]

params = get_parameters()

# Analyze the initial state and end of stimulation
for index, label in [(0, "Baseline"), (-1, "20 seconds")]:

    t = float(time[index])
    y = states[index]

    rates = calculate_rates(
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
        Ne=Ne,
        LONG="off",
    )

    glycolytic_production = (
        rates["vPGKn"] + rates["vPKn"]
    )

    mitochondrial_production = (
        3.6 * rates["vMitooutn"]
    )

    glycolytic_investment = (
        2 * rates["vHKPFKn"]
    )

    sodium_pump_consumption = rates["vPumpn"]

    basal_consumption = params["vATPasesn"]

    creatine_kinase = rates["vCKn"]

    net_flux = (
        glycolytic_production
        + mitochondrial_production
        + creatine_kinase
        - glycolytic_investment
        - sodium_pump_consumption
        - basal_consumption
    )

    print(f"\n--- {label} (t={t:.2f}s) ---")
    print(f"Glycolytic ATP production: {glycolytic_production:.6f}")
    print(f"Mitochondrial ATP production: {mitochondrial_production:.6f}")
    print(f"Glycolytic ATP investment: {glycolytic_investment:.6f}")
    print(f"Sodium pump ATP consumption: {sodium_pump_consumption:.6f}")
    print(f"Basal ATP consumption: {basal_consumption:.6f}")
    print(f"Creatine kinase contribution: {creatine_kinase:.6f}")
    print(f"Net ATP flux before buffering: {net_flux:.6f}")
    print("\nATP FLUX COMPARISON ACROSS STIMULATION CONDITIONS")
    print("-" * 75)

    print(
    f"{'Ne':>5} {'Mito ATP':>12} "
    f"{'Glycolysis':>12} {'Na Pump':>12} {'Net Flux':>12}"
)

for Ne_value in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne_value}.npz")

    time = data["time"]
    states = data["states"]

    r = calculate_rates(
        t=float(time[-1]),
        y=states[-1],
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="off",
        F0=0.012,
        species="rat",
        DCBF=0.4,
        Ne=Ne_value,
        LONG="off",
    )

    mito = 3.6 * r["vMitooutn"]
    glycolysis = r["vPGKn"] + r["vPKn"]
    pump = r["vPumpn"]

    net = (
        glycolysis
        + mito
        + r["vCKn"]
        - 2 * r["vHKPFKn"]
        - pump
        - params["vATPasesn"]
    )

    print(
        f"{Ne_value:>5} {mito:>12.6f} "
        f"{glycolysis:>12.6f} {pump:>12.6f} {net:>12.6f}"
    )
    
# Check ATP flux at the minimum ATP concentration
Ne_value = 240

data = np.load(f"stimulation_Ne_{Ne_value}.npz")

time = data["time"]
states = data["states"]

min_index = np.argmin(states[:, 14])

t = float(time[min_index])
y = states[min_index]

r = calculate_rates(
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
    Ne=Ne_value,
    LONG="off",
)

net_flux = (
    r["vPGKn"]
    + r["vPKn"]
    + 3.6 * r["vMitooutn"]
    + r["vCKn"]
    - 2 * r["vHKPFKn"]
    - r["vPumpn"]
    - params["vATPasesn"]
)

print("\nATP MINIMUM CONSISTENCY CHECK")
print(f"Ne = {Ne_value}")
print(f"ATP minimum time = {t:.4f} s")
print(f"Minimum ATP = {y[14]:.6f}")
print(f"Net ATP flux at minimum = {net_flux:.10f}")
print("\nATP MINIMUM VALIDATION ACROSS ALL CONDITIONS")
print("-" * 65)

for Ne_value in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne_value}.npz")

    time = data["time"]
    states = data["states"]

    idx = np.argmin(states[:, 14])

    t = float(time[idx])
    y = states[idx]

    r = calculate_rates(
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
        Ne=Ne_value,
        LONG="off",
    )

    net_flux = (
        r["vPGKn"]
        + r["vPKn"]
        + 3.6 * r["vMitooutn"]
        + r["vCKn"]
        - 2 * r["vHKPFKn"]
        - r["vPumpn"]
        - params["vATPasesn"]
    )

    print(
        f"Ne={Ne_value}: "
        f"ATP minimum={y[14]:.6f}, "
        f"Time={t:.4f}s, "
        f"Net flux={net_flux:.10f}"
    )
    print("\nMETABOLIC FLUXES AT ATP MINIMUM")
print("-" * 85)

print(
    f"{'Ne':>5} {'Mito ATP':>12} "
    f"{'Net Glycolysis':>16} "
    f"{'Na Pump':>12} {'Creatine Kinase':>17}"
)

for Ne_value in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne_value}.npz")

    time = data["time"]
    states = data["states"]

    idx = np.argmin(states[:, 14])

    r = calculate_rates(
        t=float(time[idx]),
        y=states[idx],
        params=params,
        TIME=0.0,
        t1=0.0,
        tstim=20.0,
        CBF="off",
        F0=0.012,
        species="rat",
        DCBF=0.4,
        Ne=Ne_value,
        LONG="off",
    )

    mito = 3.6 * r["vMitooutn"]

    net_glycolysis = (
        r["vPGKn"]
        + r["vPKn"]
        - 2 * r["vHKPFKn"]
    )

    pump = r["vPumpn"]
    ck = r["vCKn"]

    print(
        f"{Ne_value:>5} {mito:>12.6f} "
        f"{net_glycolysis:>16.6f} "
        f"{pump:>12.6f} {ck:>17.6f}"
    )
    
# Creatine kinase consistency check at PCr minimum

data = np.load("recovery_Ne_240.npz")

time = data["time"]
states = data["states"]

idx = np.argmin(states[:, 16])

t = float(time[idx])
y = states[idx]

r = calculate_rates(
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

print("\nPHOSPHOCREATINE MINIMUM CONSISTENCY CHECK")
print(f"PCr minimum time = {t:.6f} s")
print(f"Minimum PCr = {y[16]:.9f}")
print(f"ATP at PCr minimum = {y[14]:.9f}")
print(f"Creatine kinase flux = {r['vCKn']:.12f}")
print("\nCREATINE KINASE FLUX DURING RECOVERY")
print("-" * 55)

data = np.load("recovery_Ne_240.npz")
time = data["time"]
states = data["states"]

for target_time in [20, 30, 40, 50, 60, 80, 120]:

    idx = np.argmin(np.abs(time - target_time))

    r = calculate_rates(
        t=float(time[idx]),
        y=states[idx],
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

    ck_flux = r["vCKn"]

    print(
        f"Time={time[idx]:.2f}s | "
        f"CK flux={ck_flux:+.8f} | "
        f"PCr={states[idx,16]:.6f}"
    )
