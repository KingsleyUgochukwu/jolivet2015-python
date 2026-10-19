
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from parameters import get_parameters

# Load existing Python simulation results
data = np.load("recovery_Ne_240_3600s.npz")
time = data["time"]
states = data["states"]

# Use only the first 100 seconds
mask = time <= 100
time = time[mask]
states = states[mask]

# Extract NADH state variables
nadh_n_cyto = states[:, 12]
nadh_g_cyto = states[:, 13]
nadh_n_mito = states[:, 31]
nadh_g_mito = states[:, 32]

# Read original model parameters
p = get_parameters()

Vn = p["Vn"]
Vg = p["Vg"]
MVF = p["MVF"]

# MATLAB compartment weights
wnM = MVF * Vn
wgM = MVF * Vg
wnC = (1 - MVF) * Vn
wgC = (1 - MVF) * Vg

# Neuronal mitochondrial NADH (%)
nadh_n_mito_pct = (
    100 * (nadh_n_mito - nadh_n_mito[0])
    / nadh_n_mito[0]
)

# Astrocytic cytosolic NADH (%)
nadh_g_cyto_pct = (
    100 * (nadh_g_cyto - nadh_g_cyto[0])
    / nadh_g_cyto[0]
)

# Weighted whole-tissue NADH
nadh_weighted = (
    wnC * nadh_n_cyto
    + wnM * nadh_n_mito
    + wgC * nadh_g_cyto
    + wgM * nadh_g_mito
)

reference = nadh_weighted[0]

# Original MATLAB normalization
nadh_total_pct = (
    100 * (nadh_weighted - reference) / 0.03
)

# Load the original reference curves
reference_path = (
    Path.home()
    / "Downloads"
    / "Jolivet2015_MATLAB"
    / "Jolivet et al. 2015"
    / "Data repository"
    / "results-kasischke.DAT"
)

experimental = np.loadtxt(reference_path)

exp_time = experimental[:, 0]
exp_blue = experimental[:, 1]
exp_red = experimental[:, 2]

# Plot Figure 4A
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    time, nadh_n_mito_pct,
    color="blue", linewidth=1.5,
    label="Neuronal mitochondrial NADH"
)

ax.plot(
    time, nadh_g_cyto_pct,
    color="red", linewidth=1.5,
    label="Astrocytic cytosolic NADH"
)

ax.plot(
    time, nadh_total_pct,
    color="black", linewidth=1.5,
    label="Whole-tissue NADH"
)

# Original reference data
ax.plot(
    exp_time, exp_blue,
    color="blue", linestyle="--",
    linewidth=1.2, label="Reference neuronal NADH"
)

ax.plot(
    exp_time, exp_red,
    color="red", linestyle="--",
    linewidth=1.2, label="Reference astrocytic NADH"
)

# Stimulation interval
ax.axvspan(0, 20, color="lightgray", alpha=0.5)

# Original MATLAB axes
ax.set_xlim(-10, 100)
ax.set_ylim(-12, 12)

ax.set_xlabel("Time (seconds)")
ax.set_ylabel("NADH change (%)")
ax.set_title("Figure 4A — NADH dynamics")

ax.legend(fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "Figure_04A_Python.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# --------------------------------------------------
# FIGURE 4B: OXYGEN AND GLUCOSE UTILIZATION
# --------------------------------------------------

# State variables
o2_n = states[:, 18]
o2_g = states[:, 19]
o2_c = states[:, 20]

glc_n = states[:, 2]
glc_g = states[:, 3]

atp_n = states[:, 14]
atp_g = states[:, 15]

# Oxygen exchange equations
oxygen_term = (
    p["KO2"]
    * (p["HbOP"] / o2_c - 1) ** (-1 / p["nh"])
)

vO2mcn = p["PScapVn"] * (oxygen_term - o2_n)
vO2mcg = p["PScapVg"] * (oxygen_term - o2_g)

# Glucose utilization equations
vHKPFKn = (
    p["kHKPFKn"]
    * atp_n
    * glc_n / (glc_n + p["Kg"])
    / (1 + (atp_n / p["KIATP"]) ** p["nH"])
)

vHKPFKg = (
    p["kHKPFKg"]
    * atp_g
    * glc_g / (glc_g + p["Kg"])
    / (1 + (atp_g / p["KIATP"]) ** p["nH"])
)

# MATLAB Figure 4B scaling
oxygen_n = 1000 * Vn * vO2mcn
oxygen_g = 1000 * Vg * vO2mcg

glucose_n = 6000 * Vn * vHKPFKn
glucose_g = 6000 * Vg * vHKPFKg

# Plot
fig_b, ax_b = plt.subplots(figsize=(8, 5))

ax_b.plot(
    time, oxygen_n,
    color="blue", linewidth=1.5,
    label="Oxygen — neuron"
)

ax_b.plot(
    time, oxygen_g,
    color="red", linewidth=1.5,
    label="Oxygen — astrocyte"
)

ax_b.plot(
    time, glucose_n,
    color="blue", linewidth=0.8,
    label="Glucose × 6 — neuron"
)

ax_b.plot(
    time, glucose_g,
    color="red", linewidth=0.8,
    label="Glucose × 6 — astrocyte"
)

ax_b.axvspan(0, 20, color="lightgray", alpha=0.5)

ax_b.set_xlim(-10, 100)
ax_b.set_ylim(0, 30)

ax_b.set_xlabel("Time (seconds)")
ax_b.set_ylabel("Oxygen / glucose utilization (µM/s)")
ax_b.set_title("Figure 4B — Oxygen and glucose utilization")

ax_b.legend(fontsize=8)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

fig_b.tight_layout()

fig_b.savefig(
    "Figure_04B_Python.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# --------------------------------------------------
# FIGURE 4C: LACTATE TRANSPORT DYNAMICS
# --------------------------------------------------

# Lactate concentrations
lac_n = states[:, 10]
lac_g = states[:, 11]
lac_c = states[:, 22]
lac_e = states[:, 26]

# Pyruvate and cytosolic NADH
pyr_n = states[:, 8]
pyr_g = states[:, 9]

nadh_n = states[:, 12]
nadh_g = states[:, 13]

# Lactate transport rates
vLACne = p["TMaxLACne"] * (
    lac_n / (lac_n + p["KtLACne"])
    - lac_e / (lac_e + p["KtLACne"])
)

vLACge = p["TMaxLACge"] * (
    lac_g / (lac_g + p["KtLACge"])
    - lac_e / (lac_e + p["KtLACge"])
)

vLACgc = p["TMaxLACgc"] * (
    lac_g / (lac_g + p["KtLACgc"])
    - lac_c / (lac_c + p["KtLACgc"])
)

vLACec = p["TMaxLACec"] * (
    lac_e / (lac_e + p["KtLACec"])
    - lac_c / (lac_c + p["KtLACec"])
)

# Lactate dehydrogenase rates
vLDHn = (
    p["kLDHnps"] * pyr_n * nadh_n
    - p["kLDHnms"] * lac_n * (p["N"] - nadh_n)
)

vLDHg = (
    p["kLDHgps"] * pyr_g * nadh_g
    - p["kLDHgms"] * lac_g * (p["N"] - nadh_g)
)

# Apply original MATLAB compartment weighting
lac_ne = 1000 * Vn * vLACne
lac_ge = 1000 * Vg * vLACge
lac_gc = 1000 * Vg * vLACgc
lac_ec = 1000 * p["Ve"] * vLACec

ldh_n = 1000 * Vn * vLDHn
ldh_g = 1000 * Vg * vLDHg

# Plot Figure 4C
fig_c, ax_c = plt.subplots(figsize=(8, 5))

ax_c.plot(
    time, lac_ge,
    color="red", linewidth=1.5,
    label="Astrocyte → extracellular"
)

ax_c.plot(
    time, lac_ne,
    color="blue", linewidth=1.5,
    label="Neuron → extracellular"
)

ax_c.plot(
    time, lac_ec,
    color="magenta", linewidth=1.5,
    label="Extracellular → capillary"
)

# LDH rates are included as thin curves
ax_c.plot(
    time, ldh_n,
    color="blue", linewidth=0.7,
    label="Neuronal LDH"
)

ax_c.plot(
    time, ldh_g,
    color="red", linewidth=0.7,
    label="Astrocytic LDH"
)

# Original MATLAB hides lac_gc
# Therefore, we do not display it.

ax_c.axvspan(0, 20, color="lightgray", alpha=0.5)

ax_c.set_xlim(-10, 100)
ax_c.set_ylim(-8, 8)

ax_c.set_xlabel("Time (seconds)")
ax_c.set_ylabel("Lactate flux (µM/s)")
ax_c.set_title("Figure 4C — Lactate transport dynamics")

ax_c.legend(fontsize=8)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)

fig_c.tight_layout()

fig_c.savefig(
    "Figure_04C_Python.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# --------------------------------------------------
# FIGURE 4D: LACTATE CONCENTRATION CHANGES
# --------------------------------------------------

# Total compartment-weighted lactate
lac_total = (
    Vn * lac_n
    + Vg * lac_g
    + p["Ve"] * lac_e
    + p["Vcap"] * lac_c
)

# Percentage change in extracellular lactate
lac_e_pct = (
    100 * (lac_e - lac_e[0]) / lac_e[0]
)

# Percentage change in total lactate
lac_total_pct = (
    100 * (lac_total - lac_total[0])
    / lac_total[0]
)

# Plot Figure 4D
fig_d, ax_d = plt.subplots(figsize=(8, 5))

ax_d.plot(
    time, lac_e_pct,
    color="magenta",
    linewidth=1.5,
    label="Extracellular lactate"
)

ax_d.plot(
    time, lac_total_pct,
    color="black",
    linewidth=1.5,
    label="Total lactate"
)

ax_d.axvspan(0, 20, color="lightgray", alpha=0.5)

ax_d.axhline(
    0, color="lightgray", linewidth=0.8
)

ax_d.set_xlim(-10, 100)
ax_d.set_ylim(-3, 6)

ax_d.set_xlabel("Time (seconds)")
ax_d.set_ylabel("Lactate change (%)")
ax_d.set_title("Figure 4D — Lactate concentration changes")

ax_d.legend(fontsize=8)

ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)

fig_d.tight_layout()

fig_d.savefig(
    "Figure_04D_Python.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

