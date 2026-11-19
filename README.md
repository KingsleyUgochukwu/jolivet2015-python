# Python Reimplementation and Computational Analysis of the Jolivet et al. (2015) Neuron–Glia–Vasculature Model

## Project overview

This repository documents a Python implementation and computational investigation of the neuron–glia–vasculature metabolic model described by Jolivet et al. (2015). The work has two connected parts:

1. Reproduction and analysis of activity-dependent brain energy metabolism, including neuronal stimulation, ATP production and consumption, and phosphocreatine (PCr) dynamics.
2. An extension of the model to investigate transient cerebral ischemia and reperfusion, with emphasis on cellular ATP responses, extracellular glucose and lactate, and lactate transport and mass balance.

The original framework couples neuronal activity, astrocytic metabolism, substrate exchange, and vascular responses. The Python work aims to reproduce the published model behaviour and use the implementation to conduct additional, clearly identified computational experiments.

This is an independent computational research and training project developed by Kingsley Nnaemeka Ugochukwu, MSc student in Biochemistry and Molecular Biology at Obafemi Awolowo University, Nigeria.

## Scientific background and questions

Neuronal activity increases ATP demand, particularly through ion-pumping processes that restore ionic gradients. The neuron–glia–vasculature model provides a mechanistic framework for examining how neuronal and astrocytic metabolism and vascular substrate delivery interact to meet this demand.

The project addresses two groups of questions:

**Activity-dependent metabolism**
- How does increasing neuronal stimulation affect sodium accumulation, Na⁺/K⁺-ATPase activity, and neuronal ATP?
- How do mitochondrial ATP production and creatine kinase-mediated PCr buffering respond to increased energy demand?
- How similar are normalized PCr recovery trajectories across the tested stimulation conditions?

**Ischemia and reperfusion**
- How does the duration of reduced cerebral blood flow affect neuronal and astrocytic ATP depletion and recovery?
- How do extracellular glucose, lactate, and capillary oxygen vary during ischemia and reperfusion?
- How do neuronal, astrocytic, and vascular lactate transport contribute to the extracellular lactate balance?

## Original model and Python implementation

The implementation is based on the mathematical framework in:

Jolivet, R., Coggan, J. S., Allaman, I., & Magistretti, P. J. (2015). Multi-timescale modeling of activity-dependent metabolic coupling in the neuron-glia-vasculature ensemble. *PLOS Computational Biology, 11*(2), e1004036. https://doi.org/10.1371/journal.pcbi.1004036

The model implementation contains 33 state variables and uses the parameterization and coupled balance equations of the original framework. Core Python modules include:

- `parameters.py` — model parameters.
- `initial_state.py` — initial conditions.
- `rates.py` — metabolic reaction rates and flux expressions.
- `balance_equations.py` — differential equations for the model state variables.
- `conductance.py` — neuronal conductance-related calculations.
- `blood_flow.py` — vascular and blood-flow-related calculations.
- `simulation_runner.py` — numerical simulation routines.
- `simulate_full.py` — full-model simulation.
- `compare_matlab_python.py` — internal consistency audit of the Python implementation (not a quantitative MATLAB comparison).

The stiff ordinary differential equation system is integrated in Python using SciPy's BDF solver.

## Python implementation checks and validation status

The current `compare_matlab_python.py` script performs an internal audit of the Python implementation. It checks the 33-element state vector, initial conditions, resting and stimulated rate/derivative calculations, finite numerical values, and selected expected stimulation-onset behaviours.

These checks are useful for detecting implementation inconsistencies within the Python model. They do **not** constitute a quantitative comparison with the original MATLAB implementation. Direct MATLAB reference trajectories or numerical outputs are still required to establish cross-language agreement. Systematic reproduction of the published model figures and comparison with the corresponding reference results therefore remain important validation tasks.

The original MATLAB source code is not redistributed in this repository.

## Part I — Activity-dependent brain energy metabolism

### Stimulation protocol and metabolic responses

Five neuronal stimulation conditions were investigated: `Ne = 120, 180, 240, 300, 360`.

Across the tested range, stronger stimulation was associated with:
- Greater neuronal sodium accumulation.
- Increased Na⁺/K⁺-ATPase activity and pump demand.
- Greater neuronal ATP depletion during stimulation.
- Increased mitochondrial ATP production.
- Greater neuronal PCr depletion.

These are responses within the model and its selected parameterization. They illustrate the relationship between neuronal activation, ion homeostasis, ATP demand, and metabolic compensation represented by the equations.

### ATP production and consumption

The ATP-flux analysis separates contributions from glycolytic ATP production, mitochondrial ATP production, Na⁺/K⁺-ATPase consumption, basal ATP consumption, and creatine kinase-mediated buffering. This provides a flux-based interpretation of ATP concentration changes rather than relying on ATP concentration alone.

Relevant analysis scripts include:

- `analyze_experiments.py`
- `analyze_atp_fluxes.py`
- `mitochondrial_time_course.py`

### Phosphocreatine depletion and recovery

PCr is included as part of the model's cellular energy-buffering system through the creatine kinase reaction. Simulations were extended to examine PCr depletion and recovery following neuronal stimulation.

The following minima were obtained from the extended PCr time-course analyses:

| Stimulation (Ne) | Minimum PCr | Time of minimum (s) |
|---|---:|---:|
| 120 | 4.875640 | 55.12 |
| 180 | 4.786364 | 50.74 |
| 240 | 4.645318 | 47.40 |
| 300 | 4.443306 | 44.45 |
| 360 | 4.173524 | 44.35 |

Stronger stimulation produced greater PCr depletion in these simulations. To compare recovery independently of the absolute depletion magnitude, normalized recovery was calculated as:

`Recovery (%) = 100 × (PCr(t) − PCr_min) / (PCr_initial − PCr_min)`

The normalized trajectories were compared with the `Ne = 240` reference condition using root mean square deviation (RMSD):

| Ne | RMSD (percentage points) |
|---|---:|
| 120 | 0.5705 |
| 180 | 0.1669 |
| 240 | 0.0000 |
| 300 | 0.0607 |
| 360 | 0.1645 |

The small RMSD values indicate broadly similar fractional PCr recovery trajectories across the tested stimulation conditions. This is a model-specific computational observation, not a claim that biological PCr recovery is universally independent of neuronal activity.

Relevant scripts include:

- `analyze_phosphocreatine.py`
- `compare_normalized_pcr_recovery.py`
- `calculate_pcr_rmsd.py`

### Interpretation

Taken together, the stimulation simulations suggest that stronger neuronal activation increases ATP demand and PCr utilization, while mitochondrial ATP production responds to the increased energy requirement. In the model, neuronal ATP concentration can approach its reference level before the PCr pool has fully recovered. ATP availability and recovery of energy-buffering reserves therefore provide related but distinct views of the simulated metabolic response.

## Part II — Transient cerebral ischemia and reperfusion

### Research objective

The model was extended to investigate how the duration of reduced cerebral blood flow influences cellular ATP, extracellular glucose and lactate dynamics, and neuronal–astrocytic lactate transport during ischemia and subsequent reperfusion.

### Simulation protocol

The ischemia–reperfusion simulations used the following settings:

| Parameter | Setting |
|---|---:|
| Baseline cerebral blood flow | 0.012 |
| Ischemia onset | 100 s |
| Cerebral blood flow reduction | 80% |
| Cerebral blood flow during ischemia | 0.0024 |
| Ischemia durations | 300, 600, 1000, 1500, 2000 s |
| Reperfusion | Restoration of baseline flow |
| Reperfusion observation period | 1000 s |

The 1000-second ischemia condition was included as an anchor for comparison with an earlier GAMMA-MRI modelling report. The other durations are exploratory simulation conditions; they are not clinical thresholds.

### ATP depletion and recovery

At the end of ischemia, neuronal ATP depletion ranged from approximately 96.68% to 99.32% across the five simulated durations. Astrocytic ATP depletion ranged from approximately 79.07% to 93.81%.

After the fixed 1000-second reperfusion period, neuronal ATP was approximately 99.96–99.98% of its pre-ischemia baseline. Astrocytic ATP was approximately 102.01–102.04% of its pre-ischemia baseline, representing a modest model-predicted overshoot. The separate recovery metric, calculated as the fraction of the initial ATP depletion replenished, was approximately 99.96–99.98% for neuronal ATP and 102.14–102.58% for astrocytic ATP. Values above 100% indicate that the simulated ATP concentration exceeded its pre-ischemia baseline. This overshoot should not be interpreted as evidence of improved cellular health or tissue protection.

### Other metabolic variables

The simulations also examined extracellular glucose and lactate, capillary oxygen, lactate dehydrogenase reaction rates, and neuronal, astrocytic, and vascular lactate transport.

The time courses show that the metabolic response depends on both the ischemia duration and the time elapsed during reperfusion. Comparisons between different ischemia durations must therefore specify whether values are taken at the end of ischemia, at a common reperfusion time, or at a trajectory-specific peak. End-of-ischemia values from different-duration simulations are not measurements at a common elapsed time. The script `analyse_fixed_reperfusion_lactate.py` interpolates each lactate-flux time course at 0, 60, 300, 600 and 900 s after reperfusion onset and writes `lactate_fixed_reperfusion_comparison.csv`. This provides matched-time comparisons across the five ischemia durations.

### Extracellular lactate transport and mass balance

The extracellular lactate balance was decomposed into three transport contributions:

1. Neuronal-to-extracellular lactate transport.
2. Astrocytic-to-extracellular lactate transport.
3. Extracellular-to-capillary lactate transport.

Using the model's transport sign convention, the extracellular lactate balance is:

[
rac{d[LAC]_e}{dt}
=
rac{v_{LACne}}{r_{en}}
+
rac{v_{LACge}}{r_{eg}}
-
v_{LACec}
]

The analysis of the five ischemia durations showed that:

- Astrocytic transport made a positive cumulative contribution to extracellular lactate over the reperfusion period in all five simulations.
- Neuronal transport showed transient positive contributions during reperfusion in the longer ischemia simulations, although its cumulative contribution over the complete reperfusion period was negative.
- Extracellular-to-capillary transport contributed negatively to the extracellular lactate balance.
- The integrated transport contributions reproduced the simulated change in extracellular lactate concentration with very small numerical differences.

Astrocytic transport accounted for approximately 48.5–49.5% of the combined absolute integrated transport magnitude across the five conditions. This percentage describes the relative magnitude of the integrated transport contributions; it is **not** the proportion of total lactate production attributable to astrocytes.

The integrated net transport balance closely matched the actual change in extracellular lactate over the 1000-second reperfusion interval. This agreement is a numerical mass-balance consistency check within the model, rather than independent experimental validation.

### Comparison with experimental studies

The model's glucose and lactate responses were compared descriptively with published transient cerebral ischemia studies, including:

- Frykholm et al. (2005), a study of cerebral blood flow, oxygen metabolism, and extracellular glucose and lactate during middle cerebral artery occlusion and reperfusion in nonhuman primates.
- Havelund et al. (2019), a study of brain metabolites during mitochondrial inhibition and transient cerebral ischemia in rats.

This comparison was exploratory and was not intended as formal model validation. Differences in species, ischemia protocols, brain regions, sampling methods, and observation periods limit direct quantitative comparison.

### Figures and results

The analysis generated figures covering:
- ATP, glucose, lactate, and oxygen time courses.
- Exploratory comparison with published experimental glucose and lactate observations.
- Cellular lactate dehydrogenase fluxes.
- Lactate transport fluxes.
- Extracellular lactate balance and mass-balance consistency.

### Reproducibility

The main ischemia–reperfusion simulation and analysis entry points are:

```bash
# Run an ischemia–reperfusion simulation
python simulate_ischemia.py

# Run the ischemia-duration sweep (one duration per run)
python simulate_ischemia_sweep.py 300
python simulate_ischemia_sweep.py 600
python simulate_ischemia_sweep.py 1000
python simulate_ischemia_sweep.py 1500
python simulate_ischemia_sweep.py 2000

# Analyse the duration sweep and lactate fluxes
python analyse_ischemia_sweep.py
python analyse_lactate_fluxes.py

# Compare lactate and transport contributions at matched reperfusion times
python analyse_fixed_reperfusion_lactate.py

# Analyse extracellular lactate mass balance and contributions
python analyse_lactate_mass_balance.py
python analyse_reperfusion_contribution_peaks.py
python analyse_cumulative_lactate_contributions.py
python analyse_relative_lactate_contributions.py
```

The analysis scripts depend on the corresponding simulation output files. Generate the required simulation datasets before running downstream analyses. The exact execution order and input filenames should be checked in each script.

## Limitations and future directions

The results presented here are predictions from a mechanistic computational model. They do not independently establish neuronal survival, tissue protection, or clinical outcomes.

Key limitations include:
- The current comparison script performs internal consistency checks on the Python implementation; quantitative comparison against MATLAB reference outputs has not yet been completed.
- Additional stimulation and ischemia–reperfusion simulations are model-based predictions and have not been independently validated against new experimental data.
- Results depend on the model structure, parameterization, initial conditions, and numerical integration settings.
- The current model does not explicitly represent glycogen metabolism. Adding glycogen metabolism or the pentose phosphate pathway would require further model development and evaluation.

Potential future work includes parameter sensitivity analysis, additional quantitative comparisons with experimental measurements, and carefully evaluated extensions to the metabolic pathways represented in the model.

## References

Frykholm, P., Hillered, L., Långström, B., Persson, L., Valtysson, J., & Enblad, P. (2005). Relationship between cerebral blood flow and oxygen metabolism, and extracellular glucose and lactate concentrations during middle cerebral artery occlusion and reperfusion: A microdialysis and positron emission tomography study in nonhuman primates. *Journal of Neurosurgery, 102*(6), 1076–1084. https://doi.org/10.3171/jns.2005.102.6.1076

Havelund, J. F., Nygaard, K. H., Nielsen, T. H., Nordström, C.-H., Poulsen, F. R., Færgeman, N. J., Forsse, A., & Gramsbergen, J. B. (2019). In vivo microdialysis of endogenous and ¹³C-labeled TCA metabolites in rat brain: Reversible and persistent effects of mitochondrial inhibition and transient cerebral ischemia. *Metabolites, 9*(10), 204. https://doi.org/10.3390/metabo9100204

Jolivet, R. (2023). *D5.1: Report on biophysical model* (Version 2.0). GAMMA-MRI, Horizon 2020 Research and Innovation Programme, Grant Agreement No. 964644. Submitted 4 December 2023.

Jolivet, R., Coggan, J. S., Allaman, I., & Magistretti, P. J. (2015). Multi-timescale modeling of activity-dependent metabolic coupling in the neuron-glia-vasculature ensemble. *PLOS Computational Biology, 11*(2), e1004036. https://doi.org/10.1371/journal.pcbi.1004036

## Attribution

The original mathematical model and scientific framework were developed by Jolivet and colleagues. This repository documents my independent Python implementation, internal consistency checks, and subsequent computational analyses. The original MATLAB source code is not redistributed in this repository.

## Author

**Kingsley Nnaemeka Ugochukwu**  
MSc Student, Biochemistry and Molecular Biology  
Obafemi Awolowo University, Ile-Ife, Nigeria

Research interests: Computational Biology, Systems Neuroscience, Brain Energy Metabolism, and Mathematical Modeling.

GitHub: https://github.com/KingsleyUgochukwu
