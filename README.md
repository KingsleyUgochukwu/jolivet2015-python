# Python Reimplementation and Computational Analysis of the Jolivet et al. (2015) Neuron–Glia–Vasculature Model

## Project Overview

This project presents a Python reimplementation and computational investigation of the neuron–glia–vasculature metabolic model described by Jolivet et al. (2015).

The original model integrates neuronal electrical activity, astrocytic metabolism, energy substrate exchange, and vascular responses to investigate activity-dependent metabolic coupling in the brain.

My objective is to reproduce the published mathematical framework in Python, evaluate its numerical behavior against the original MATLAB implementation, and investigate how changes in neuronal stimulation influence ATP availability, ATP-consuming processes, mitochondrial ATP production, and phosphocreatine (PCr) recovery.

This project forms part of my independent training in computational biology and systems neuroscience, building on my academic background in Biochemistry and Molecular Biology.

## Scientific Background

Neuronal activity creates an increased demand for ATP, particularly through ion transport mechanisms required to restore ionic gradients following electrical activity.

The Jolivet et al. (2015) framework provides a mechanistic representation of the interactions between neuronal metabolism, astrocytic metabolism, and vascular substrate delivery.

The present computational investigation focuses on three connected questions:

1. How does increasing neuronal stimulation affect ATP availability and Na⁺/K⁺-ATPase activity?
2. How do mitochondrial ATP production and creatine kinase-mediated phosphocreatine buffering respond to increased energy demand?
3. Does the fractional recovery of neuronal phosphocreatine depend strongly on stimulation intensity?

## Python Implementation

The project implements the model's coupled differential equations and associated metabolic rate expressions in Python.

Key components include:

- `parameters.py` — model parameters.
- `initial_state.py` — initial conditions for the state variables.
- `rates.py` — metabolic reaction rates and flux expressions.
- `balance_equations.py` — differential equations governing model dynamics.
- `conductance.py` — neuronal conductance-related calculations.
- `blood_flow.py` — vascular and blood-flow-related calculations.
- `simulation_runner.py` — numerical simulation routines.
- `simulate_full.py` — full-model simulation.
- `compare_matlab_python.py` — comparison with MATLAB reference outputs.

The model is numerically integrated using SciPy's stiff ODE solver (`BDF`).

## MATLAB–Python Validation

A major priority was to establish whether the Python implementation reproduces the original model's behavior before conducting additional analyses.

For a selected set of 14 response and endpoint metrics under the reference simulation conditions, the comparison yielded:

- Mean percentage error: 0.0185%
- Maximum percentage error: 0.1507%

These results indicate close agreement for the selected metrics. They should not be interpreted as proof of complete numerical equivalence across all state variables, transient peaks, or simulation settings.

## Stimulation-Dependent Metabolic Responses

Five neuronal stimulation conditions were investigated:

`Ne = 120, 180, 240, 300, 360`

Across this range, stronger stimulation was associated with:

- Greater neuronal sodium accumulation.
- Increased Na⁺/K⁺-ATPase demand.
- Greater neuronal ATP depletion during stimulation.
- Increased mitochondrial ATP production.
- Greater neuronal phosphocreatine depletion.

The findings are consistent with the model's representation of activity-dependent energetic demand and metabolic compensation.

### ATP Production and Consumption

The analysis separates contributions from glycolytic ATP production, mitochondrial ATP production, Na⁺/K⁺-ATPase consumption, basal ATP consumption, and creatine kinase-mediated buffering.

The purpose is to interpret ATP concentration changes in terms of the metabolic fluxes represented in the mathematical model, rather than ATP concentration alone.

Relevant scripts:

- `analyze_experiments.py`
- `analyze_atp_fluxes.py`
- `mitochondrial_time_course.py`

## Phosphocreatine Depletion and Recovery

Phosphocreatine provides an important component of the modeled cellular energy-buffering system through the creatine kinase reaction.

Simulations were extended to examine the recovery of neuronal phosphocreatine following stimulation.

For the five stimulation conditions, the neuronal phosphocreatine minima were approximately:

| Stimulation (Ne) | Minimum PCr | Time of minimum (s) |
|---|---:|---:|
| 120 | 4.875640 | 55.12 |
| 180 | 4.786364 | 50.74 |
| 240 | 4.645318 | 47.40 |
| 300 | 4.443306 | 44.45 |
| 360 | 4.173524 | 44.35 |

Although stronger stimulation produced greater phosphocreatine depletion, normalized recovery trajectories were similar across the tested conditions.

To compare recovery independently of the absolute depletion magnitude, the following normalization was used:

`Recovery (%) = 100 × (PCr(t) − PCr_min) / (PCr_initial − PCr_min)`

The normalized recovery trajectories were compared against the reference condition (`Ne = 240`) using root mean square deviation (RMSD).

| Ne | RMSD (percentage points) |
|---|---:|
| 120 | 0.5705 |
| 180 | 0.1669 |
| 240 | 0.0000 |
| 300 | 0.0607 |
| 360 | 0.1645 |

These small differences suggest that the model exhibits broadly similar fractional PCr recovery trajectories across the investigated stimulation range.

This is a computational observation within the tested model configuration, not a claim that biological recovery kinetics are universally independent of neuronal activity.

Relevant scripts:

- `analyze_phosphocreatine.py`
- `compare_pcr_recovery.py`
- `compare_normalized_pcr_recovery.py`
- `calculate_pcr_rmsd.py`

## Research Interpretation

The simulations suggest that stronger neuronal activation increases ATP demand and phosphocreatine utilization, while mitochondrial ATP production responds to support energy requirements.

An additional observation is that neuronal ATP concentration can approach its reference level before the phosphocreatine pool has fully recovered.

This difference highlights the importance of considering both ATP availability and the recovery of intracellular energy-buffering reserves when interpreting the model's metabolic response.

The project also illustrates how mathematical modeling can be used to investigate relationships between ion homeostasis, energy demand, mitochondrial metabolism, and creatine kinase dynamics.

## Reproducibility

The project uses Python with NumPy, SciPy, pandas, and Matplotlib.

Simulation and analysis entry points include:

```bash
python run_experiments.py
python analyze_experiments.py
python analyze_atp_fluxes.py
python mitochondrial_time_course.py
python analyze_phosphocreatine.py
python compare_normalized_pcr_recovery.py
python calculate_pcr_rmsd.py
```

Some analysis scripts require simulation output files generated by earlier runs. Numerical settings, input file requirements, and execution order should be checked before attempting full reproduction.

## Limitations and Future Directions

The current work is based on a published mechanistic model and its specified parameterization.

Important limitations include numerical solver sensitivity, dependence on model assumptions, and the absence of independent experimental validation of the additional simulation results.

Future work could include parameter sensitivity analysis, investigation of metabolic perturbations, and comparisons with additional experimental measurements.


## Transient Cerebral Ischemia and Reperfusion Investigation

Building on the original model reproduction, I extended the Python implementation to investigate neuronal and astrocytic energy metabolism during transient cerebral ischemia and reperfusion.

### Research Objective

To investigate how the duration of reduced cerebral blood flow influences ATP availability, extracellular glucose and lactate dynamics, and neuronal–astrocytic lactate transport during ischemia and subsequent reperfusion.

### Simulation Protocol

- Baseline cerebral blood flow: 0.012.
- Ischemia onset: 100 seconds.
- Cerebral blood flow reduction: 80%.
- Ischemia durations investigated: 300, 600, 1000, 1500 and 2000 seconds.
- Reperfusion: restoration of baseline cerebral blood flow.
- Reperfusion observation period: 1000 seconds.

The 1000-second condition was included as an anchor for comparison with an earlier GAMMA-MRI modelling report. The additional durations were exploratory simulation conditions, not clinical thresholds.

### Metabolic Response Analysis

The simulations were used to examine:

- Neuronal and astrocytic ATP depletion and recovery.
- Extracellular glucose and lactate concentration dynamics.
- Capillary oxygen dynamics.
- Cellular lactate dehydrogenase reaction rates.
- Neuronal, astrocytic and vascular lactate transport.
- Changes in extracellular lactate balance during reperfusion.

Across the five simulated durations, neuronal ATP depletion at the end of ischemia ranged from approximately 96.68% to 99.32%, while astrocytic ATP depletion ranged from 79.07% to 93.81%.

After the fixed 1000-second reperfusion period, neuronal ATP approached its pre-ischemia baseline in all five simulations. Astrocytic ATP showed a modest overshoot relative to baseline in the model.

### Extracellular Lactate Mass-Balance Analysis

The extracellular lactate balance was decomposed into three transport contributions:

1. Neuronal-to-extracellular lactate transport.
2. Astrocytic-to-extracellular lactate transport.
3. Extracellular-to-capillary lactate transport.

The model's extracellular lactate balance is represented by:


The extracellular lactate balance is represented by:

$$
\frac{d[LAC]_e}{dt}
=
\frac{v_{LACne}}{r_{en}}
+
\frac{v_{LACge}}{r_{eg}}
-
v_{LACec}
$$


The analysis revealed:

- Astrocytic transport made a positive cumulative contribution to extracellular lactate across all five simulations.
- Neuronal transport exhibited transient positive contributions during reperfusion in the longer ischemia simulations, although its cumulative contribution over the complete reperfusion period remained negative.
- Extracellular-to-capillary transport contributed negatively to the extracellular lactate balance.
- The integrated transport contributions reproduced the simulated change in extracellular lactate concentration with very small numerical differences.

Astrocytic transport represented approximately 48.4–49.5% of the combined absolute integrated transport magnitude across the five conditions. These percentages describe relative transport magnitudes, not the proportion of total lactate production.

### Comparison with Experimental Studies

The model's glucose and lactate responses were compared descriptively with published transient cerebral ischemia studies, including:

- Frykholm et al. (2005), investigating transient middle cerebral artery occlusion in macaques.
- Havelund et al. (2019), investigating transient focal ischemia in rats.

The comparison was exploratory and was not intended as formal model validation. Differences in species, ischemia protocols, brain regions, sampling methods and observation periods limit direct quantitative comparisons.

### Figures and Results

- Figures 1–5: ATP, glucose, lactate, oxygen and metabolic time-course responses.
- Figure 6: Exploratory comparison with experimental glucose and lactate observations.
- Figure 7: Cellular lactate dehydrogenase fluxes.
- Figure 8: Lactate transport fluxes.
- Figure 9: Extracellular lactate balance.
- Figure 10: Numerical consistency of the extracellular lactate mass-balance analysis.

### Reproducibility

The ischemia–reperfusion investigation can be explored using:

```bash

# Run the original ischemia-reperfusion simulation
python simulate_ischemia.py

# Run ischemia simulations for five durations
python simulate_ischemia_sweep.py 300
python simulate_ischemia_sweep.py 600
python simulate_ischemia_sweep.py 1000
python simulate_ischemia_sweep.py 1500
python simulate_ischemia_sweep.py 2000

# Analyse ischemia and metabolic responses
python analyse_ischemia_sweep.py
python analyse_lactate_fluxes.py

# Analyse extracellular lactate balance
python analyse_lactate_mass_balance.py
python analyse_reperfusion_mass_balance.py
python analyse_reperfusion_contribution_peaks.py
python analyse_cumulative_lactate_contributions.py
python analyse_relative_lactate_contributions.py

```

The scripts use the model implementation and simulation output files in this repository. Some analyses require the corresponding simulation datasets to have been generated first.

### Limitations

These results are predictions from a mechanistic computational model. They do not independently establish neuronal survival, tissue protection or clinical outcomes.

The original model does not explicitly represent glycogen metabolism. Further extensions, including glycogen metabolism and the pentose phosphate pathway, would require additional model development and evaluation.

The current investigation is intended as a foundation for further comparison with experimental measurements and discussion of potential model extensions.

###  References

Frykholm, P., Hillered, L., Långström, B., Persson, L., Valtysson, J., & Enblad, P. (2005). Relationship between cerebral blood flow and oxygen metabolism, and extracellular glucose and lactate concentrations during middle cerebral artery occlusion and reperfusion: A microdialysis and positron emission tomography study in nonhuman primates. *Journal of Neurosurgery, 102*(6), 1076–1084. https://doi.org/10.3171/jns.2005.102.6.1076

Havelund, J. F., Nygaard, K. H., Nielsen, T. H., Nordström, C.-H., Poulsen, F. R., Færgeman, N. J., Forsse, A., & Gramsbergen, J. B. (2019). In vivo microdialysis of endogenous and ¹³C-labeled TCA metabolites in rat brain: Reversible and persistent effects of mitochondrial inhibition and transient cerebral ischemia. *Metabolites, 9*(10), 204. https://doi.org/10.3390/metabo9100204

Jolivet, R. (2023). *D5.1: Report on biophysical model* (Version 2.0). GAMMA-MRI, Horizon 2020 Research and Innovation Programme, Grant Agreement No. 964644. Submitted 4 December 2023.

## Original Publication and Attribution

Jolivet, R., Coggan, J. S., Allaman, I., & Magistretti, P. J. (2015). Multi-timescale modeling of activity-dependent metabolic coupling in the neuron-glia-vasculature ensemble. *PLOS Computational Biology, 11*(2), e1004036.

https://doi.org/10.1371/journal.pcbi.1004036

Attribution: The original mathematical model and its scientific framework were developed by Jolivet and colleagues. This repository documents my independent Python implementation, numerical comparisons, and subsequent computational analyses. The original MATLAB source code is not redistributed in this repository.

## Author
Kingsley Nnaemeka Ugochukwu
MSc Student, Biochemistry and Molecular Biology  
Obafemi Awolowo University, Ile-Ife, Nigeria

Research interests: Computational Biology, Systems Neuroscience, Brain Energy Metabolism, and Mathematical Modeling.

GitHub: https://github.com/KingsleyUgochukwu