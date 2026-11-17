
import pandas as pd
import matplotlib.pyplot as plt

# Load our model-derived summary
model = pd.read_csv("metabolite_response_summary.csv")

# Calculate representative model responses
model_glucose = model["glucose_minimum_percent_baseline"].mean()
model_lactate = model["lactate_maximum_percent_baseline"].mean()

# Approximate experimental reference values
studies = ["Our Python model", "Havelund (2019)", "Frykholm (2005)"]

glucose = [model_glucose, 50.0, 5.26]
lactate = [model_lactate, 250.0, 270.0]

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Glucose comparison
axes[0].bar(studies, glucose)
axes[0].set_title("Minimum Extracellular Glucose")
axes[0].set_ylabel("Glucose (% of baseline)")
axes[0].set_ylim(0, 110)
axes[0].tick_params(axis="x", rotation=20)
axes[0].grid(axis="y", alpha=0.3)

for i, value in enumerate(glucose):
    label = "<50%" if i == 1 else f"{value:.1f}%"
    axes[0].text(i, value + 2, label, ha="center")

# Lactate comparison
axes[1].bar(studies, lactate)
axes[1].set_title("Peak Extracellular Lactate")
axes[1].set_ylabel("Lactate (% of baseline)")
axes[1].tick_params(axis="x", rotation=20)
axes[1].grid(axis="y", alpha=0.3)

for i, value in enumerate(lactate):
    axes[1].text(i, value + 8, f"{value:.1f}%", ha="center")

fig.suptitle(
    "Comparison of Model Predictions and Experimental Findings",
    fontsize=14
)

plt.tight_layout()

plt.savefig(
    "figure6_model_experimental_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figure 6 saved successfully.")
print(f"Mean model glucose minimum: {model_glucose:.2f}%")
print(f"Mean model lactate peak: {model_lactate:.2f}%")
