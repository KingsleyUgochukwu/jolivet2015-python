
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(9, 5))

print("\nNEURONAL PHOSPHOCREATINE ANALYSIS")
print("-" * 65)

for Ne in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne}.npz")

    time = data["time"]
    states = data["states"]

    # Neuronal phosphocreatine (PCr)
    PCr = states[:, 16]

    initial_PCr = PCr[0]
    minimum_PCr = np.min(PCr)

    depletion_percent = (
        (initial_PCr - minimum_PCr)
        / initial_PCr
    ) * 100

    print(
        f"Ne={Ne}: "
        f"Initial PCr={initial_PCr:.6f}, "
        f"Minimum PCr={minimum_PCr:.6f}, "
        f"PCr depletion={depletion_percent:.2f}%"
    )

    plt.plot(time, PCr, label=f"Ne = {Ne}")

plt.xlabel("Time (seconds)")
plt.ylabel("Neuronal phosphocreatine concentration")
plt.title("Phosphocreatine Dynamics During Neuronal Stimulation")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig("Phosphocreatine_time_course.png", dpi=300)
plt.show()

print("\nSaved Phosphocreatine_time_course.png")
