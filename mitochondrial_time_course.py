
import numpy as np
import matplotlib.pyplot as plt

from parameters import get_parameters

p = get_parameters()

plt.figure(figsize=(9, 5))

for Ne in [120, 180, 240, 300, 360]:

    data = np.load(f"stimulation_Ne_{Ne}.npz")

    time = data["time"]
    y = data["states"]

    # Extract neuronal state variables
    ATPn = y[:, 14]
    O2n = y[:, 18]
    NADHn = y[:, 31]

    # Calculate ADP using the original model equation
    qAK = 0.92
    A = 2.212

    ADPn = (
        ATPn / 2
        * (
            -qAK
            + np.sqrt(
                qAK**2
                + 4 * qAK * (A / ATPn - 1)
            )
        )
    )

    # Calculate mitochondrial oxidative metabolism
    vMitooutn = (
        p["VMaxMitooutn"]
        * O2n / (O2n + p["KO2Mito"])
        * ADPn / (ADPn + p["KmADPn"])
        * NADHn / (p["KmNADHn"] + NADHn)
    )

    # ATP production term in the neuronal ATP balance
    mito_ATP = 3.6 * vMitooutn

    plt.plot(time, mito_ATP, label=f"Ne = {Ne}")

    print(
        f"Ne={Ne}: "
        f"Baseline={mito_ATP[0]:.6f}, "
        f"Maximum={np.max(mito_ATP):.6f}, "
        f"Peak time={time[np.argmax(mito_ATP)]:.4f}s, "
        f"At 20s={mito_ATP[-1]:.6f}"
    )

plt.xlabel("Time (seconds)")
plt.ylabel("Mitochondrial ATP production term")
plt.title("Neuronal Mitochondrial ATP Production During Stimulation")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig("Mitochondrial_ATP_time_course.png", dpi=300)
plt.show()

print("\nSaved Mitochondrial_ATP_time_course.png")
