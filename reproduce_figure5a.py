import numpy as np
import matplotlib.pyplot as plt
from blood_flow import blood_flow

t = np.linspace(0, 100, 2001)

cbf = np.array([
    blood_flow(
        ti, t1=0, tstim=60,
        cbf="on", F0=0.012,
        species="rat", DCBF=0.4
    )
    for ti in t
])

fig, ax = plt.subplots(figsize=(8, 5))

ax.axvspan(0, 60, color="lightgray", alpha=0.5)
ax.plot(t, cbf, color="black", linewidth=2)

ax.set_xlim(0, 100)
ax.set_xlabel("Time (s)")
ax.set_ylabel("CBF (s^-1)")
ax.set_title("Figure 5A - Cerebral Blood Flow")
ax.grid(alpha=0.2)

fig.tight_layout()
fig.savefig("Figure_05A_Python.png", dpi=300)
plt.show()

