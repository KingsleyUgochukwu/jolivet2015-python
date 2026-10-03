import numpy as np


def excitatory_conductance(t, t1, tstim, Ne=240, LONG="off"):
    """
    Translation of MATLAB ExcitatoryConductance.m
    """

    if LONG == "off" and tstim == 20:
        finf = 3.5
        f0 = 32.0
        tau = 2.5

    elif LONG == "off" and tstim == 60:
        finf = 2.5
        f0 = 23.0
        tau = 2.5

    else:
        finf = 3.7
        f0 = finf
        tau = 2.5

    gbar = 7.8e-06

    ge = (
        Ne
        * gbar
        * (
            finf
            + (f0 - finf) * np.exp(-t / tau)
        )
    )

    if t1 <= t <= (tstim + t1):
        return ge

    return 0.0


def inhibitory_conductance(t, t1, tstim):
    """
    Translation of MATLAB InhibitoryConductance.m.

    The original 2015 implementation returns zero
    during and outside stimulation.
    """

    return 0.0


def mean_field_jstim(x):
    """
    Translation of MATLAB MeanFieldJstim.m
    """

    if x == 0:
        return 0.0

    elif 0 < x < 0.004:
        A = -0.00024
        B = 13.005

        return A + B * x

    elif 0.004 <= x < 0.008:
        A1 = 0.00043
        A2 = 0.30818
        x0 = 0.00696
        dx = 0.00182

        return A2 + (A1 - A2) / (
            1 + np.exp((x - x0) / dx)
        )

    else:
        A1 = -0.06444
        A2 = 4622.7114
        x0 = 2322.61814
        power = 0.77215

        return A2 + (A1 - A2) / (
            1 + (x / x0) ** power
        )
if __name__ == "__main__":

    # Before stimulation
    ge_before = excitatory_conductance(
        t=5,
        t1=10,
        tstim=20,
        Ne=240,
        LONG="off",
    )

    # During stimulation
    ge_during = excitatory_conductance(
        t=10,
        t1=10,
        tstim=20,
        Ne=240,
        LONG="off",
    )

    # Inhibitory conductance
    gi = inhibitory_conductance(
        t=10,
        t1=10,
        tstim=20,
    )

    print("Excitatory conductance before stimulation:", ge_before)
    print("Excitatory conductance during stimulation:", ge_during)
    print("Inhibitory conductance:", gi)