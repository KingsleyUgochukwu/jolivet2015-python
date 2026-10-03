import numpy as np


def blood_flow(t, t1, tstim, cbf, F0, species, DCBF=0.0):
    """
    Blood-flow function translated from BloodFlow.m
    in the Jolivet et al. (2015) MATLAB model.

    Parameters
    ----------
    t : float
        Current simulation time.
    t1 : float
        Beginning of stimulation.
    tstim : float
        Duration of stimulation.
    cbf : str
        'on' or 'off'.
    F0 : float
        Baseline blood flow.
    species : str
        'rat' or 'human'.
    DCBF : float
        Relative change in cerebral blood flow.
        Used by the human model.
    """

    lag = 1.0

    # ---------------------------------------------
    # Rat
    # ---------------------------------------------
    if species == "rat":

        def builtin1(x):
            return 1.1 + 1.5 * (
                np.exp(-x / 5.0)
                - np.exp(-x / 2.0)
            )

        def builtin2(x):
            return 1.5 * (
                -1.0 / 5.0 * np.exp(-x / 5.0)
                + 1.0 / 2.0 * np.exp(-x / 2.0)
            )

        def builtin3(x, lag_value, taurise):
            return np.exp(
                (x - lag_value) / taurise
            )

        lastF0 = F0 * builtin1(tstim - lag)

        taurise = 0.1 / builtin2(0)

        if cbf == "on":

            if (t <= tstim + t1) and (t >= t1 + lag):

                y = F0 * builtin1(
                    t - t1 - lag
                )

            elif t > tstim + t1:

                y = (
                    F0
                    + (lastF0 - F0)
                    * np.exp(
                        -(t - tstim - t1) / 5.0
                    )
                )

            else:

                y = F0 * (
                    1
                    + 0.1
                    * builtin3(
                        t - t1,
                        lag,
                        taurise
                    )
                )

        else:
            y = F0

    # ---------------------------------------------
    # Human
    # ---------------------------------------------
    elif species == "human":

        rt = 30.0
        lag0 = 2.0
        lag1 = 10.0

        if cbf == "on":

            if (
                t <= rt + t1 + lag0
                and t >= t1 + lag0
            ):

                y = F0 * (
                    1
                    + DCBF
                    * (t - t1 - lag0) / rt
                )

            elif (
                t <= t1 + tstim + rt + lag1
                and t >= t1 + tstim + lag1
            ):

                y = F0 * (
                    1
                    + DCBF
                    - DCBF
                    * (t - t1 - tstim - lag1)
                    / rt
                )

            elif (
                t >= rt + t1 + lag0
                and t <= t1 + tstim + lag1
            ):

                y = (1 + DCBF) * F0

            else:
                y = F0

        else:
            y = F0

    else:
        raise ValueError(
            "Unknown species. Use 'rat' or 'human'."
        )

    return y
if __name__ == "__main__":

    test_flow = blood_flow(
        t=0,
        t1=10,
        tstim=20,
        cbf="off",
        F0=1.0,
        species="rat"
    )

    print("Test blood flow:", test_flow)