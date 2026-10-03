import numpy as np

from parameters import get_parameters
from rates import calculate_rates

def balance_equations(
    t,
    y,
    params=None,
    TIME=0.0,
    t1=0.0,
    tstim=20.0,
    CBF="off",
    F0=0.012,
    species="rat",
    DCBF=0.4,
    Ne=240,
    LONG="off",
):
    """
    Differential equations for the Jolivet et al. (2015) model.

    Translation of MATLAB balanceequ.m.
    """

    if params is None:
        params = get_parameters()

    r = calculate_rates(
    t=t,
    y=y,
    params=params,
    TIME=TIME,
    t1=t1,
    tstim=tstim,
    CBF=CBF,
    F0=F0,
    species=species,
    DCBF=DCBF,
    Ne=Ne,
    LONG=LONG,
)
    # -------------------------------------------------
    # AMP-ATP relationship
    # Original MATLAB: AMPATP.m
    # -------------------------------------------------

    un = (
        params["qAK"] ** 2
        + 4 * params["qAK"]
        * (params["A"] / y[14] - 1)
    )

    ug = (
        params["qAK"] ** 2
        + 4 * params["qAK"]
        * (params["A"] / y[15] - 1)
    )

    dAMPdATPn = (
        -1
        + params["qAK"] / 2
        - 0.5 * np.sqrt(un)
        + params["qAK"] * params["A"]
        / (y[14] * np.sqrt(un))
    )

    dAMPdATPg = (
        -1
        + params["qAK"] / 2
        - 0.5 * np.sqrt(ug)
        + params["qAK"] * params["A"]
        / (y[15] * np.sqrt(ug))
    )

    # 33 state-variable derivatives
    dydt = np.zeros(33)

    # -------------------------------------------------
    # States 1-2: Na+
    # -------------------------------------------------

    # MATLAB dY(1)
    dydt[0] = (
        r["vLeakNan"]
        - 3 * r["vPumpn"]
        + r["vnstim"]
    )

    # MATLAB dY(2)
    dydt[1] = (
        r["vLeakNag"]
        - 3 * r["vPumpg"]
        + r["vgstim"]
    )

    # -------------------------------------------------
    # States 3-4: Glucose
    # -------------------------------------------------

    # MATLAB dY(3)
    dydt[2] = (
        r["vGLCen"]
        - r["vHKPFKn"]
    )

    # MATLAB dY(4)
    dydt[3] = (
        r["vGLCcg"]
        + r["vGLCeg"]
        - r["vHKPFKg"]
    )

    # -------------------------------------------------
    # States 5-6: GAP
    # -------------------------------------------------

    # MATLAB dY(5)
    dydt[4] = (
        2 * r["vHKPFKn"]
        - r["vPGKn"]
    )

    # MATLAB dY(6)
    dydt[5] = (
        2 * r["vHKPFKg"]
        - r["vPGKg"]
    )

    # -------------------------------------------------
    # States 7-8: PEP
    # -------------------------------------------------

    # MATLAB dY(7)
    dydt[6] = (
        r["vPGKn"]
        - r["vPKn"]
    )

    # MATLAB dY(8)
    dydt[7] = (
        r["vPGKg"]
        - r["vPKg"]
    )
    # -------------------------------------------------
    # States 9-10: Pyruvate
    # -------------------------------------------------

    # MATLAB dY(9)
    dydt[8] = (
        r["vPKn"]
        - r["vLDHn"]
        - r["vMitoinn"]
    )

    # MATLAB dY(10)
    dydt[9] = (
        r["vPKg"]
        - r["vLDHg"]
        - r["vMitoing"]
    )

    # -------------------------------------------------
    # States 11-12: Lactate
    # -------------------------------------------------

    # MATLAB dY(11)
    dydt[10] = (
        r["vLDHn"]
        - r["vLACne"]
    )

    # MATLAB dY(12)
    dydt[11] = (
        r["vLDHg"]
        - r["vLACge"]
        - r["vLACgc"]
    )

    # -------------------------------------------------
    # States 13-14: Cytosolic NADH
    # -------------------------------------------------

    # MATLAB dY(13)
    dydt[12] = (
        1 / (1 - params["MVF"])
        * (
            r["vPGKn"]
            - r["vLDHn"]
            - r["vShuttlen"]
        )
    )

    # MATLAB dY(14)
    dydt[13] = (
        1 / (1 - params["MVF"])
        * (
            r["vPGKg"]
            - r["vLDHg"]
            - r["vShuttleg"]
        )
    )
    # -------------------------------------------------
    # States 15-16: ATP
    # -------------------------------------------------

    # MATLAB dY(15): neuronal ATP
    dydt[14] = (
        -2 * r["vHKPFKn"]
        + r["vPGKn"]
        + r["vPKn"]
        - params["vATPasesn"]
        - r["vPumpn"]
        + 3.6 * r["vMitooutn"]
        + r["vCKn"]
    ) / (1 - dAMPdATPn)

    # MATLAB dY(16): glial ATP
    dydt[15] = (
        -2 * r["vHKPFKg"]
        + r["vPGKg"]
        + r["vPKg"]
        - params["vATPasesg"]
        - (7 / 4) * r["vPumpg"]
        + (3 / 4) * params["vPumpg0"]
        + 3.6 * r["vMitooutg"]
        + r["vCKg"]
    ) / (1 - dAMPdATPg)

    # -------------------------------------------------
    # States 17-18: Phosphocreatine
    # -------------------------------------------------

    # MATLAB dY(17)
    dydt[16] = -r["vCKn"]

    # MATLAB dY(18)
    dydt[17] = -r["vCKg"]

    # -------------------------------------------------
    # States 19-20: Oxygen
    # -------------------------------------------------

    # MATLAB dY(19)
    dydt[18] = (
        r["vO2mcn"]
        - 0.6 * r["vMitooutn"]
    )

    # MATLAB dY(20)
    dydt[19] = (
        r["vO2mcg"]
        - 0.6 * r["vMitooutg"]
    )
    # -------------------------------------------------
    # States 21-25: Capillary / vascular compartment
    # Original MATLAB: balanceequ.m
    # -------------------------------------------------

    if CBF == "on":

        # MATLAB dY(21): capillary oxygen
        dydt[20] = (
            r["vO2c"]
            - (1 / params["rcn"]) * r["vO2mcn"]
            - (1 / params["rcg"]) * r["vO2mcg"]
        )

        # MATLAB dY(22): capillary glucose
        dydt[21] = (
            r["vGLCc"]
            - (1 / params["rce"]) * r["vGLCce"]
            - (1 / params["rcg"]) * r["vGLCcg"]
        )

        # MATLAB dY(23): capillary lactate
        dydt[22] = (
            r["vLACc"]
            + (1 / params["rce"]) * r["vLACec"]
            + (1 / params["rcg"]) * r["vLACgc"]
        )

        # MATLAB dY(24): venous blood volume
        dydt[23] = (
            r["BF"]
            - params["F0"]
            * (y[23] / params["Vv0"]) ** (1 / params["alphav"])
        ) / (
            1
            + params["F0"]
            * params["tauv"]
            / params["Vv0"]
            * (y[23] / params["Vv0"]) ** (-1 / 2)
        )

        # Venous outflow
        Fout = (
            params["F0"]
            * (
                (y[23] / params["Vv0"]) ** (1 / params["alphav"])
                + dydt[23]
                * params["tauv"]
                / params["Vv0"]
                * (y[23] / params["Vv0"]) ** (-1 / 2)
            )
        )

        # MATLAB dY(25)
        dydt[24] = (
            r["BF"] * (params["O2a"] - r["O2cbar"])
            - Fout * y[24] / y[23]
        )

    else:

        # Original MATLAB CBF-off branch
        dydt[20] = 0.0
        dydt[21] = 0.0
        dydt[22] = 0.0

        dydt[23] = (
            r["BF"]
            - params["F0"]
            * (y[23] / params["Vv0"]) ** (1 / params["alphav"])
        ) / (
            1
            + params["F0"]
            * params["tauv"]
            / params["Vv0"]
            * (y[23] / params["Vv0"]) ** (-1 / 2)
        )

        Fout = 0.0
        dydt[24] = 0.0

    # -------------------------------------------------
    # States 26-27: Extracellular compartment
    # -------------------------------------------------

    # MATLAB dY(26): extracellular glucose
    dydt[25] = (
        r["vGLCce"]
        - (1 / params["reg"]) * r["vGLCeg"]
        - (1 / params["ren"]) * r["vGLCen"]
    )

    # MATLAB dY(27): extracellular lactate
    dydt[26] = (
        (1 / params["ren"]) * r["vLACne"]
        + (1 / params["reg"]) * r["vLACge"]
        - r["vLACec"]
    )
    # -------------------------------------------------
    # States 28-31: Hodgkin-Huxley electrophysiology
    # Original MATLAB: balanceequ.m
    # -------------------------------------------------

    if LONG == "on":

        # MATLAB dY(28)-dY(31)
        dydt[27] = 0.0
        dydt[28] = 0.0
        dydt[29] = 0.0
        dydt[30] = 0.0

    else:

        # MATLAB dY(28): neuronal membrane potential
        dydt[27] = (
            -r["IL"]
            - r["INa"]
            - r["IK"]
            - r["ICa"]
            - r["ImAHP"]
            - r["dIPump"]
            + r["Isyne"]
            + r["Isyni"]
        ) / params["Cm"]

        # MATLAB dY(29): h gating variable
        dydt[28] = (
            params["phih"]
            * (r["hinf"] - y[28])
            / r["tauh"]
        )

        # MATLAB dY(30): n gating variable
        dydt[29] = (
            params["phin"]
            * (r["ninf"] - y[29])
            / r["taun"]
        )

        # MATLAB dY(31): intracellular calcium
        dydt[30] = (
            -params["SmVn"] / params["F"] * r["ICa"]
            - (y[30] - params["Ca0"]) / params["tauCa"]
        )

    # -------------------------------------------------
    # States 32-33: Mitochondrial NADH
    # -------------------------------------------------

    # MATLAB dY(32): neuronal mitochondrial NADH
    dydt[31] = (
        1 / params["MVF"]
        * (
            4 * r["vMitoinn"]
            - r["vMitooutn"]
            + r["vShuttlen"]
        )
    )

    # MATLAB dY(33): glial mitochondrial NADH
    dydt[32] = (
        1 / params["MVF"]
        * (
            4 * r["vMitoing"]
            - r["vMitooutg"]
            + r["vShuttleg"]
        )
    )
    return dydt
if __name__ == "__main__":


    from initial_state import initial_state

    params = get_parameters()

    dydt = balance_equations(
        t=0.0,
        y=initial_state,
        params=params,
    )

    print("Number of derivatives:", len(dydt))

    print("dY1  neuronal Na+:", dydt[0])
    print("dY2  glial Na+:", dydt[1])
    print("dY3  neuronal glucose:", dydt[2])
    print("dY4  glial glucose:", dydt[3])
    print("dY5  neuronal GAP:", dydt[4])
    print("dY6  glial GAP:", dydt[5])
    print("dY7  neuronal PEP:", dydt[6])
    print("dY8  glial PEP:", dydt[7])
    print("dY9  neuronal pyruvate:", dydt[8])
    print("dY10 glial pyruvate:", dydt[9])
    print("dY11 neuronal lactate:", dydt[10])
    print("dY12 glial lactate:", dydt[11])
    print("dY13 neuronal cytosolic NADH:", dydt[12])
    print("dY14 glial cytosolic NADH:", dydt[13])
    print("dY15 neuronal ATP:", dydt[14])
    print("dY16 glial ATP:", dydt[15])
    print("dY17 neuronal phosphocreatine:", dydt[16])
    print("dY18 glial phosphocreatine:", dydt[17])
    print("dY19 neuronal oxygen:", dydt[18])
    print("dY20 glial oxygen:", dydt[19])
    print("dY21 capillary oxygen:", dydt[20])
    print("dY22 capillary glucose:", dydt[21])
    print("dY23 capillary lactate:", dydt[22])
    print("dY24 venous blood volume:", dydt[23])
    print("dY25 deoxyhemoglobin:", dydt[24])
    print("dY26 extracellular glucose:", dydt[25])
    print("dY27 extracellular lactate:", dydt[26])
    print("dY28 membrane potential:", dydt[27])
    print("dY29 h gating variable:", dydt[28])
    print("dY30 n gating variable:", dydt[29])
    print("dY31 intracellular calcium:", dydt[30])
    print("dY32 neuronal mitochondrial NADH:", dydt[31])
    print("dY33 glial mitochondrial NADH:", dydt[32])