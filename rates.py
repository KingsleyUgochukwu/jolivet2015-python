import numpy as np

from parameters import get_parameters
from blood_flow import blood_flow
from conductance import (
    excitatory_conductance,
    inhibitory_conductance,
    mean_field_jstim,
)

def calculate_rates(
    t,
    y,
    params=None,
    TIME=0.0,
    t1=10.0,
    tstim=20.0,
    CBF="off",
    F0=0.012,
    species="rat",
    DCBF=0.0,
    Ne=240,
    LONG="off",
):
    """
    Calculate metabolic rates for the Jolivet et al. (2015) model.

    Parameters
    ----------
    t : float
        Simulation time.
    y : array-like
        Vector containing the 33 model state variables.
    params : dict, optional
        Model parameters.

    Returns
    -------
    rates : dict
        Dictionary containing calculated reaction/transport rates.
    """

    if params is None:
        params = get_parameters()

    p = params
    r = {}

    # Baseline neuronal sodium concentration
    # Original MATLAB: Na0 = Y0(1)
    Na0 = y[0]

    r["vLeakNan"] = (
        p["SmVn"] * p["gNan"] / p["F"]
        * (p["RTF"] * np.log(p["Nae"] / y[0]) - y[27])
    )

    r["vLeakNag"] = (
        p["SmVg"] * p["gNag"] / p["F"]
        * (p["RTF"] * np.log(p["Nae"] / y[1]) - p["Vm"])
    )

    # -------------------------------------------------
    # Na+/K+-ATPase
    # -------------------------------------------------

    r["vPumpn"] = (
        p["SmVn"] * p["kPumpn"] * y[14] * y[0]
        / (1 + y[14] / p["KmPump"])
    )

    r["vPumpg"] = (
        p["SmVg"] * p["kPumpg"] * y[15] * y[1]
        / (1 + y[15] / p["KmPump"])
    )

    # -------------------------------------------------
    # Glucose transport
    # -------------------------------------------------

    r["vGLCen"] = p["TMaxGLCen"] * (
        y[25] / (y[25] + p["KtGLCen"])
        - y[2] / (y[2] + p["KtGLCen"])
    )

    r["vGLCeg"] = p["TMaxGLCeg"] * (
        y[25] / (y[25] + p["KtGLCeg"])
        - y[3] / (y[3] + p["KtGLCeg"])
    )

    r["vGLCcg"] = p["TMaxGLCcg"] * (
        y[21] / (y[21] + p["KtGLCcg"])
        - y[3] / (y[3] + p["KtGLCcg"])
    )

    r["vGLCce"] = p["TMaxGLCce"] * (
        y[21] / (y[21] + p["KtGLCce"])
        - y[25] / (y[25] + p["KtGLCce"])
    )

    # -------------------------------------------------
    # Hexokinase / phosphofructokinase
    # -------------------------------------------------

    r["vHKPFKn"] = (
        p["kHKPFKn"]
        * y[14]
        * y[2] / (y[2] + p["Kg"])
        * 1 / (1 + (y[14] / p["KIATP"]) ** p["nH"])
    )

    r["vHKPFKg"] = (
        p["kHKPFKg"]
        * y[15]
        * y[3] / (y[3] + p["Kg"])
        * 1 / (1 + (y[15] / p["KIATP"]) ** p["nH"])
    )
        # -------------------------------------------------
    # ADP concentrations
    # Original MATLAB: Scripts/ADP.m
    # -------------------------------------------------

    qAK = 0.92
    A = 2.212

    ADPn = (
        y[14] / 2
        * (
            -qAK
            + np.sqrt(
                qAK * qAK
                + 4 * qAK * (A / y[14] - 1)
            )
        )
    )

    ADPg = (
        y[15] / 2
        * (
            -qAK
            + np.sqrt(
                qAK * qAK
                + 4 * qAK * (A / y[15] - 1)
            )
        )
    )

    r["ADPn"] = ADPn
    r["ADPg"] = ADPg

    # -------------------------------------------------
    # Phosphoglycerate kinase (PGK)
    # -------------------------------------------------

    r["vPGKn"] = (
        p["kPGKn"]
        * y[4]
        * ADPn
        * (p["N"] - y[12])
        / y[12]
    )

    r["vPGKg"] = (
        p["kPGKg"]
        * y[5]
        * ADPg
        * (p["N"] - y[13])
        / y[13]
    )

    # -------------------------------------------------
    # Pyruvate kinase (PK)
    # -------------------------------------------------

    r["vPKn"] = p["kPKn"] * y[6] * ADPn
    r["vPKg"] = p["kPKg"] * y[7] * ADPg

    # -------------------------------------------------
    # Lactate dehydrogenase (LDH)
    # -------------------------------------------------

    r["vLDHn"] = (
        p["kLDHnps"] * y[8] * y[12]
        - p["kLDHnms"] * y[10] * (p["N"] - y[12])
    )

    r["vLDHg"] = (
        p["kLDHgps"] * y[9] * y[13]
        - p["kLDHgms"] * y[11] * (p["N"] - y[13])
    )

    # -------------------------------------------------
    # Lactate transport
    # -------------------------------------------------

    r["vLACne"] = p["TMaxLACne"] * (
        y[10] / (y[10] + p["KtLACne"])
        - y[26] / (y[26] + p["KtLACne"])
    )

    r["vLACge"] = p["TMaxLACge"] * (
        y[11] / (y[11] + p["KtLACge"])
        - y[26] / (y[26] + p["KtLACge"])
    )

    r["vLACgc"] = p["TMaxLACgc"] * (
        y[11] / (y[11] + p["KtLACgc"])
        - y[22] / (y[22] + p["KtLACgc"])
    )

    r["vLACec"] = p["TMaxLACec"] * (
        y[26] / (y[26] + p["KtLACec"])
        - y[22] / (y[22] + p["KtLACec"])
    )
        # -------------------------------------------------
    # Mitochondrial metabolism
    # Original MATLAB: rates.m
    # -------------------------------------------------

    # Pyruvate entry into neuronal mitochondria
    r["vMitoinn"] = (
        p["VMaxMitoinn"]
        * y[8] / (y[8] + p["KmMito"])
        * (p["N"] - y[31])
        / (p["N"] - y[31] + p["KmNADn"])
    )

    # Neuronal mitochondrial oxidative metabolism
    r["vMitooutn"] = (
        p["VMaxMitooutn"]
        * y[18] / (y[18] + p["KO2Mito"])
        * ADPn / (ADPn + p["KmADPn"])
        * y[31] / (p["KmNADHn"] + y[31])
    )

    # Pyruvate entry into glial mitochondria
    r["vMitoing"] = (
        p["VMaxMitoing"]
        * y[9] / (y[9] + p["KmMito"])
        * (p["N"] - y[32])
        / (p["N"] - y[32] + p["KmNADg"])
    )

    # Glial mitochondrial oxidative metabolism
    r["vMitooutg"] = (
        p["VMaxMitooutg"]
        * y[19] / (y[19] + p["KO2Mito"])
        * ADPg / (ADPg + p["KmADPg"])
        * y[32] / (p["KmNADHg"] + y[32])
    )

    # -------------------------------------------------
    # Cytosolic-mitochondrial NADH shuttle
    # -------------------------------------------------

    Rn_minus = y[12] / (p["N"] - y[12])
    Rn_plus = (p["N"] - y[31]) / y[31]

    Rg_minus = y[13] / (p["N"] - y[13])
    Rg_plus = (p["N"] - y[32]) / y[32]

    r["Rn_minus"] = Rn_minus
    r["Rn_plus"] = Rn_plus
    r["Rg_minus"] = Rg_minus
    r["Rg_plus"] = Rg_plus

    r["vShuttlen"] = (
        p["TnNADH"]
        * Rn_minus / (p["MnCyto"] + Rn_minus)
        * Rn_plus / (p["MnMito"] + Rn_plus)
    )

    r["vShuttleg"] = (
        p["TgNADH"]
        * Rg_minus / (p["MgCyto"] + Rg_minus)
        * Rg_plus / (p["MgMito"] + Rg_plus)
    )

    # -------------------------------------------------
    # Creatine kinase
    # -------------------------------------------------

    r["vCKn"] = (
        p["kCKnps"] * y[16] * ADPn
        - p["kCKnms"] * (p["C"] - y[16]) * y[14]
    )

    r["vCKg"] = (
        p["kCKgps"] * y[17] * ADPg
        - p["kCKgms"] * (p["C"] - y[17]) * y[15]
    )
        # -------------------------------------------------
    # Oxygen exchange
    # Original MATLAB: rates.m
    # -------------------------------------------------

    r["vO2mcn"] = (
        p["PScapVn"]
        * (
            p["KO2"]
            * (p["HbOP"] / y[20] - 1) ** (-1 / p["nh"])
            - y[18]
        )
    )

    r["vO2mcg"] = (
        p["PScapVg"]
        * (
            p["KO2"]
            * (p["HbOP"] / y[20] - 1) ** (-1 / p["nh"])
            - y[19]
        )
    )
        # -------------------------------------------------
    # Cerebral blood flow
    # -------------------------------------------------

    BF = blood_flow(
        t=t + TIME,
        t1=t1,
        tstim=tstim,
        cbf=CBF,
        F0=F0,
        species=species,
        DCBF=DCBF,
    )

    r["BF"] = BF

    # -------------------------------------------------
    # Capillary exchange
    # -------------------------------------------------

    r["vO2c"] = (
        2 * BF / p["Vcap"]
        * (p["O2a"] - y[20])
    )

    r["vGLCc"] = (
        2 * BF / p["Vcap"]
        * (p["GLCa"] - y[21])
    )

    r["vLACc"] = (
        2 * BF / p["Vcap"]
        * (p["LACa"] - y[22])
    )

    r["O2cbar"] = 2 * y[20] - p["O2a"]
    # -------------------------------------------------
    # Hodgkin-Huxley gating kinetics
    # Original MATLAB: rates.m
    # -------------------------------------------------

    Vm = y[27]   # MATLAB Y(28)

    r["alpham"] = (
        -0.1 * (Vm + 33)
        / (np.exp(-0.1 * (Vm + 33)) - 1)
    )

    r["betam"] = (
        4 * np.exp(-(Vm + 58) / 12)
    )

    r["alphah"] = (
        0.07 * np.exp(-(Vm + 50) / 10)
    )

    r["betah"] = (
        1 / (np.exp(-0.1 * (Vm + 20)) + 1)
    )

    r["alphan"] = (
        -0.01 * (Vm + 34)
        / (np.exp(-0.1 * (Vm + 34)) - 1)
    )

    r["betan"] = (
        0.125 * np.exp(-(Vm + 44) / 25)
    )

    # Steady-state gating values
    r["minf"] = r["alpham"] / (r["alpham"] + r["betam"])
    r["ninf"] = r["alphan"] / (r["alphan"] + r["betan"])
    r["hinf"] = r["alphah"] / (r["alphah"] + r["betah"])

    # Gating time constants
    r["taun"] = (
        1 / (r["alphan"] + r["betan"]) * 1e-03
    )

    r["tauh"] = (
        1 / (r["alphah"] + r["betah"]) * 1e-03
    )
       # -------------------------------------------------
    # Neuronal ionic currents
    # Original MATLAB: rates.m
    # -------------------------------------------------

    # Leak reversal potential
    r["EL"] = (
        p["gKpas"] * p["EK"] / (p["gKpas"] + p["gNan"])
        + p["gNan"] / (p["gKpas"] + p["gNan"])
        * p["RTF"] * np.log(p["Nae"] / y[0])
    )

    # Leak current
    r["IL"] = p["gL"] * (Vm - r["EL"])

    # Sodium current
    r["INa"] = (
        p["gNa"]
        * r["minf"] ** 3
        * y[28]
        * (
            Vm
            - p["RTF"] * np.log(p["Nae"] / y[0])
        )
    )

    # Potassium current
    r["IK"] = (
        p["gK"]
        * y[29] ** 4
        * (Vm - p["EK"])
    )

    # Calcium activation
    r["mCa"] = (
        1 / (1 + np.exp(-(Vm + 20) / 9))
    )

    # Calcium current
    r["ICa"] = (
        p["gCa"]
        * r["mCa"] ** 2
        * (Vm - p["ECa"])
    )

    # Calcium-dependent after-hyperpolarization current
    r["ImAHP"] = (
        p["gmAHP"]
        * y[30] / (y[30] + p["KD"])
        * (Vm - p["EK"])
    )

    # Pump-related current
    r["dIPump"] = (
        p["F"]
        * p["kPumpn"]
        * y[14]
        * (y[0] - Na0)
        / (1 + y[14] / p["KmPump"])
    )
    # -------------------------------------------------
    # Synaptic currents and stimulation
    # Original MATLAB: rates.m
    # -------------------------------------------------

    ge = excitatory_conductance(
        t=t + TIME,
        t1=t1,
        tstim=tstim,
        Ne=Ne,
        LONG=LONG,
    )

    gi = inhibitory_conductance(
        t=t + TIME,
        t1=t1,
        tstim=tstim,
    )

    r["ge"] = ge
    r["gi"] = gi

    # Excitatory synaptic current
    r["Isyne"] = (
        -ge * (Vm - p["Ee"])
    )

    # Inhibitory synaptic current
    r["Isyni"] = (
        -gi * (Vm - p["Ei"])
    )

    # Neuronal stimulation rate
    if LONG == "off":

        r["vnstim"] = (
            p["SmVn"] / p["F"]
            * (
                (2 / 3) * r["Isyne"]
                - r["INa"]
            )
        )

    else:

        r["vnstim"] = mean_field_jstim(ge)

    # Glial stimulation rate
    r["vgstim"] = (
        p["SmVg"] / p["F"]
        * (2 / 3)
        * p["glia"]
        * ge
    )
    return r
if __name__ == "__main__":

    from initial_state import initial_state

    params = get_parameters()
    rates = calculate_rates(
    t=0.0,
    y=initial_state,
    params=params,
    TIME=0.0,
    t1=0.0,
    tstim=20.0,
    CBF="off",
    F0=0.012,
    species="rat",
    DCBF=0.4,
    Ne=240,
    LONG="off",
)

    print("Number of rates calculated:", len(rates))
    print("Neuronal Na leak:", rates["vLeakNan"])
    print("Neuronal Na/K pump:", rates["vPumpn"])
    print("Neuronal glucose transport:", rates["vGLCen"])
    print("Neuronal HK/PFK:", rates["vHKPFKn"])
    print("Neuronal ADP:", rates["ADPn"])
    print("Glial ADP:", rates["ADPg"])
    print("Neuronal PGK:", rates["vPGKn"])
    print("Neuronal PK:", rates["vPKn"])
    print("Neuronal LDH:", rates["vLDHn"])
    print("Neuron-to-extracellular lactate:", rates["vLACne"])
    print("Neuronal mitochondrial input:", rates["vMitoinn"])
    print("Neuronal mitochondrial output:", rates["vMitooutn"])
    print("Neuronal NADH shuttle:", rates["vShuttlen"])
    print("Neuronal creatine kinase:", rates["vCKn"])
    print("Neuron-capillary oxygen flux:", rates["vO2mcn"])
    print("Glia-capillary oxygen flux:", rates["vO2mcg"])
    print("Blood flow:", rates["BF"])
    print("Capillary oxygen supply:", rates["vO2c"])
    print("Capillary glucose supply:", rates["vGLCc"])
    print("Capillary lactate exchange:", rates["vLACc"])
    print("Mean capillary oxygen:", rates["O2cbar"])
    print("Alpha m:", rates["alpham"])
    print("Beta m:", rates["betam"])
    print("Alpha h:", rates["alphah"])
    print("Beta h:", rates["betah"])
    print("Alpha n:", rates["alphan"])
    print("Beta n:", rates["betan"])
    print("m infinity:", rates["minf"])
    print("h infinity:", rates["hinf"])
    print("n infinity:", rates["ninf"])
    print("Tau h:", rates["tauh"])
    print("Tau n:", rates["taun"])
    print("Leak reversal potential:", rates["EL"])
    print("Leak current:", rates["IL"])
    print("Sodium current:", rates["INa"])
    print("Potassium current:", rates["IK"])
    print("Calcium activation mCa:", rates["mCa"])
    print("Calcium current:", rates["ICa"])
    print("mAHP current:", rates["ImAHP"])
    print("Pump-related current:", rates["dIPump"])
    print("Excitatory conductance:", rates["ge"])
    print("Inhibitory conductance:", rates["gi"])
    print("Excitatory synaptic current:", rates["Isyne"])
    print("Inhibitory synaptic current:", rates["Isyni"])
    print("Neuronal stimulation rate:", rates["vnstim"])
    print("Glial stimulation rate:", rates["vgstim"])
    print("\n--- ELECTROPHYSIOLOGY DIAGNOSTIC ---")

    print("IL:", rates["IL"])
    print("INa:", rates["INa"])
    print("IK:", rates["IK"])
    print("ICa:", rates["ICa"])
    print("ImAHP:", rates["ImAHP"])
    print("dIPump:", rates["dIPump"])
    print("Isyne:", rates["Isyne"])
    print("Isyni:", rates["Isyni"])

    current_residual = (
        -rates["IL"]
        - rates["INa"]
        - rates["IK"]
        - rates["ICa"]
        - rates["ImAHP"]
        - rates["dIPump"]
        + rates["Isyne"]
        + rates["Isyni"]
    )

    print("Current residual:", current_residual)
    print("Cm:", params["Cm"])
    print("Predicted dVm/dt:", current_residual / params["Cm"])