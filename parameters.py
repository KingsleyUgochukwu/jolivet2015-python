# Parameters for the Jolivet et al. (2015) model
# Translated from MATLAB: Parameters/parameters.m


def get_parameters():
    """
    Return model parameters in a Python dictionary.
    """

    p = {}


    # -------------------------------------------------
    # Compartment volumes (per unit tissue volume)
    # -------------------------------------------------

    p["Ve"] = 0.2
    p["Vcap"] = 0.0055
    p["Vg"] = 0.25
    p["Vn"] = 0.45
    p["MVF"] = 0.07

    # Volume ratios
    p["ren"] = p["Ve"] / p["Vn"]
    p["reg"] = p["Ve"] / p["Vg"]
    p["rce"] = p["Vcap"] / p["Ve"]
    p["rcn"] = p["Vcap"] / p["Vn"]
    p["rcg"] = p["Vcap"] / p["Vg"]

    # -------------------------------------------------
    # Surface-to-volume ratios and sodium conductances
    # -------------------------------------------------

    p["SmVn"] = 2.5e04       # cm^-1
    p["SmVg"] = 2.5e04       # cm^-1
    p["gNan"] = 0.0155       # mS cm^-2
    p["gNag"] = 0.00623      # mS cm^-2
    p["gKpas"] = 0.232       # mS cm^-2

    # -------------------------------------------------
    # Physical constants and parameters
    # -------------------------------------------------

    p["R"] = 8.314510        # J mol^-1 K^-1
    p["F"] = 9.64853e04      # C mol^-1
    p["RTF"] = 26.73         # mV, corresponding to 37 C
    p["Nae"] = 150           # mmol/L
    p["Vm"] = -70            # mV

    # -------------------------------------------------
    # Na+/K+-ATPase
    # -------------------------------------------------

    p["kPumpn"] = 2.49e-06
    p["kPumpg"] = 4.64e-07
    p["KmPump"] = 0.5
    p["vPumpg0"] = 0.0708
        # -------------------------------------------------
    # Glucose exchange constants
    # -------------------------------------------------

    p["TMaxGLCen"] = 0.041
    p["TMaxGLCce"] = 0.239
    p["TMaxGLCeg"] = 0.147
    p["TMaxGLCcg"] = 0.00164

    p["KtGLCen"] = 8
    p["KtGLCeg"] = 8
    p["KtGLCcg"] = 8
    p["KtGLCce"] = 8

    # -------------------------------------------------
    # Hexokinase-phosphofructokinase system
    # -------------------------------------------------

    p["kHKPFKn"] = 0.050435
    p["kHKPFKg"] = 0.185
    p["KIATP"] = 1
    p["nH"] = 4
    p["Kg"] = 0.05

    # -------------------------------------------------
    # Phosphoglycerate kinase
    # -------------------------------------------------

    p["kPGKn"] = 3.97
    p["kPGKg"] = 401.7
    p["N"] = 0.212

    # -------------------------------------------------
    # Pyruvate kinase
    # -------------------------------------------------

    p["kPKn"] = 36.7
    p["kPKg"] = 135.2
        # -------------------------------------------------
    # Lactate exchange constants
    # -------------------------------------------------

    # FourCIN is a simulation-level multiplier in the
    # original MATLAB model. Baseline value = 1.
    p["FourCIN"] = 1.0

    p["TMaxLACne"] = p["FourCIN"] * 23.5
    p["TMaxLACge"] = 107
    p["TMaxLACec"] = 0.30
    p["TMaxLACgc"] = 2.43e-03

    p["KtLACne"] = 0.74
    p["KtLACge"] = 3.5
    p["KtLACgc"] = 1.0
    p["KtLACec"] = 1.0

    # -------------------------------------------------
    # Lactate dehydrogenase (LDH)
    # -------------------------------------------------

    p["kLDHnps"] = 78.1
    p["kLDHgps"] = 1.71
    p["kLDHnms"] = 0.768
    p["kLDHgms"] = 0.099

    # -------------------------------------------------
    # NADH shuttles
    # -------------------------------------------------

    p["TnNADH"] = 9446
    p["TgNADH"] = 134.2

    p["MnCyto"] = 4.653e-08
    p["MgCyto"] = 2.614e-04
    p["MnMito"] = 3.666e05
    p["MgMito"] = 9.620e03

    p["KmNADn"] = 4.85e-01
    p["KmNADg"] = 38.24
    p["KmNADHn"] = 4.54e-02
    p["KmNADHg"] = 2.66e-02
        # -------------------------------------------------
    # Constant ATP consumption rates
    # -------------------------------------------------

    p["vATPasesn"] = 0.1388
    p["vATPasesg"] = 0.1351

    # -------------------------------------------------
    # Mitochondrial respiration
    # -------------------------------------------------

    p["KmMito"] = 0.04
    p["KO2Mito"] = 0.001

    p["KmADPn"] = 3.328e-03
    p["KmADPg"] = 4.989e-04

    p["VMaxMitooutn"] = 0.1610
    p["VMaxMitooutg"] = 0.0627

    p["VMaxMitoinn"] = 0.147
    p["VMaxMitoing"] = 5.31

    # -------------------------------------------------
    # Creatine kinase
    # -------------------------------------------------

    p["kCKnps"] = 8.85e-02
    p["kCKgps"] = 1.0e-03

    p["kCKnms"] = 0.00057
    p["kCKgms"] = 0.000007

    p["C"] = 10
        # -------------------------------------------------
    # Oxygen exchange constants
    # -------------------------------------------------

    p["PScapVn"] = 1.6608
    p["PScapVg"] = 0.8736
    p["KO2"] = 0.0361
    p["HbOP"] = 8.6
    p["nh"] = 2.73

    # -------------------------------------------------
    # Blood flow contributions
    # -------------------------------------------------

    p["O2a"] = 8.35
    p["GLCa"] = 4.75
    p["LACa"] = 0.498

    # -------------------------------------------------
    # Venous flow
    # -------------------------------------------------

    p["tauv"] = 35
    p["alphav"] = 0.5
    # -------------------------------------------------
    # Ratio of excitatory conductance
    # -------------------------------------------------

    p["glia"] = 50

    # -------------------------------------------------
    # Neuronal electrophysiology parameters
    # -------------------------------------------------

    p["Cm"] = 1e-03
    p["gL"] = 0.02
    p["gNa"] = 40
    p["gK"] = 18
    p["gCa"] = 0.02
    p["gmAHP"] = 6.5

    # Calcium-dependent parameters
    p["KD"] = 30e-03
    p["tauCa"] = 150e-03
    p["Ca0"] = 0.5e-04

    # Reversal potentials
    p["EK"] = -80
    p["ECa"] = 120
    p["Ee"] = 0
    p["Ei"] = -80

    # Gating kinetics
    p["phih"] = 4
    p["phin"] = 4
      # -------------------------------------------------
    # Stabilized coefficients
    # Original MATLAB: NewCoeffsStab.m
    #
    # IMPORTANT:
    # These values overwrite selected parameters above,
    # exactly as NewCoeffsStab.m does in the MATLAB model.
    # -------------------------------------------------

    p["kCKnps"] = 4.33e-02
    p["kCKgps"] = 1.35e-03

    p["KmNADn"] = 4.09e-01
    p["KmNADg"] = 4.03e01
    p["KmNADHn"] = 4.44e-02
    p["KmNADHg"] = 2.69e-02

    p["kLDHnps"] = 7.23e01
    p["kLDHgps"] = 1.59

    p["MnCyto"] = 4.9e-08
    p["MgCyto"] = 2.5e-04
    p["MnMito"] = 3.93e05
    p["MgMito"] = 1.06e04

    p["KmADPn"] = 3.41e-03
    p["KmADPg"] = 4.83e-04

    p["TMaxLACgc"] = 2.59e-03

    # Constrained parameters
    p["TMaxGLCen"] = 0.041
    p["TMaxGLCce"] = 0.239
    p["TMaxGLCeg"] = 0.147
    p["TMaxGLCcg"] = 0.0016

    p["kHKPFKn"] = 0.0504
    p["kHKPFKg"] = 0.185

    p["kPGKn"] = 3.97
    p["kPGKg"] = 401.7

    p["kPKn"] = 36.7
    p["kPKg"] = 135.2

    p["kCKnms"] = 0.00028
    p["kCKgms"] = 0.00001

    p["PScapVn"] = 1.66
    p["PScapVg"] = 0.87

    p["VMaxMitooutn"] = 0.164
    p["VMaxMitooutg"] = 0.064

    p["TnNADH"] = 10330
    p["TgNADH"] = 150

    p["VMaxMitoinn"] = 0.1303
    p["VMaxMitoing"] = 5.7

    p["vATPasesn"] = 0.1695
    p["vATPasesg"] = 0.1404

    p["kLDHnms"] = 0.72
    p["kLDHgms"] = 0.071

    p["TMaxLACne"] = p["FourCIN"] * 24.3
    p["TMaxLACge"] = 106.1
    p["TMaxLACec"] = 0.25

    p["LACa"] = 0.506

    # Pump and conductance corrections
    p["kPumpn"] = 2.2e-06
    p["gNan"] = 0.0136
    p["gNag"] = 0.0061
    p["gKpas"] = 0.2035
    p["kPumpg"] = 4.5e-07
    p["vPumpg0"] = 0.0687
    
    # -------------------------------------------------
    # Adenine nucleotide / adenylate kinase parameters
    # Original MATLAB: ADP.m and AMPATP.m
    # -------------------------------------------------

    p["qAK"] = 0.92
    p["A"] = 2.212       # mmol/L
    
    # -------------------------------------------------
    # Baseline cerebral blood-flow / vascular parameters
    # Original MATLAB simulation.m
    # -------------------------------------------------

    p["F0"] = 0.012
    p["Vv0"] = 0.021
    return p


# Simple test
if __name__ == "__main__":
    params = get_parameters()
    

    print("Number of parameters loaded:", len(params))
    print("Neuronal volume Vn:", params["Vn"])
    print("Extracellular/neuron volume ratio ren:", params["ren"])
    print("Faraday constant F:", params["F"])
    print("TMaxLACne:", params["TMaxLACne"])
    print("Neuronal LDH forward constant:", params["kLDHnps"])
    print("Neuronal NADH shuttle capacity:", params["TnNADH"])
    print("Stabilized kCKnps:", params["kCKnps"])
    print("Stabilized TMaxLACne:", params["TMaxLACne"])
    print("Stabilized neuronal ATPase:", params["vATPasesn"])
    print("Stabilized neuronal mitochondrial output:", params["VMaxMitooutn"])