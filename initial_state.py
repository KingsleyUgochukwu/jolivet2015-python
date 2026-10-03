import numpy as np

# Initial steady-state values from:
# Jolivet et al. (2015) MATLAB model
# States/initstate.data
#
# IMPORTANT:
# The order of these 33 variables must not be changed because
# the ODE system refers to each state by its position.

initial_state = np.array([
    7.9717478717769090,       # 1  Na_n
    15.105962705116731,       # 2  Na_g
    1.1952401106798076,       # 3  GLC_n
    1.1916831724306007,       # 4  GLC_g
    0.0046031315201704066,    # 5  GAP_n
    0.004592153756399458,     # 6  GAP_g
    0.016409369231409576,     # 7  PEP_n
    0.014203935093932973,     # 8  PEP_g
    0.17389794176222959,      # 9  PYR_n
    0.16297929416942439,      # 10 PYR_g
    0.5981237834039330,       # 11 LAC_n
    0.6001383440530418,       # 12 LAC_g
    0.006243658877436395,     # 13 NADH_n
    0.10386865243459381,      # 14 NADH_g
    2.1973762129517347,       # 15 ATP_n
    2.1951002997965290,       # 16 ATP_g
    4.9460223423181864,       # 17 PCr_n
    4.9242059083825644,       # 18 PCr_g
    0.027584951432999300,     # 19 O2_n
    0.027594056537295686,     # 20 O2_g
    6.9808721131857396,       # 21 O2_c
    4.500386628495174,        # 22 GLC_c
    0.5488507807424990,       # 23 LAC_c
    0.021,                    # 24 V_v
    0.057503371246197574,     # 25 dHb
    2.479690857581554,        # 26 GLC_e
    0.5991310930946384,       # 27 LAC_e
    -73.59279632137904,       # 28 V_m
    0.9937194195642355,       # 29 h
    0.018509072079738374,     # 30 n
    0.00005100694314844131,   # 31 Ca_n
    0.12345447794054716,      # 32 NADH_n_mito
    0.12483363148236443       # 33 NADH_g_mito
])

# State names in exactly the same order
state_names = [
    "Na_n", "Na_g",
    "GLC_n", "GLC_g",
    "GAP_n", "GAP_g",
    "PEP_n", "PEP_g",
    "PYR_n", "PYR_g",
    "LAC_n", "LAC_g",
    "NADH_n", "NADH_g",
    "ATP_n", "ATP_g",
    "PCr_n", "PCr_g",
    "O2_n", "O2_g", "O2_c",
    "GLC_c", "LAC_c",
    "V_v", "dHb",
    "GLC_e", "LAC_e",
    "V_m", "h", "n", "Ca_n",
    "NADH_n_mito", "NADH_g_mito"
]

# Safety check
assert len(initial_state) == 33
assert len(state_names) == 33

print("Number of state variables:", len(initial_state))
print("Initial neuronal ATP:", initial_state[14])
print("Initial membrane potential:", initial_state[27])