import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from brayton_cycle import *

# Your existing parameters and computed states

m_dot = 65  # kg/s
m_al = 8.75  # kg/s
q_preheater_per_al = 1e6  # J/kg
q_reheat_per_al = 0
q_boiler_per_al = 12.5e6 - q_reheat_per_al
Q_preheater_total = q_preheater_per_al * m_al
q_reheat_total = q_reheat_per_al * m_al
Q_boiler_total = q_boiler_per_al * m_al
 
P_condenser = 0.01e6     
P_stage_1 = 22.2e6          
P_stage_2 = 12e6
P_bleed_1 = 10e6
P_bleed_2 = 4e6           
T_boiler_out = 773.15
eff_turbine = 0.85
eff_pump = 0.9
fluid = 'Water'

# --- State 1: After condenser (saturated liquid) ---
P1 = P_condenser
Q1 = 0
H1 = PropsSI('H','P',P1,'Q',Q1, fluid)
S1 = PropsSI('S','P',P1,'Q',Q1, fluid)
T1 = PropsSI('T','P',P1,'Q',Q1, fluid)

# --- State 2: After compressor 1 ---
P2 = P_bleed_2
H2s = PropsSI('H','P',P2,'S',S1, fluid)
H2 = H1 + (H2s - H1) / eff_pump
S2 = PropsSI('S','P',P2,'H',H2, fluid)
T2 = PropsSI('T','P',P2,'H',H2, fluid)

# --- State 3: After OFWH 1 ---
P3 = P_bleed_2
Q3 = 0
H3 = PropsSI('H','P',P3,'Q',Q3,fluid)
S3 = PropsSI('S','P',P3,'Q',Q3,fluid)
T3 = PropsSI('T','P',P3,'Q',Q3,fluid)

# --- State 4: After compressor 2 ---
P4 = P_bleed_1
H4s = PropsSI('H','P',P4,'S',S3,fluid)
H4 = H3 + (H4s - H3) / eff_pump
S4 = PropsSI('S','P',P4,'H',H4,fluid)
T4 = PropsSI('T','P',P4,'H',H4,fluid)

# --- State 5: OFWH 2 ---
P5 = P_bleed_1
Q5 = 0
H5 = PropsSI('H','P',P5,'Q',Q5,fluid)
S5 = PropsSI('S','P',P5,'Q',Q5,fluid)
T5 = PropsSI('T','P',P5,'Q',Q5,fluid)

# --- State 6: After compressor 3 ---
P6 = P_stage_1
H6s = PropsSI('H','P',P6,'S',S5,fluid)
H6 = H5 + (H6s - H5) / eff_pump
S6 = PropsSI('S','P',P6,'H',H6,fluid)
T6 = PropsSI('T','P',P6,'H',H6,fluid)

# --- State 7: After boiler (to turbine inlet) ---
# Add Q_boiler/m_dot to H3
H7 = H6 + Q_boiler_total / m_dot
P7 = P_stage_1
T7 = PropsSI('T','P',P7,'H',H7, fluid)
S7 = PropsSI('S','P',P7,'H',H7, fluid)

# --- State 8: After first turbine ---
P8 = P_stage_2
H8s = PropsSI('H','S',S7,'P',P8, fluid)
H8 = H7 + (H8s - H7) * eff_turbine
T8 = PropsSI('T','P',P8,'H',H8, fluid)
S8 = PropsSI('S','P',P8,'H',H8, fluid)

# ---- State 9: After reheat between turbine 1 and turbine 2 -----
H9 = H8 + Q_preheater_total / m_dot
P9 = P_stage_2
S9 = PropsSI('S','P',P9,'H',H9,fluid)
T9 = PropsSI('T','P',P9,'H',H9,fluid)

# --- State 10I: First Bled Water ---
P10I = P_bleed_1
H10Is = PropsSI('H','P',P10I,'S',S9,fluid)
H10I = H9 - (H9 - H10Is) * eff_turbine
S10I = PropsSI('S','P',P10I,'H',H10I,fluid)
T10I = PropsSI('T','P',P10I,'H',H10I,fluid)

# --- State 10II: Second Bled Water ---
P10II = P_bleed_2
H10IIs = PropsSI('H','P',P10II,'S',S10I,fluid)
H10II = H10I - (H10I - H10IIs) * eff_turbine
S10II = PropsSI('S','P',P10II,'H',H10II,fluid)
T10II = PropsSI('T','P',P10II,'H',H10II,fluid)

# --- State 10: After second turbine, before condenser ---
P10 = P_condenser
H10s = PropsSI('H','S',S10II,'P',P10, fluid)
# Notice the sign in the standard isentropic step is H1s - H6, but you added (H6 - H1s)/eff_turbine 
H10 = H10II - (H10II - H10s) * eff_turbine
T10 = PropsSI('T','P',P10,'H',H10, fluid)
S10 = PropsSI('S','P',P10,'H',H10, fluid)

#--------------------------%--------------------------
#              Mass Flow Rate Calculations
#--------------------------%--------------------------

y_frac = (H5 - H4) / (H10I - H4)
x_frac = (1 - y_frac) * (H2 + H3) / (H10II - H2)
m_frac = 1-x_frac-y_frac
assert(m_frac + y_frac + x_frac == 1,"Issue with Mass Fraction: maybe check enthalpy compatibility?")


# Collect all states in a dict
# Convert P -> MPa (divide by 1e6), H -> kJ/kg (divide by 1e3), S -> kJ/(kg.K) (divide by 1e3)
# Round to 2 decimals
states_data = [
    {
        'State': '1',
        'P (MPa)': round(P1/1e6, 2),
        'T (C)': round(T1 - 273.15, 2),
        'H (kJ/kg)': round(H1/1e3, 2),
        'S (kJ/kg.K)': round(S1/1e3, 2),
    },
    {
        'State': '2',
        'P (MPa)': round(P2/1e6, 2),
        'T (C)': round(T2 - 273.15, 2),
        'H (kJ/kg)': round(H2/1e3, 2),
        'S (kJ/kg.K)': round(S2/1e3, 2),
    },
    {
        'State': '3',
        'P (MPa)': round(P3/1e6, 2),
        'T (C)': round(T3 - 273.15, 2),
        'H (kJ/kg)': round(H3/1e3, 2),
        'S (kJ/kg.K)': round(S3/1e3, 2),
    },
    {
        'State': '4',
        'P (MPa)': round(P4/1e6, 2),
        'T (C)': round(T4 - 273.15, 2),
        'H (kJ/kg)': round(H4/1e3, 2),
        'S (kJ/kg.K)': round(S4/1e3, 2),
    },
    {
        'State': '5',
        'P (MPa)': round(P5/1e6, 2),
        'T (C)': round(T5 - 273.15, 2),
        'H (kJ/kg)': round(H5/1e3, 2),
        'S (kJ/kg.K)': round(S5/1e3, 2),
    },
    {
        'State': '6',
        'P (MPa)': round(P6/1e6, 2),
        'T (C)': round(T6 - 273.15, 2),
        'H (kJ/kg)': round(H6/1e3, 2),
        'S (kJ/kg.K)': round(S6/1e3, 2),
    },
    {
        'State': '7',
        'P (MPa)': round(P7/1e6, 2),
        'T (C)': round(T7 - 273.15, 2),
        'H (kJ/kg)': round(H7/1e3, 2),
        'S (kJ/kg.K)': round(S7/1e3, 2),
    },
    {
        'State': '8',
        'P (MPa)': round(P8/1e6, 2),
        'T (C)': round(T8 - 273.15, 2),
        'H (kJ/kg)': round(H8/1e3, 2),
        'S (kJ/kg.K)': round(S8/1e3, 2),
    },
    {
        'State': '9',
        'P (MPa)': round(P9/1e6, 2),
        'T (C)': round(T9 - 273.15, 2),
        'H (kJ/kg)': round(H9/1e3, 2),
        'S (kJ/kg.K)': round(S9/1e3, 2),
    },
    {
        'State': '10',
        'P (MPa)': round(P10/1e6, 2),
        'T (C)': round(T10 - 273.15, 2),
        'H (kJ/kg)': round(H10/1e3, 2),
        'S (kJ/kg.K)': round(S10/1e3, 2),
    },
    {
        'State': '10I',
        'P (MPa)': round(P10I/1e6, 2),
        'T (C)': round(T10I - 273.15, 2),
        'H (kJ/kg)': round(H10I/1e3, 2),
        'S (kJ/kg.K)': round(S10I/1e3, 2),
    },
    {
        'State': '10II',
        'P (MPa)': round(P10II/1e6, 2),
        'T (C)': round(T10II - 273.15, 2),
        'H (kJ/kg)': round(H10II/1e3, 2),
        'S (kJ/kg.K)': round(S10II/1e3, 2),
    }
]

# ============= Work math ================
# Work output
W_turb1 = (H7 - H8) * m_dot
W_turb2 = m_dot * ((H9 - H10I) + (H10I - H10II)*(1-y_frac) + (H10II - H10)*(1-x_frac-y_frac))

W_pump_1 = (H2 - H1) * (1-x_frac-y_frac) * m_dot
W_pump_2 = (H4 - H3) * (1-y_frac) * m_dot
W_pump_3 = (H6 - H5) * m_dot

W_net_water = W_turb1 + W_turb2 - W_pump_1 - W_pump_2 - W_pump_3


# ============ Heat math ===============
# Heat input (should match given totals)
Q_pre = Q_preheater_total
Q_boiler = Q_boiler_total
Q_in = Q_pre + Q_boiler

eta_th_water = W_net_water / Q_in
eta_total = (W_net_water) / Q_in

# ============ Cooling math ===============
orc_fluid = 'Isobutane'
P_b = 1.6e6
P_a = 0.1013e6
T_river = 288.19
P_river = 0.1013e6
m_dot_orc = 32
m_dot_river = 1.5*m_al
H_river = PropsSI('H','T',T_river,'P',P_river,fluid)
H_al_water = H_river / m_dot_river
T_al_water = PropsSI('T','H',H_al_water,'P',P_river,fluid)

# ================= H2 Theoretical Work ================
# we need ceramic coating and really good coating
# Out of the reactor the H2 is at 300bar and 500C
M_Al = 26.98
M_H2 = 2.016
LHV_H2 = 120e6   # J/kg
eff_fc = 0.5     # Fuel-cell efficiency

kg_H2_per_kg_Al = (3*M_H2)/(2*M_Al)  # ~0.1117
m_H2 = m_al * kg_H2_per_kg_Al       # mass flow of H2 (kg/s)

PH2 = 3e7
TH2 = 500 + 273.15 # K
eff_h2 = 0.85
PH2_out = 100000 # 1 bar, yeah I'm that optimistic

HH2 = PropsSI('H', 'P', PH2, 'T', TH2, 'Hydrogen')
SH2 = PropsSI('S', 'P', PH2, 'T', TH2, 'Hydrogen')


h_s_out = PropsSI('H','P',PH2_out,'S', SH2, 'Hydrogen')
h_out   = HH2 - eff_fc*(HH2 - h_s_out)

H2T_out = PropsSI('T', 'H', h_out, 'P', PH2_out, 'Hydrogen')

w_TurbH2 = (HH2 - h_out) * m_H2

# -- Brayton-cycle parameters (replace with real design if needed) --
T_in_C = 200.58         # Inlet H2 temperature (°C) - or use your known T_in
T_in   = T_in_C + 273.15 # K
P_in_bar = 1.0
P_in = P_in_bar*1e5     # 1 bar -> Pa

# For the simplest case, assume you have enough air to keep T3 within reason:
m_air = 80             # [kg/s], adjust to control turbine inlet temp
pr = 10.0               # overall pressure ratio
cp_air = 1005.0         # J/(kg*K)
gamma_air = 1.4

eta_comp = 0.88         # compressor isentropic efficiency
eta_turb = 0.90         # turbine isentropic efficiency

# -- Get the net Brayton output power [W] for the hydrogen cycle --
W_br, T2_br, T3_br, T4_br = brayton_cycle_h2(m_H2, T_in, P_in,
                                            pr,
                                            gamma_air, 
                                            LHV_H2, m_air,
                                            eta_comp=eta_comp,
                                            eta_turb=eta_turb)

# Expansion of H2 through a gas turbine
'''HH1 = PropsSI('H','T',T_boiler_out,'P',30e6,'hydrogen')
SH1 = PropsSI('H','T',T_boiler_out,'P',30e6,'hydrogen')
HH2s = PropsSI('H','S',SH1,'P',0.1013e6,'hydrogen')
HH2 = HH1 + (HH2s - HH1) * eff_turbine
W_dot_H2_brayton = (HH1 - HH2) * m_H2'''

print(f"New brayton shaft power {W_br:.2f}")
ACC_power_output = W_turb1 + W_turb2 + w_TurbH2 + W_br

work_kWh_per_tonne = (W_net_water + W_br + w_TurbH2) * 1000 / (3.6e6 *m_al)
print(f"Work per tonne of Al: {work_kWh_per_tonne:.2f} kWh/tonne")
print(f"Thermal efficiency of Steam Rankine Cycle: {eta_th_water*100:.2f} %")
#print(f"Thermal efficiency of Organic Rankine Cycle: {eta_th_orc*100:.2f} %")
print(f"Thermal efficiency of Total Cycle: {eta_total*100:.2f} %")
print(f"Actual Power Output {(ACC_power_output/1e6):.2f} MW")

df_states = pd.DataFrame(states_data)
print("\n=== Thermodynamic States ===")
print(df_states)

# If you want to display as a Matplotlib table (like an image):
fig, ax = plt.subplots()
ax.axis('tight')
ax.axis('off')
table = ax.table(
    cellText=df_states.values,
    colLabels=df_states.columns,
    loc='center'
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(df_states.columns))))
fig.tight_layout()
plt.show()
