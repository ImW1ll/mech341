import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from brayton import *
from sensible_enthalpy import sensible_enthalpy_PT

# Your existing parameters and computed states

m_dot = 65  # kg/s
m_al = 9  # kg/s
q_preheater_per_al = 1e6  # J/kg
q_reheat_per_al = 0
q_boiler_per_al = 12.5e6 - q_reheat_per_al # Requires revision for exactitude
Q_preheater_total = q_preheater_per_al * m_al
q_reheat_total = q_reheat_per_al * m_al
Q_boiler_total = q_boiler_per_al * m_al
 
P_boiler   = 25e6   # Steam generation pressure
P_stage_1_out = 8e6
P_bleed_1  = 10e6    # Pressure at which turbine 1 bleeds steam
P_condenser = 0.01e6 # Expansion continues to condenser pressure         
T_boiler_out = 773.15 # well at least we need to limit it
eff_turbine = 0.85
eff_pump = 0.9
fluid = 'Water'

#--------------------------%--------------------------
#                    Rankine Cycle
#--------------------------%--------------------------

# --- State 1: After condenser (saturated liquid) ---
P1 = P_condenser
Q1 = 0
H1 = PropsSI('H','P',P1,'Q',Q1, fluid)
S1 = PropsSI('S','P',P1,'Q',Q1, fluid)
T1 = PropsSI('T','P',P1,'Q',Q1, fluid)
Sens1 = sensible_enthalpy_PT(P1, T1)

# --- State 2: After compressor 1 ---
P2 = P_bleed_1
H2s = PropsSI('H','P',P2,'S',S1, fluid)
H2 = H1 + (H2s - H1) / eff_pump
S2 = PropsSI('S','P',P2,'H',H2, fluid)
T2 = PropsSI('T','P',P2,'H',H2, fluid)
Sens2 = sensible_enthalpy_PT(P2, T2)

# --- State 3: After OFWH 1 ---
P3 = P_bleed_1
Q3 = 0
H3 = PropsSI('H','P',P3,'Q',Q3,fluid)
S3 = PropsSI('S','P',P3,'Q',Q3,fluid)
T3 = PropsSI('T','P',P3,'Q',Q3,fluid)
Sens3 = sensible_enthalpy_PT(P3, T3)

# --- State 4: After compressor 2 ---
P4 = P_boiler
H4s = PropsSI('H','P',P4,'S',S3,fluid)
H4 = H3 + (H4s - H3) / eff_pump
S4 = PropsSI('S','P',P4,'H',H4,fluid)
T4 = PropsSI('T','P',P4,'H',H4,fluid)
Sens4 = sensible_enthalpy_PT(P4, T4)

# --- State 5: After boiler (to turbine inlet) ---
# Add Q_boiler/m_dot to H3
H5 = H4 + Q_boiler_total / m_dot
P5 = P_boiler
T5 = PropsSI('T','P',P5,'H',H5, fluid)
S5 = PropsSI('S','P',P5,'H',H5, fluid)
Sens5 = sensible_enthalpy_PT(P5, T5)
assert T5 < T_boiler_out, "Temperature out of boiler is too big bozzo"

# --- State 6: After first turbine ---
P6 = P_stage_1_out
H6s = PropsSI('H','S',S5,'P',P6, fluid)
H6 = H5 + (H6s - H5) * eff_turbine
T6 = PropsSI('T','P',P6,'H',H6, fluid)
S6 = PropsSI('S','P',P6,'H',H6, fluid)
Sens6 = sensible_enthalpy_PT(P6, T6)

# --- State 6I: First Bled Water ---
P6I = P_bleed_1
H6Is = PropsSI('H','P',P6I,'S',S5,fluid)
H6I = H5 - (H5 - H6Is) * eff_turbine
S6I = PropsSI('S','P',P6I,'H',H6I,fluid)
T6I = PropsSI('T','P',P6I,'H',H6I,fluid)
Sens6I = sensible_enthalpy_PT(P6I, T6I)

# Bleeding fraction
x_frac = (H3 - H2) / (H6I - H2)

# ---- State 7: After reheat between turbine 1 and turbine 2 -----
H7 = H6 + Q_preheater_total / (m_dot * (1 - x_frac))
P7 = P_bleed_1
S7 = PropsSI('S','P',P7,'H',H7,fluid)
T7 = PropsSI('T','P',P7,'H',H7,fluid)
Sens7 = sensible_enthalpy_PT(P7, T7)

"""
# --- State 10II: Second Bled Water ---
P10II = P_bleed_2
H10IIs = PropsSI('H','P',P10II,'S',S10I,fluid)
H10II = H10I - (H10I - H10IIs) * eff_turbine
S10II = PropsSI('S','P',P10II,'H',H10II,fluid)
T10II = PropsSI('T','P',P10II,'H',H10II,fluid)
Sens10II = sensible_enthalpy_PT(P10II, T10II)

"""

# --- State 8: After second turbine, before condenser ---
P8 = P_condenser
H8s = PropsSI('H','S',S7,'P',P8, fluid)
# Notice the sign in the standard isentropic step is H1s - H6, but you added (H6 - H1s)/eff_turbine 
H8 = H7 - (H7 - H8s) * eff_turbine
T8 = PropsSI('T','P',P8,'H',H8, fluid)
S8 = PropsSI('S','P',P8,'H',H8, fluid)
Sens8 = sensible_enthalpy_PT(P8, T8)

#--------------------------%--------------------------
#              Mass Flow Rate Calculations
#--------------------------%--------------------------

m_frac = 1 - x_frac

# ============= Work math ================
# Work output
W_turb1 = m_dot * ((H5 - H6I) + m_frac * (H6I - H6))
W_turb2 = m_dot * m_frac * (H7 - H8)
W_pump_1 = (H2 - H1) * m_frac * m_dot
W_pump_2 = (H4 - H3) * m_dot

W_net_water = W_turb1 + W_turb2 - W_pump_1  - W_pump_2

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

#--------------------------%--------------------------
#              Hydrogen Related Work
#--------------------------%--------------------------

# --- SATP Conditions ---

gas = 'Hydrogen'
P_SATP = 0.1013e6 # Pa
T_SATP = 298.19 # K

# --- Out of the reactor the H2 is at 300bar and 500C ---
M_Al = 26.98
M_H2 = 2.016
LHV_H2 = 120e6   # J/kg
eff_fc = 0.5     # Fuel-cell efficiency

kg_H2_per_kg_Al = (3*M_H2)/(2*M_Al)  # ~0.1117
m_H2 = m_al * kg_H2_per_kg_Al       # mass flow of H2 (kg/s)

PH2 = 3e7
TH2 = 500 + 273.15 # K
eff_h2 = 0.85
PH2_out = 300000 # 1 bar, yeah I'm that optimistic

HH2 = PropsSI('H', 'P', PH2, 'T', TH2, gas)
SH2 = PropsSI('S', 'P', PH2, 'T', TH2, gas)


h_s_out = PropsSI('H','P',PH2_out,'S', SH2, gas)
h_out   = HH2 - eff_fc*(HH2 - h_s_out)

H2T_out = PropsSI('T', 'H', h_out, 'P', PH2_out, gas)

w_TurbH2 = (HH2 - h_out) * m_H2

# -- Brayton-cycle parameters (replace with real design if needed) --
T_in_C = 200.58        # Inlet H2 temperature (°C) - or use your known T_in
T_in = T_in_C + 273.15 # K
P_in_bar = 1.0
P_in = P_in_bar*1e5    # 1 bar -> Pa

# For the simplest case, assume you have enough air to keep T3 within reason:
m_air = 150            # [kg/s], adjust to control turbine inlet temp
pr = 10.0              # overall pressure ratio
cp_air = 1005.0        # J/(kg*K) Idealized at the moment, check error of idealization
gamma_air = 1.4

eta_comp = 0.88        # compressor isentropic efficiency
eta_turb = 0.90        # turbine isentropic efficiency

# Yeah I know I made it simpler
brayton_states_and_work = brayton_cycle_h2(m_H2,T_boiler_out,PH2,2)
W_br = brayton_states_and_work['Work']

# --- Supplementary Work ---

TB1 = T_boiler_out
PB1 = 30e6
SB1 = PropsSI('S','T',TB1,'P',PB1,gas)
HB1 = PropsSI('H','T',TB1,'P',PB1,gas)
""" I've commented this out because I wanna push asap and don't me have time to fix
PB2 = 3e6 # As instructed by the professor
TB2 = T2_br
HB2 = PropsSI('H','T',TB2,'P',PB2,gas)
SB2 = PropsSI('S','T',TB2,'P',PB2,gas)

PB3 = PB2
TB3 = T3_br
HB3 = PropsSI('H','T',TB3,'P',PB3,gas)
SB3 = PropsSI('S','T',TB3,'P',PB3,gas)

PB4 = P_in
TB4 = T4_br
HB4 = PropsSI('H','T',TB4,'P',PB4,gas)
SB4 = PropsSI('S','T',TB4,'P',PB4,gas)

these states are also in jail
{
        'State': 'B2',
        'P (MPa)': round(PB2/1e6,2),
        'T (C)': round(T2_br - 273.19,2),
        'H (kJ/kg)': round(HB2/1e3,2),
        'S (kJ/kg.K)': round(SB2/1e3,2),
    },
    {
        'State': 'B3',
        'P (MPa)': round(PB3/1e6,2),
        'T (C)': round(T3_br - 273.19,2),
        'H (kJ/kg)': round(HB3/1e3,2),
        'S (kJ/kg.K)': round(SB3/1e3,2),
    },
    {
        'State': 'B4',
        'P (MPa)': round(PB4/1e6/1e6,2),
        'T (C)': round(T4_br - 273.19,2),
        'H (kJ/kg)': round(HB4/1e3,2),
        'S (kJ/kg.K)': round(SB4/1e3,2) 
    },

"""

# ================= Output results =======================

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
        'Sens. H': Sens1 / 1000,
    },
    {
        'State': '2',
        'P (MPa)': round(P2/1e6, 2),
        'T (C)': round(T2 - 273.15, 2),
        'H (kJ/kg)': round(H2/1e3, 2),
        'S (kJ/kg.K)': round(S2/1e3, 2),
        'Sens. H': Sens2 / 1000,
    },
    {
        'State': '3',
        'P (MPa)': round(P3/1e6, 2),
        'T (C)': round(T3 - 273.15, 2),
        'H (kJ/kg)': round(H3/1e3, 2),
        'S (kJ/kg.K)': round(S3/1e3, 2),
        'Sens. H': Sens3 / 1000,
    },
    {
        'State': '4',
        'P (MPa)': round(P4/1e6, 2),
        'T (C)': round(T4 - 273.15, 2),
        'H (kJ/kg)': round(H4/1e3, 2),
        'S (kJ/kg.K)': round(S4/1e3, 2),
        'Sens. H': Sens4 / 1000,
    },
    {
        'State': '5',
        'P (MPa)': round(P5/1e6, 2),
        'T (C)': round(T5 - 273.15, 2),
        'H (kJ/kg)': round(H5/1e3, 2),
        'S (kJ/kg.K)': round(S5/1e3, 2),
        'Sens. H': Sens5 / 1000,
    },
    {
        'State': '6',
        'P (MPa)': round(P6/1e6, 2),
        'T (C)': round(T6 - 273.15, 2),
        'H (kJ/kg)': round(H6/1e3, 2),
        'S (kJ/kg.K)': round(S6/1e3, 2),
        'Sens. H': Sens6 / 1000,
    },
    {
        'State': '6I',
        'P (MPa)': round(P6I/1e6, 2),
        'T (C)': round(T6I - 273.15, 2),
        'H (kJ/kg)': round(H6I/1e3, 2),
        'S (kJ/kg.K)': round(S6I/1e3, 2),
        'Sens. H': Sens3 / 1000,
    },
    {
        'State': '7',
        'P (MPa)': round(P7/1e6, 2),
        'T (C)': round(T7 - 273.15, 2),
        'H (kJ/kg)': round(H7/1e3, 2),
        'S (kJ/kg.K)': round(S7/1e3, 2),
        'Sens. H': Sens7 / 1000,
    },
    {
        'State': '8',
        'P (MPa)': round(P8/1e6, 2),
        'T (C)': round(T8 - 273.15, 2),
        'H (kJ/kg)': round(H8/1e3, 2),
        'S (kJ/kg.K)': round(S8/1e3, 2),
        'Sens. H': Sens8 / 1000,
    },
    
]

print(f"New brayton shaft power {W_br:.2f}")
ACC_power_output = W_turb1 + W_turb2 + W_br

work_kWh_per_tonne = (W_net_water + W_br + w_TurbH2) * 1000 / (3.6e6 *m_al)
print(f"Work per tonne of Al: {work_kWh_per_tonne:.2f} kWh/tonne")
print(f"Thermal efficiency of Steam Rankine Cycle: {eta_th_water*100:.2f} %")
print(f"Bleeding fraction after first Turbine: {(x_frac*100):.2f} %")
print(f"Thermal efficiency of Total Cycle: {eta_total*100:.2f} %")
print(f"Actual Power Output {(ACC_power_output/1e6):.2f} MW")

df_states = pd.DataFrame(states_data)
print("\n=== Thermodynamic States ===")
print(df_states)

# If you want to display as a Matplotlib table:
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