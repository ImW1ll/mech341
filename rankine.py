import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from brayton_cycle import *

# =============== GIVEN PARAMETERS ======================
m_dot = 8.5*10      # kg/s
m_al = 2*10           # kg/s
q_preheater_per_al = 1.5e6  # J/kg
q_reheat_per_al    = 0
q_boiler_per_al    = 12.5e6 - q_reheat_per_al
Q_preheater_total  = q_preheater_per_al * m_al
q_reheat_total     = q_reheat_per_al * m_al
Q_boiler_total     = q_boiler_per_al * m_al

P_condenser = 0.01e6
P_boiler    = 23e6
P_t1out     = 5e6
T_boiler_out = 773.15   # 500 °C
eff_turbine = 0.85
eff_pump    = 0.85
fluid       = 'Water'

# =============== HELPER FUNCTIONS ======================
def pump_enthalpy(h_in, s_in, P_out, eta):
    h_s_out = PropsSI('H', 'P', P_out, 'S', s_in, fluid)
    h_out   = h_in + (h_s_out - h_in)/eta
    return h_out

def turbine_enthalpy(h_in, s_in, P_out, eta):
    h_s_out = PropsSI('H','P',P_out,'S', s_in, fluid)
    h_out   = h_in - eta*(h_in - h_s_out)
    return h_out

# --- State 2: After condenser (saturated liquid) ---
P2 = P_condenser
H2 = PropsSI('H','P',P2,'Q',0, fluid)
S2 = PropsSI('S','P',P2,'Q',0, fluid)
T2 = PropsSI('T','P',P2,'Q',0, fluid)

# --- State 3: After pump from P_condenser -> P_t1out ---
P3 = P_t1out
H3 = pump_enthalpy(H2, S2, P3, eff_pump)
S3 = PropsSI('S','P',P3,'H',H3, fluid)
T3 = PropsSI('T','P',P3,'H',H3, fluid)

# --- State 4: After boiler (to turbine inlet) ---
# Add Q_boiler_total/m_dot to H3
H4 = H3 + Q_boiler_total / m_dot
P4 = P_boiler
# We find T4, S4 from P4, H4
T4 = PropsSI('T','P',P4,'H',H4, fluid)
S4 = PropsSI('S','P',P4,'H',H4, fluid)

# --- State 5: After first turbine (expand from P4 -> P_t1out) ---
H5 = turbine_enthalpy(H4, S4, P_t1out, eff_turbine)
P5 = P_t1out
T5 = PropsSI('T','P',P5,'H',H5, fluid)
S5 = PropsSI('S','P',P5,'H',H5, fluid)

# --- State 6: After "reheat" between turbine 1 and turbine 2 -----
# Adds Q_preheater_total / m_dot
H6 = H5 + Q_preheater_total / m_dot
P6 = P5
T6 = PropsSI('T','P',P6,'H',H6, fluid)
S6 = PropsSI('S','P',P6,'H',H6, fluid)

# --- State 1: After second turbine, before condenser ---
# Expand from P6 -> P_condenser
P1 = P_condenser
H1 = turbine_enthalpy(H6, S6, P1, eff_turbine) 
T1 = PropsSI('T','P',P1,'H',H1, fluid)
S1 = PropsSI('S','P',P1,'H',H1, fluid)

# ============= CREATE DICTIONARY OF STATES ============
states_data = [
    {
        'State': '1',
        'P (MPa)': round(P1/1e6, 2),
        'T (°C)':  round(T1 - 273.15, 2),
        'H (kJ/kg)':round(H1/1e3, 2),
        'S (kJ/kg.K)':round(S1/1e3, 2),
    },
    {
        'State': '2',
        'P (MPa)': round(P2/1e6, 2),
        'T (°C)':  round(T2 - 273.15, 2),
        'H (kJ/kg)':round(H2/1e3, 2),
        'S (kJ/kg.K)':round(S2/1e3, 2),
    },
    {
        'State': '3',
        'P (MPa)': round(P3/1e6, 2),
        'T (°C)':  round(T3 - 273.15, 2),
        'H (kJ/kg)':round(H3/1e3, 2),
        'S (kJ/kg.K)':round(S3/1e3, 2),
    },
    {
        'State': '4',
        'P (MPa)': round(P4/1e6, 2),
        'T (°C)':  round(T4 - 273.15, 2),
        'H (kJ/kg)':round(H4/1e3, 2),
        'S (kJ/kg.K)':round(S4/1e3, 2),
    },
    {
        'State': '5',
        'P (MPa)': round(P5/1e6, 2),
        'T (°C)':  round(T5 - 273.15, 2),
        'H (kJ/kg)':round(H5/1e3, 2),
        'S (kJ/kg.K)':round(S5/1e3, 2),
    },
    {
        'State': '6',
        'P (MPa)': round(P6/1e6, 2),
        'T (°C)':  round(T6 - 273.15, 2),
        'H (kJ/kg)':round(H6/1e3, 2),
        'S (kJ/kg.K)':round(S6/1e3, 2),
    }
]

# ============== H2 turbine ====================
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


# ============= Work and Efficiency Math ================
W_turb1 = (H4 - H5) * m_dot
W_turb2 = (H6 - H1) * m_dot
W_pump  = (H3 - H2) * m_dot

W_net = W_turb1 + W_turb2 - W_pump + w_TurbH2

Q_pre = Q_preheater_total
Q_boiler = Q_boiler_total
Q_in = Q_pre + Q_boiler


# ============== H2 TURBINE (BRAYTON CYCLE) ====================
# From stoichiometry or reaction, the H2 mass flow per kg Al:
M_Al = 26.98
M_H2 = 2.016
kg_H2_per_kg_Al = (3*M_H2)/(2*M_Al)  # ~0.1117

m_H2 = m_al * kg_H2_per_kg_Al       # [kg/s] mass flow of hydrogen
LHV_H2 = 120e6                      # J/kg (hydrogen lower heating value)

# -- Brayton-cycle parameters (replace with real design if needed) --
T_in_C = 200.58         # Inlet H2 temperature (°C) - or use your known T_in
T_in   = T_in_C + 273.15 # K
P_in_bar = 1.0
P_in = P_in_bar*1e5     # 1 bar -> Pa

# For the simplest case, assume you have enough air to keep T3 within reason:
m_air = 5.0             # [kg/s], adjust to control turbine inlet temp
pr = 10.0               # overall pressure ratio
cp_air = 1005.0         # J/(kg*K)
gamma_air = 1.4

eta_comp = 0.88         # compressor isentropic efficiency
eta_turb = 0.90         # turbine isentropic efficiency

# -- Get the net Brayton output power [W] for the hydrogen cycle --
W_br, T2_br, T3_br, T4_br = brayton_cycle_h2(m_H2, T_in, P_in,
                                            pr,
                                            cp_air, gamma_air, 
                                            LHV_H2, m_air,
                                            eta_comp=eta_comp,
                                            eta_turb=eta_turb)

# ============= Work and Efficiency Math (Rankine + H2) ================
W_turb1 = (H4 - H5) * m_dot
W_turb2 = (H6 - H1) * m_dot
W_pump  = (H3 - H2) * m_dot

# Combine Rankine net with Brayton net
#   (Rankine net = W_turb1 + W_turb2 - W_pump)
W_net_rankine = W_turb1 + W_turb2 - W_pump

# Add the hydrogen Brayton cycle
W_net_total = W_net_rankine + W_br

Q_pre = Q_preheater_total
Q_boiler = Q_boiler_total
Q_in = Q_pre + Q_boiler

eta_th = W_net_rankine / Q_in  # if you only want the Rankine's "thermal" efficiency

# For "total" output from both cycles per mass of Al, you might do:
work_kWh_per_tonne = (W_net_total)*1000 / (3.6e6*m_al)


# =========== Print Results ===========
print(f"T out after H2 expansion = {(H2T_out - 273.15):.2f} C")
print(f"Work turbine 1 = {W_turb1:.2f}  W")
print(f"Work turbine 2 = {W_turb2:.2f}  W")
print(f"Pump work      = {W_pump:.2f}  W")
print(f"Net Rankine Work = {W_net:.2f}  W")
print(f"H2 Turbine Work = {w_TurbH2:.2} W")
print(f"Rankine Thermal Efficiency = {eta_th*100:.2f} %")

print(f"Total (Rankine) Work per tonne Al = {work_kWh_per_tonne:.2f} kWh/tonne")

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
