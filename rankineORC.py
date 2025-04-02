import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Your existing parameters and computed states

m_dot = 8.5  # kg/s
m_al = 2  # kg/s
q_preheater_per_al = 1e6  # J/kg
q_reheat_per_al = 0
q_boiler_per_al = 12.5e6 - q_reheat_per_al
Q_preheater_total = q_preheater_per_al * m_al
q_reheat_total = q_reheat_per_al * m_al
Q_boiler_total = q_boiler_per_al * m_al
 
P_condenser = 0.1013e6     
P_boiler = 5e6          
P_t1out = 2.75e6           
T_boiler_out = 773.15   
eff_turbine = 0.8
eff_pump = 0.8
fluid = 'Water'

# --- State 2: After condenser (saturated liquid) ---
P2 = P_condenser
Q2 = 0
H2 = PropsSI('H','P',P2,'Q',Q2, fluid)
S2 = PropsSI('S','P',P2,'Q',Q2, fluid)
T2 = PropsSI('T','P',P2,'Q',Q2, fluid)

# --- State 3: After compressor 1 ---
P3 = P_t1out
H3s = PropsSI('H','P',P3,'S',S2, fluid)
H3 = H2 + (H3s - H2) / eff_pump
S3 = PropsSI('S','P',P3,'H',H3, fluid)
T3 = PropsSI('T','P',P3,'H',H3, fluid)

# --- State 4: After boiler (to turbine inlet) ---
# Add Q_boiler/m_dot to H3
H4 = H3 + Q_boiler_total / m_dot
P4 = P_boiler
T4 = PropsSI('T','P',P4,'H',H4, fluid)
S4 = PropsSI('S','P',P4,'H',H4, fluid)

# --- State 5: After first turbine ---
H5s = PropsSI('H','S',S4,'P',P_t1out, fluid)
H5 = H4 + (H5s - H4) * eff_turbine
P5 = P_t1out
T5 = PropsSI('T','P',P5,'H',H5, fluid)
S5 = PropsSI('S','P',P5,'H',H5, fluid)

# ---- State 6: After reheat between turbine 1 and turbine 2 -----
H6 = H5 + Q_preheater_total / m_dot
P6 = P5
S6 = PropsSI('S','P',P6,'H',H6,fluid)
T6 = PropsSI('T','P',P6,'H',H6,fluid)

# --- State 1: After second turbine, before condenser ---
P1 = P_condenser
H1s = PropsSI('H','S',S6,'P',P1, fluid)
# Notice the sign in the standard isentropic step is H1s - H6, but you added (H6 - H1s)/eff_turbine 
H1 = H6 - (H6 - H1s)*eff_turbine
T1 = PropsSI('T','P',P1,'H',H1, fluid)
S1 = PropsSI('S','P',P1,'H',H1, fluid)

#--------------------------%--------------------------
#               Organic Rankine Cycle
#--------------------------%--------------------------

orc_fluid = 'Isobutane'
P_b = 1.6e6
P_a = 0.1013e6
T_river = 288.19
P_river = 0.1013e6
m_dot_orc = 32
m_dot_river = 1.5*m_al

# --- State a: After condensation ---
#H_water_in = PropsSI('H','T',T_river,'P',P_river,fluid)
#H_water_out = PropsSI('H','T',298.19,'P',P_river,fluid)
#Q_bc = m_dot_river*(H_water_out - H_water_in)
#Hc = Hb - Q_bc/m_dot_orc
Qa = 0
Pa = P_a
Ha = PropsSI('H','Q',Qa,'P',Pa,orc_fluid)
Sa = PropsSI('S','Q',Qa,'P',Pa,orc_fluid)
Ta = PropsSI('T','Q',Qa,'P',Pa,orc_fluid)

# --- State b: After pump ---
Pb = P_b
Hbs = PropsSI('H','P',Pb,'S',Sa,orc_fluid)
Hb = Ha + (Hbs - Ha) / eff_pump
Tb = PropsSI('T','H',Hb,'P',Pb,orc_fluid)
Sb = PropsSI('S','H',Hb,'P',Pb,orc_fluid)

# --- State c: After the water condenser, which is the isobutane boiler ---
#q_orc_boiler = (H1 - H2) * m_dot
Tc = T1
Pc = P_b
Hc = PropsSI('H','T',Tc,'P',Pc,orc_fluid)
Sc = PropsSI('S','T',Tc,'P',Pc,orc_fluid)

# --- State d: After the first expansion cycle ---
Pd = P_a
Hds = PropsSI('H','S',Sc,'P',Pd,orc_fluid)
Hd = Hc - (Hc - Hds)*eff_turbine
Td = PropsSI('T','H',Hd,'P',Pd,orc_fluid)
Sd = PropsSI('S','H',Hd,'P',Pd,orc_fluid)

#Currently, we have imaginary cooling. This can be resolved simply by adding a cooling tower of sorts.

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
        'State': 'a',
        'P (MPa)': round(Pa/1e6, 2),
        'T (C)': round(Ta - 273.15, 2),
        'H (kJ/kg)': round(Ha/1e3, 2),
        'S (kJ/kg.K)': round(Sa/1e3, 2),
    },
    {
        'State': 'b',
        'P (MPa)': round(Pb/1e6, 2),
        'T (C)': round(Tb - 273.15, 2),
        'H (kJ/kg)': round(Hb/1e3, 2),
        'S (kJ/kg.K)': round(Sb/1e3, 2),
    },
    {
        'State': 'c',
        'P (MPa)': round(Pc/1e6, 2),
        'T (C)': round(Tc - 273.15, 2),
        'H (kJ/kg)': round(Hc/1e3, 2),
        'S (kJ/kg.K)': round(Sc/1e3, 2),
    },
    {
        'State': 'd',
        'P (MPa)': round(Pd/1e6, 2),
        'T (C)': round(Td - 273.15, 2),
        'H (kJ/kg)': round(Hd/1e3, 2),
        'S (kJ/kg.K)': round(Sd/1e3, 2),
    },
]

# ============= Work math ================
# Work output
W_turb1 = (H4 - H5) * m_dot
W_turb2 = (H6 - H1) * m_dot
W_turb_orc = (Hc - Hd) * m_dot_orc

W_pump = (H3 - H2) * m_dot
W_pump_orc = (Hb - Ha) * m_dot_orc

W_net_water = W_turb1 + W_turb2 - W_pump
W_net_orc = W_turb_orc - W_pump_orc
ACC_power_output = W_turb1 + W_turb2 + W_turb_orc

# ============ Heat math ===============
# Heat input (should match given totals)
Q_pre = Q_preheater_total
Q_boiler = Q_boiler_total
Q_in = Q_pre + Q_boiler

eta_th_water = W_net_water / Q_in
orc_Q_in = (Hc - Hb) * m_dot_orc
eta_th_orc = W_net_orc / orc_Q_in
eta_total = (W_net_water + W_net_orc) / Q_in

# ============ Cooling math ===============

Q_out_orc = (Hd - Ha) * m_dot_orc
H_river = PropsSI('H','T',T_river,'P',P_river,fluid)
H_al_water = H_river + Q_out_orc / m_dot_river
T_al_water = PropsSI('T','H',H_al_water,'P',P_river,fluid)

# ================= H2 Theoretical Work ================
# ==== H2 output
M_Al = 26.98      # g/mol
M_H2 = 2.016      # g/mol
LHV_H2 = 120e6    # J/kg heating value of H2
eff_fc = 0.5 # According to us dept of energy

# H2 produced per kg of Al:
kg_H2_per_kg_Al = (3 * M_H2) / (2 * M_Al)  # around 0.1117 kg

# If m_al is in [kg/s], then H2 mass flow is:
m_H2 = m_al * kg_H2_per_kg_Al  # kg/s of H2, which is still 0.1117kg/s

# Chemical energy flow rate in hydrogen (W = J/s)
Q_dot_H2 = m_H2 * LHV_H2  *eff_fc      # J/s

# Convert to kW
Q_dot_H2_kW = Q_dot_H2 / 1000  # kW 

work_kWh_per_tonne = (W_net_water + W_net_orc + Q_dot_H2) * 1000 / (3.6e6 *m_al)
print(f"Work per tonne of Al: {work_kWh_per_tonne:.2f} kWh/tonne")
print(f"Thermal efficiency of Steam Rankine Cycle: {eta_th_water*100:.2f} %")
print(f"Thermal efficiency of Organic Rankine Cycle: {eta_th_orc*100:.2f} %")
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
