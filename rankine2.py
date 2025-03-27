"""
Need to redo a solid rewriting of the cycle


"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# === Parameters ===
m_dot = 10  # kg/s of water/steam 2000MW, 100MW, 22.5kg/s
m_al = 1   # kg/s of aluminum reacting


# Heat available per kg of aluminum
q_preheater_per_al = 0.7e6  # J/kg Al, I reduced this because it didn't make sense if it was 1.7MJ
q_reheat_per_al = 0 # J/kg
q_boiler_per_al = 12.5e6 - q_reheat_per_al # Heat to cool reactor to keep it at 500C at SS


# Total heat available per second
Q_preheater_total = q_preheater_per_al * m_al
q_reheat_total = q_reheat_per_al * m_al
Q_boiler_total = q_boiler_per_al * m_al # W

# Cycle pressures
P_condenser = 0.01e6     # Pa 
P_boiler = 3e6          # Pa, should be at operating temp of the turb 1
P_t1out = 0.1e6           # Pa, pressure after first turbine

# Other parameters
T_boiler_out = 773.15   # K, 500C which is the reactor temperature at ss
eff_turbine = 0.8 # from online this should make sense
eff_pump = 0.8
fluid = 'Water'

"""
This should work to get a temp of 500C at the exit of the boiler
"""

# Stage after condenser (2)
P2 = P_condenser
Q2 = 0
H2 = PropsSI('H', 'P', P2, 'Q', Q2, fluid)
S2 = PropsSI('S', 'P', P2, 'Q', Q2, fluid)

# Compressor
P3 = P_boiler
H3s = PropsSI('H', 'P', P3, 'S', S2, fluid)
H3 = H2 + (H3s - H2) / eff_pump

# Stage after boiler
#T5 = T_boiler_out # 500C
P5 = P_boiler # 40bar
H5 = H3 + Q_boiler_total / m_dot
#H5 = PropsSI('H', 'P', P5, 'T', T5, fluid)
S5 = PropsSI('S', 'P', P5, 'T', T5, fluid)

# Stage after first turbine
H6s = PropsSI('H', 'S', S5, 'P', P_t1out, fluid) # Pressure here is arbritary
H6 = H5 + (H6s - H5) / eff_turbine
P6 = P_t1out
S6 = PropsSI('S', 'P', P6, 'H', H6, fluid)

# Stage after second turbine - before condenser
P1 = P_condenser # Arbritrary pressure
H1s = PropsSI('H', 'P', P1, 'S', S6, fluid)
H1 = H6 + (H6 - H1s) / eff_turbine
T1 = PropsSI('T', 'H', H1, 'P', P1, fluid)



# Stage before boiler, reheater
H4 = H5 - Q_boiler_total / m_dot
P4 = P5 # We assume that the boiler is an isobaric process

# Find how much energy we need to reheat
q_reheat = (H4 - H3) * m_dot



# ============= Work math ================
# Work output
W_turb1 = (H5 - H6) * m_dot
W_turb2 = (H6 - H1) * m_dot

W_pump = (H3 - H2) * m_dot

W_net = W_turb1 + W_turb2 - W_pump


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


work_kWh_per_tonne = (W_net) * 1000 / 3.6e6
print(f"Work per tonne of Al: {work_kWh_per_tonne:.2f} kWh/tonne")

