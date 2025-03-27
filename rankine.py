import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# All the shit we can play with

"""
Stuff I'm noting while playing with this
So UK power plants produce about 450kg/s according to wikipedia for a 2000MW plant
Interpolate to our 100MW, this is aboud 22kg/s, by increments of 5kg/s, we gain around 300kwh/tonne
however, at 5kg/s, the heat transfer in the boiler is weird  and the final temperature is 1600K which doesn't make any sense
what I think we could do is take the steam we don't bleed and reheat using a small fraction of the boiler

We can gain net power by doing P_condenser to 0.01e6

So, atm the temp of the superheated vapour is 900k

I also forced a H at 6 which is then used to find the bleed ratio

"""

# === Parameters ===
m_dot = 18  # kg/s of water/steam
m_al = 1   # kg/s of aluminum reacting

# Heat available per kg of aluminum
q_preheater_per_al = 0.7e6  # J/kg Al, I reduced this because it didn't make sense if it was 1.7MJ
q_boiler_per_al = 12.5e6  - 0  # J/kg Al
q_reheat_per_al = 1e6 # J/kg, so about 2.5e6

# Total heat available per second
Q_preheater_total = q_preheater_per_al * m_al
Q_boiler_total = q_boiler_per_al * m_al
Q_reheat_total = q_reheat_per_al * m_al

# Cycle pressures
P_condenser = 0.1e6     # Pa 
P_boiler = 4e6          # Pa
P_bleed = 1e6           # Pa

# Other parameters
T_boiler_out = 773.15   # K, 500C which is the reactor temperature at ss
eff_turbine = 0.8 # from online this should make sense
eff_pump = 0.8
fluid = 'Water'

# ----------------------------------------------- Logic of the cycle condenser - pump- regen - pre-heat - boiler - turbine 1 - turbine 2

# === State 1: After condenser ===
P1 = P_condenser
H1 = PropsSI('H', 'P', P1, 'Q', 0, fluid)
S1 = PropsSI('S', 'P', P1, 'Q', 0, fluid)

# === State 2: After pump === 
"""
Problem here as we bump up the pressure in the boiler phase
We would need a compressor somewhere
"""

P2 = P_bleed
H2s = PropsSI('H', 'P', P2, 'S', S1, fluid)
H2 = H1 + (H2s - H1) / eff_pump # Isentropic eff
S2 = PropsSI('S', 'P', P2, 'H', H2, fluid)
T2 = PropsSI('T', 'P', P2, 'H', H2, fluid)

# === State 3: After regenerator ===, fix this after, no estimate
print(T_boiler_out)
H3 = PropsSI('H', 'P', P2, 'T', T_boiler_out - 200, fluid)  # Estimate
T3 = PropsSI('T', 'P', P2, 'H', H3, fluid)
S3 = PropsSI('S', 'P', P2, 'H', H3, fluid)

# === State 4: After preheater ===, we gotta use a q that makes it a sat vapour
H4_ideal = PropsSI('H', 'P', P2, 'Q', 1, fluid)
Q_reheat_total = (H4_ideal - H3) * m_dot
T4 = PropsSI('T', 'P', P2, 'H', H4_ideal, fluid)
S4 = PropsSI('S', 'P', P2, 'H', H4_ideal, fluid)

# Update the actual boiler energy
q_reheat_per_al = Q_reheat_total / m_al
print(f"The energy used to eat to sat vapour {q_reheat_per_al:.2f}")
q_boiler_per_al = 12.5e6 - q_reheat_per_al
Q_boiler_total = q_boiler_per_al * m_al

# ====== Compressor between pre-heat and boiler
P_comp = P_boiler
H_comps = PropsSI('H', 'P', P_comp, 'S', S4, fluid)
H_comp = H4_ideal + (H_comps - H4_ideal) / eff_pump
S_comp = PropsSI('S', 'P', P_comp, 'H', H_comp, fluid)

# === State 5: After boiler ===
H5 = H_comp + Q_boiler_total / m_dot
P5 = P_boiler

try:
    T5 = PropsSI('T', 'P', P5, 'H', H5, fluid)
    S5 = PropsSI('S', 'P', P5, 'H', H5, fluid)
except ValueError:
    print("X Error: Boiler outlet enthalpy too high. Increase m_dot or reduce Q_boiler_total. Or reduce cc temperature, but modify q out")
    exit()

# === Turbine 1: Expand to bleed pressure ===
H6s = PropsSI('H', 'P', P_bleed, 'S', S5, fluid)
H6 = H5 - eff_turbine * (H5 - H6s)
T6 = PropsSI('T', 'P', P_bleed, 'H', H6, fluid)
S6 = PropsSI('S', 'P', P_bleed, 'H', H6, fluid)

# Bleed part goes to regenerator — assume mass fraction x
x_bleed = (H3 - H2) / (H6 - H2)

"""
# Reheat stage between Turbine 1 and Turbine 2, we'll call it 6.5
H6_5 = H6 + Q_reheat_total / ((1-x_bleed)*m_dot)
P6_5 = P_bleed
S6_5 = PropsSI('S', 'P', P6_5, 'H', H6_5, fluid)"""

# === Turbine 2: Remaining steam to condenser ===
H7s = PropsSI('H', 'P', P_condenser, 'S', S6, fluid)
H7 = H6 - eff_turbine * (H6 - H7s)
T7 = PropsSI('T', 'P', P_condenser, 'H', H7, fluid)
S7 = PropsSI('S', 'P', P_condenser, 'H', H7, fluid)

# === Energetics ===
# Turbine work
W_turb1 = (H5 - H6) * m_dot
W_turb2 = (1 - x_bleed) * m_dot * (H6 - H7)

# Pump work
W_pump = m_dot * (H2 - H1)

# Heat input (should match given totals)
Q_pre = Q_preheater_total
Q_boiler = Q_boiler_total
Q_in = Q_pre + Q_boiler

# Net work & efficiency
W_net = W_turb1 + W_turb2 - W_pump
eta_th = W_net / Q_in


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


# === Output ===
print("\n--- Regenerative Rankine Cycle ---")
print(f"Temp after boiler:  {T5:.2f} K")
print(f"Added Q during preheat {q_boiler_per_al/(1e6):.2f}")
print(f"Aluminum mass flow: {m_al:.2f} kg/s")
print(f"Steam mass flow:    {m_dot:.2f} kg/s")
print(f"Turbine 1 work:     {W_turb1/1000:.2f} kW")
print(f"Turbine 2 work:     {W_turb2/1000:.2f} kW")
print(f"Pump work:          {W_pump/1000:.2f} kW")
print(f"Heat input:         {Q_in/1000:.2f} kW")
print(f"Net work:           {W_net/1000:.2f} kW")
print(f"Thermal efficiency: {eta_th*100:.2f} %")
print(f"Bleed fraction:     {x_bleed*100:.1f} %")
print(f"Energy in H2 Fuel Cell: {Q_dot_H2_kW:.2f} kW")
work_kWh_per_tonne = (W_net + Q_dot_H2) * 1000 / 3.6e6
print(f"Work per tonne of Al: {work_kWh_per_tonne:.2f} kWh/tonne")

# === T-s Diagram ===
T_vals = [PropsSI('T', 'H', h, 'P', p, fluid) for h, p in zip(
    [H1, H2, H3, H4_ideal, H5, H6, H7, H1],
    [P1, P2, P2, P2, P5, P_bleed, P_bleed, P1, P1]
)]
S_vals = [PropsSI('S', 'H', h, 'P', p, fluid) for h, p in zip(
    [H1, H2, H3, H4_ideal, H5, H6, H7, H1],
    [P1, P2, P2, P2, P5, P_bleed, P_bleed, P1, P1]
)]

plt.figure(figsize=(8, 6))
plt.plot(S_vals, T_vals, 'o-', lw=2)

# Annotate each state
for i, (s, t) in enumerate(zip(S_vals, T_vals), start=1):
    plt.text(s + 5, t, str(i), fontsize=10, verticalalignment='bottom')

plt.xlabel("Entropy [J/kg·K]")
plt.ylabel("Temperature [K]")
plt.title("Regenerative Rankine Cycle T–s Diagram")
plt.grid(True)
plt.tight_layout()
plt.show()
