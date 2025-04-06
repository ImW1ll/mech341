# We should maybe write here what we did for the reactor

"""
The TA basically did the same thing that I did, but went way overboard to make it super complicated...

"""
# Okay so now we fix this shit
import math as mt
import numpy as np
from CoolProp.CoolProp import PropsSI

def enthalpy_aluminum(T, P_actual=None):
    # Constants
    cp_aluminum = 900  # J/kg/K, approx solid
    hf0_aluminum = 0   # J/kg 
    T_ref = 298.15     # Reference temperature [K]
    P_ref = 101325     # Pa

    # Sensible enthalpy [J/kg]
    hs = cp_aluminum * (T - T_ref)

    # Total enthalpy
    h = hf0_aluminum + hs

    if P_actual is not None and P_actual != P_ref:
        specific_volume = 1/2700  # m³/kg
        deltaH_pressure = specific_volume * (P_actual - P_ref)  # J/kg
        h += deltaH_pressure

    return h

def enthalpy_alumina_lowT(T, P_actual=None):
    """
    So, for some fucking reason, he uses Nasa on alumina but not on Aluminum because the hf0 on alumina is not zero

    This is for 200 to 500 max
    """
    # Constants
    R_u = 8.314462618  # J/mol/K
    M_alumina = 0.1019613  # kg/mol
    P_ref = 101325  # Pa

    # Formation enthalpy (per kg)
    hf0_alumina_mol = -1675000  # J/mol
    hf0_alumina = hf0_alumina_mol / M_alumina  # J/kg

    # NASA polynomial coefficients (valid 500–1200 K)
    a1 = -5391550
    a2 = 103667.7
    a3 = -817.323
    a4 = 3.388259
    a5 = -0.00751
    a6 = 8.66e-6
    a7 = -4.1e-9
    b1 = -666013

    # Enthalpy calculation [J/mol]
    h_mol = R_u * (
        - a1 / T
        + a2 * np.log(T)
        + a3 * T
        + (a4 * T**2) / 2
        + (a5 * T**3) / 3
        + (a6 * T**4) / 4
        + (a7 * T**5) / 5
        + b1
    )

    # Convert to [J/kg]
    h_total = h_mol / M_alumina

    if P_actual is not None and P_actual != P_ref:
        specific_volume_molar = 2.56e-5  # [m³/mol]
        deltaH_pressure_mol = specific_volume_molar * (P_actual - P_ref)  # [J/mol]
        deltaH_pressure_kg = deltaH_pressure_mol / M_alumina  # [J/kg]
        h_total += deltaH_pressure_kg

    # Sensible enthalpy = total - formation enthalpy
    h_sensible = h_total - hf0_alumina

    return h_sensible

def enthalpy_alumina_highT(T, P_actual=None):
    """
    This is for 500K to 1200K, so stay here all the time actually
    """
    # Constants
    R_u = 8.314462618  # J/mol/K
    M_alumina = 0.1019613  # kg/mol
    P_ref = 101325  # Pa

    # Formation enthalpy (per kg)
    hf0_alumina_mol = -1675000  # J/mol
    hf0_alumina = hf0_alumina_mol / M_alumina  # J/kg

    # NASA polynomial coefficients (valid 500–1200 K)
    a1 = -604209
    a2 = 0
    a3 = 14.75481
    a4 = 0.000827
    a5 = 0
    a6 = 0
    a7 = 0
    b1 = -207924

    # Enthalpy calculation [J/mol]
    h_mol = R_u * (
        - a1 / T
        + a2 * np.log(T)
        + a3 * T
        + (a4 * T**2) / 2
        + (a5 * T**3) / 3
        + (a6 * T**4) / 4
        + (a7 * T**5) / 5
        + b1
    )

    # Convert to [J/kg]
    h_total = h_mol / M_alumina

    if P_actual is not None and P_actual != P_ref:
        specific_volume_molar = 2.56e-5  # [m³/mol]
        deltaH_pressure_mol = specific_volume_molar * (P_actual - P_ref)  # [J/mol]
        deltaH_pressure_kg = deltaH_pressure_mol / M_alumina  # [J/kg]
        h_total += deltaH_pressure_kg

    # Sensible enthalpy = total - formation enthalpy
    h_sensible = h_total - hf0_alumina

    return h_sensible



# ==================== Energy dissipated from reaction ========================


hf0 = -1675700

def aluminum_heat_rate_out(m_Al,m_H2O_rankine): 
    
    p_reactor = 3e7
    t_reactor = 500 + 273.19 #K
    p_rankine = 25e6
    t_ambient = 298.13
    p_ambient = 101.3e3 # Pa

    # ========= Enthalpy of Reactants =========:
    
    # --- Al ---
    M_Al = 26.98 # kg/kmol
    T_Al = t_ambient # K at SATP
    P_Al = p_reactor
    H_Al_0 = 0
    H_Al = enthalpy_aluminum(T_Al) # J/kg
    H_Al_ref = 0
    H_s = H_Al - H_Al_ref
    H_dot_Al_s = H_s * m_Al
    H_dot_Al_f0 = 0

    # --- Reactant H2O ---
    M_H2O = 18.015 # kg/kmol
    T_H2O = t_ambient
    P_H2O = p_reactor
    m_H2O = m_Al * (3 * M_H2O) / (2 * M_Al)
    H_H2O_0 = -285830 * 1000 / M_H2O # J/kg, from the book
    H_H2O = PropsSI('H','T',T_H2O,'P',P_H2O)
    H_H2O_ref = PropsSI('H','T',t_ambient,'P',p_ambient,'Water') # J/kg
    H_H2O_s = H_H2O - H_H2O_ref
    H_dot_H2O_s = H_H2O_s * m_H2O
    H_dot_H2O_f0 = H_H2O_0 * m_H2O

    # --- Rankine H2O ---
    T_H2O_R = 2 #!!!!!!
    P_H2O_R = p_rankine
    H_H2O_R = PropsSI('H','T',T_H2O_R,'P',P_H2O_R,'Water')
    H_H2O_R_s = H_H2O_R - H_H2O_ref
    H_dot_H2O_R_s = H_H2O_R_s * m_H2O_rankine
    H_dot_H2O_R_f0 = H_H2O_0 * m_H2O_rankine

    # ========= Enthalpy of Products =========:

    # --- Al2O3 ---
    M_Al2O3 = 101.96 # kg/kmol
    T_Al2O3 = t_reactor
    m_Al2O3 = m_Al * (M_Al2O3) / (2 * M_Al)
    H_Al2O3_0 = -1675700 / M_Al2O3 # J/kg
    H_Al2O3 = enthalpy_alumina_highT(T_Al2O3, 2.5e7)
    H_Al2O3_ref = enthalpy_alumina_lowT(t_ambient)
    H_Al2O3_s = H_Al2O3 - H_Al2O3_ref
    H_dot_Al2O3_s = H_Al2O3_s * m_Al2O3
    H_dot_Al2O3_f0 = H_Al2O3_0 * m_Al2O3

    # --- H2 ---
    M_H2 = 2.015 # kg/kmol
    T_H2 = t_reactor
    P_H2 = p_reactor # Assume stoichiometric conversion
    m_H2 = m_Al * (3 * M_H2) / (2 * M_Al)
    H_H2_0 = 0
    H_H2 = PropsSI('H','T',T_H2,'P',P_H2,'Hydrogen')
    H_H2_ref = PropsSI('H','T',t_ambient,'P',p_ambient,'Hydrogen')    
    H_H2_s = H_H2 - H_H2_ref
    H_dot_H2_s = H_H2_s * m_H2
    H_dot_H2_f0 = H_H2_0 * m_H2

    # --- Rankine H2O (out) ---
    T_H2O_RO = t_reactor
    P_H2O_RO = p_rankine
    H_H2O_RO = PropsSI('H','T',t_reactor,'P',p_rankine,'Water')
    H_H2O_RO_s = H_H2O_RO - H_H2O_ref
    H_dot_H2O_RO_s = H_H2O_RO_s * m_H2O_rankine
    H_dot_H2O_RO_f0 = H_H2O_0 * m_H2O_rankine

    # ========= First Law Analysis =========:

    H_dot_in_s = H_dot_Al_s + H_dot_H2O_s
    H_dot_in_f0 = H_dot_Al_f0 + H_dot_H2O_f0

    H_dot_out_s = H_dot_Al2O3_s + H_dot_H2_s
    H_dot_out_f0 = H_dot_Al2O3_f0 + H_dot_H2_f0

    return H_dot_out_s - H_dot_in_s # this should be the net work






# checking reheat values
if __name__ == "__main__":
    # First we need to check reheat in alumina
    HAL_high  = enthalpy_alumina_highT(500+273.15, 3e7)
    HAL_low  = enthalpy_alumina_highT(330+273.15, 3e7)
    reheat_q = (HAL_high - HAL_low)*17.005 # that's mass flow alumina for 9 kg of aluminium
    print(reheat_q/1e6)