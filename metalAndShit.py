# We should maybe write here what we did for the reactor

"""
The TA basically did the same thing that I did, but went way overboard to make it super complicated...

"""
# Okay so now we fix this shit
import math as mt
import numpy as np

def enthalpy_alumina_highT(T, P_actual=None):
    """
    WE'RE FUCKING USING A NASA EQUATION
    """
    # 500-1200K range
    hf0 = -1675700  
    a1 = -604209
    a2 = 0
    a3 = 14.75481
    a4 = 0.000827
    a5 = 0
    a6 = 0
    a7 = 0
    b1 = -207924

    # Enthalpy calculation
    h = (
        hf0
        - a1 / T
        + a2 * np.log(T)
        + a3 * T
        + (a4 * T**2) / 2
        + (a5 * T**3) / 3
        + (a6 * T**4) / 4
        + (a7 * T**5) / 5
        + b1
    )
    if P_actual is not None:
        # Add pressure correction given by teacher :)
        specific_volume = 2.5575e-5  # from elgoogle
        P_ref = 101325  # Pa
        deltaH_pressure = specific_volume * (P_actual - P_ref)  # J/mol
        h += deltaH_pressure

    return h  # [J/mol] !



def enthalpy_alumina_lowT(T, P_actual=None):
    hf0 = -1675700
    a1 = -604209
    a2 = 0
    a3 = 14.75481
    a4 = 0.000827
    a5 = 0
    a6 = 0
    a7 = 0
    b1 = -207924

    # Calculate enthalpy at temperature T
    h = (
        hf0
        - a1 / T
        + a2 * np.log(T)
        + a3 * T
        + (a4 * T**2) / 2
        + (a5 * T**3) / 3
        + (a6 * T**4) / 4
        + (a7 * T**5) / 5
        + b1
    )

    if P_actual is not None:
        # Add pressure correction given by teacher :)
        specific_volume = 2.5575e-5  # from elgoogle
        P_ref = 101325  # Pa
        deltaH_pressure = specific_volume * (P_actual - P_ref)  # J/mol
        h += deltaH_pressure

    return h  # [J/mol]

def enthalpy_aluminum(T, P_actual=None):
    hf0 = 0
    a1 = -62518.1
    a2 = 634.3934
    a3 = -0.71319
    a4 = 0.010887
    a5 = -1.5E-05
    a6 = 9.96E-09
    a7 = -1.8E-12 
    b1 = -3985.44

    # Calculate enthalpy at temperature T
    h = (
        hf0
        - a1 / T
        + a2 * np.log(T)
        + a3 * T
        + (a4 * T**2) / 2
        + (a5 * T**3) / 3
        + (a6 * T**4) / 4
        + (a7 * T**5) / 5
        + b1
    )

    if P_actual is not None:
        # Add pressure correction given by teacher :)
        specific_volume = 2.5575e-5  # from elgoogle
        P_ref = 101325  # Pa
        deltaH_pressure = specific_volume * (P_actual - P_ref)  # J/mol
        h += deltaH_pressure

    return h  # [J/mol]

# =================== Enthalpy of reactants
# Water



# ALUMINUMMMMMMMM
P_reactor = 300e5 # PA






# ==================== Energy dissipated from reaction ========================
molar_mass_Al = 26.98  # g/mol
energy_per_2mol_Al = 818.21  # kJ

energy_per_mol_Al = energy_per_2mol_Al / 2  

mol_per_kg = 1000 / molar_mass_Al
energy_per_kg_Al_kJ = energy_per_mol_Al * mol_per_kg

energy_per_kg_Al_MJ = energy_per_kg_Al_kJ / 1000 #15.16
#print(energy_per_kg_Al_MJ)
# So we have the energy if we keep everything at STAP which is 25 + 273.15 = 298.15K
# But our alumina, our H2 is at 500C which is 773.15K, so we need to look into the energy that will be used to heat them up

# Heating up my stuff
hAlumina_low = enthalpy_alumina_lowT(25+273.15)
hALumina_high = enthalpy_alumina_highT(500+273.15)
delHAlumina = hALumina_high - hAlumina_low
print(delHAlumina)
