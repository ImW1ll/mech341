import math
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

# =============== Helper function ========================

def cp_air_kj_per_kmolK(T):
    a = 28.11
    b = 0.1967e-2
    c = 0.4802e-5
    d = -1.966e-9   
    return a + b*T + c*(T**2) + d*(T**3)

def cp_air_j_per_kgK(T):
    # Molar mass of air [kg/kmol] cos the book give it in kmol?...
    cp_kj_per_kmolK = cp_air_kj_per_kmolK(T)
    cp_j_per_kgK = (cp_kj_per_kmolK * 1000.0) / 0.02897
    return cp_j_per_kgK

def cp_N2_j_per_kgK(T):
    if 100 <= T <= 500:
        a = 28.98641
        b = 1.853978
        c = -9.648459
        d = 16.63537
        e = 0.000117
    elif 500 <= T <= 2000:
        a = 19.50583
        b = 19.88705
        c = -8.598535
        d = 1.369784
        e = 0.527601
    else:
        raise AssertionError("Dude you stupid")
    return (a + b*T + c*T**2 + d*T**3 + e/T**2) / 0.0280134

def cp_O2_j_per_kgK(T):
    if 100 <= T <= 700:
        a = 31.32234
        b = -20.23531
        c = 57.86644
        d = -36.50624
        e = -0.007374
    elif 700 <= T <= 2000:
        a = 30.03235
        b = 8.772972
        c = -3.988133
        d = 0.788313
        e = -0.741599
    else:
        AssertionError("Dude you stupid part 2")
    return (a + b*T + c*T**2 + d*T**3 + e/T**2) / 0.032

# =============== Main Body ========================

# Notes: MOLE fraction of air: 0.21 O2: 0.79 N2. MASS fraction of air: 0.232 O2: 0.768 N2

def brayton_cycle_h2(
    m_H2,               # mass flow hydrogen
    TB0,                # Temperature of Hydrogen before heat exchange (K)
    PB0,                # Pressure of hydrogen before heat exchange (Pa)
    lmdba=2.1,                # Air-to-fuel ratio for mixing
    fuel_name = 'Hydrogen',
    eta_turb = 0.88,
    eta_comp = 0.88,
    ):

    # --- Molar Masses --- 
    M_H2 = 0.002016  # kg/mol
    M_O2 = 0.032
    M_H2O = 0.018015 
    M_N2 = 0.0280134

    # --- State B0: Before the reheat
    SB0 = PropsSI('S','P',PB0,'T',TB0,fuel_name)
    HB0 = PropsSI('H', 'P', PB0, 'T', TB0, fuel_name)

    # --- State B1: Hydrogen exiting reheat ---
    TB1 = 340+273.15
    PB1 = PB0
    SB1 = PropsSI('S','T',TB1,'P',PB1,fuel_name)
    HB1 = PropsSI('H','T',TB1,'P',PB1,fuel_name)
    # reheat = m_H2*(HB0 - HB1) 2370615.719104707, so way down on this one


    # --- State B2: Hydrogen after expansion ---
    PB2 = 3e6                                      # Teacher said 25-30bar with 1800K makes sense
    HB2s = PropsSI('H','P',PB2,'S',SB1,fuel_name)
    HB2 = HB1 + (HB2s - HB1) * eta_turb
    SB2 = PropsSI('S','P',PB2,'H',HB2,fuel_name)
    TB2 = PropsSI('T','P',PB2,'H',HB2,fuel_name)

    # --- State B3: Air before compression ---
    # Air is ~21% O₂, ~78% N₂ by mole fraction, STAP
    PB3 = 1e5 
    TB3 = 25 + 273.15
    PB3_O2 = 0.21 * PB3 # Partials
    PB3_N2 = 0.79 * PB3
    HB3_O2 = PropsSI('H', 'P', PB3_O2, 'T', TB3, 'Oxygen')
    HB3_N2 = PropsSI('H', 'P', PB3_N2, 'T', TB3, 'Nitrogen')
    SB3_O2 = PropsSI('S', 'P', PB3_O2, 'T', TB3, 'Oxygen')
    SB3_N2 = PropsSI('S', 'P', PB3_N2, 'T', TB3, 'Nitrogen')
    HB3 = 0.232 * HB3_O2 + 0.768 * HB3_N2
    SB3 = 0.232 * SB3_O2 + 0.768 * SB3_N2

    # --- State B4: Air After compression ---
    # No mass flow here, It's found through the lmdba controller at state 5 :)
    PB4 = PB2
    PB4_O2 = 0.21 * PB4
    PB4_N2 = 0.79 * PB4
    HB4s_O2 = PropsSI('H', 'P', PB4_O2, 'S', SB3_O2, 'Oxygen')
    HB4s_N2 = PropsSI('H', 'P', PB4_N2, 'S', SB3_N2, 'Nitrogen')
    HB4_O2 = HB3_O2 + (HB4s_O2 - HB3_O2) / eta_comp
    HB4_N2 = HB3_N2 + (HB4s_N2 - HB3_N2) / eta_comp
    SB4_O2 = PropsSI('S', 'P', PB4_O2, 'H', HB4_O2, 'Oxygen')
    SB4_N2 = PropsSI('S', 'P', PB4_N2, 'H', HB4_N2, 'Nitrogen')
    HB4 = 0.232 * HB4_O2 + 0.768 * HB4_N2
    SB4 = 0.232 * SB4_O2 + 0.768 * SB4_N2
    TB4 = HB4 / cp_air_j_per_kgK(799.5) # Guessed and checked

    # --- State B5: After Combustion --- we cooking here
    n_H2 = m_H2 / M_H2
    n_O2_stoich = 0.5 * n_H2
    n_O2_supplied = lmdba * n_O2_stoich
    n_O2_excess = n_O2_supplied - n_O2_stoich
    n_H2O = n_H2  # watah produced
    n_N2 = n_O2_supplied * (79.0 / 21.0) # 3.76

    # Convert mol/s to kg/s
    m_O2_supplied = n_O2_supplied * M_O2
    m_N2 = n_N2 * M_N2
    m_H2_in = n_H2 * M_H2

    # Enthalpy of reactants (sensible only)
    H_reactants = (m_O2_supplied * HB4_O2) + (m_N2 * HB4_N2) + (m_H2_in * HB2)

    # teacher gave equation with watah at liquid :(
    # If final is vapor, 2.42 not 2.46 as given by teacher
    DeltaH_comb_per_mol = 2.42e5 # that is j/mol
    H_rxn = n_H2 * DeltaH_comb_per_mol

    n_total = n_H2O + n_O2_excess + n_N2
    X_H2O = n_H2O / n_total
    X_O2  = n_O2_excess / n_total
    X_N2  = n_N2 / n_total

     # Partial pressures again
    PB5 = PB4
    PB5_H2O = X_H2O * PB5
    PB5_O2  = X_O2  * PB5
    PB5_N2  = X_N2  * PB5

    def energy_balance(T_guess):
        h_H2O = PropsSI('H', 'P', PB5_H2O, 'T', T_guess, 'Water')
        h_O2  = PropsSI('H', 'P', PB5_O2, 'T', T_guess, 'Oxygen')
        h_N2  = PropsSI('H', 'P', PB5_N2, 'T', T_guess, 'Nitrogen')
        
        m_H2O = n_H2O * 0.018015      # kg/s
        m_O2_excess = n_O2_excess * M_O2
        m_N2_out = n_N2 * M_N2
        
        # Sensible enthalpy of products:
        H_products = (m_H2O * h_H2O) + (m_O2_excess * h_O2) + (m_N2_out * h_N2)
        
        # Final Balance
        return H_products - (H_reactants + H_rxn)

    T_B5_guess = 1800
    TB5 = fsolve(energy_balance, T_B5_guess)[0]
    print(f"Temperature of the products after combustion: {TB5:.1f} K")
    
   

    # Enthalpy
    HB5_H2O = PropsSI('H','P',PB5_H2O,'T',TB5,'Water')
    HB5_O2  = PropsSI('H','P',PB5_O2,'T',TB5,'Oxygen')
    HB5_N2  = PropsSI('H','P',PB5_N2,'T',TB5,'Nitrogen') 
    HB5 = (HB5_H2O * X_H2O * M_H2O + HB5_N2 * X_N2 * M_N2 + HB5_O2 * X_O2 * M_O2) / (X_H2O * M_H2O + X_N2 * M_N2 + X_O2 * M_O2)   
    SB5_H2O = PropsSI('S','P',PB5_H2O,'T',TB5,'Water')
    SB5_O2  = PropsSI('S','P',PB5_O2,'T',TB5,'Oxygen')
    SB5_N2  = PropsSI('S','P',PB5_N2,'T',TB5,'Nitrogen')
    SB5 = (SB5_H2O * X_H2O * M_H2O + SB5_N2 * X_N2 * M_N2 + SB5_O2 * X_O2 * M_O2) / (X_H2O * M_H2O + X_N2 * M_N2 + X_O2 * M_O2)
    
    # --- State B6: After Expansion --- 

    PB6 = 101325 # 1atm basically
    PB6_H2O = X_H2O * PB6
    PB6_O2  = X_O2  * PB6
    PB6_N2  = X_N2  * PB6

    # Ideal isentropic expansion (still use entropy)
    HB6s_H2O = PropsSI('H', 'P', PB6_H2O, 'S', SB5_H2O, 'Water')
    HB6s_O2  = PropsSI('H', 'P', PB6_O2, 'S', SB5_O2,  'Oxygen')
    HB6s_N2  = PropsSI('H', 'P', PB6_N2, 'S', SB5_N2,  'Nitrogen')
    
    # Real expansion with efficiency
    HB6_H2O = HB5_H2O + eta_turb * (HB6s_H2O - HB5_H2O)
    HB6_O2  = HB5_O2  + eta_turb * (HB6s_O2  - HB5_O2)
    HB6_N2  = HB5_N2  + eta_turb * (HB6s_N2  - HB5_N2)
    HB6 = (HB6_H2O * X_H2O * M_H2O + HB6_N2 * X_N2 * M_N2 + HB6_O2 * X_O2 * M_O2) / (X_H2O * M_H2O + X_N2 * M_N2 + X_O2 * M_O2)
    SB6_H2O = PropsSI('S','P',PB6_H2O,'H',HB6_H2O,'Water')
    SB6_O2  = PropsSI('S','P',PB6_O2,'H',HB6_O2,'Oxygen')
    SB6_N2  = PropsSI('S','P',PB6_O2,'H',HB6_N2,'Nitrogen')
    SB6 = (SB6_H2O * X_H2O * M_H2O + SB6_N2 * X_N2 * M_N2 + SB6_O2 * X_O2 * M_O2) / (X_H2O * M_H2O + X_N2 * M_N2 + X_O2 * M_O2)
    TB6 = HB6 / cp_air_j_per_kgK(600) # idk if this is right

    # And then idk mother nature will manage the rest

    # --- Work! Work! Work! ---
    m_H2O = n_H2O * 0.018015      # kg/s)
    m_O2_excess = n_O2_excess * M_O2
    m_N2_out = n_N2 * M_N2
    m_air = m_O2_supplied + m_N2

    #print(f"HB4 {HB4}")
    #print(f"HB3 {HB3}")

    W_comp = m_air * (HB4 - HB3)
    print(f"Compression work: {W_comp:.2f}")
    W_turb = m_H2O * (HB5_H2O - HB6_H2O) + m_O2_excess*(HB5_O2 - HB6_O2) + m_N2_out*(HB5_N2 - HB6_N2)
    W_turb2 = m_H2*(HB1-HB2)
    print(f"Turbine work out: {W_turb:.2f}")
    work_comp = (W_turb + W_turb2) / W_comp
    print(f"Work over comp work: {work_comp:.2f}")

    W_net = W_turb - W_comp + W_turb2

    States = dict()
    States['B0'] = (PB0,TB0,HB0,SB0)
    States['B1'] = (PB1,TB1,HB1,SB1)
    States['B2'] = (PB2,TB2,HB2,SB2)
    States['B3'] = (PB3,TB3,HB3,SB3)
    States['B4'] = (PB4,TB4,HB4,SB4)
    States['B5'] = (PB5,TB5,HB5,SB5)
    States['B6'] = (PB6,TB6,HB6,SB6)
    States['Work'] = W_net

    return States

 # Testing Unit
if __name__ == "__main__":

    m_H2 = 1.0087472201630838

    # Run the Brayton cycle function
    Dog = brayton_cycle_h2(
        m_H2,500+273.15,3e7,1.92
        )
