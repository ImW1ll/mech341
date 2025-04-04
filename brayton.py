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
    M_air = 28.97
    cp_kj_per_kmolK = cp_air_kj_per_kmolK(T)
    cp_j_per_kgK = (cp_kj_per_kmolK * 1000.0) / M_air
    return cp_j_per_kgK


# =============== Main Body ========================

def brayton_cycle_h2(
    m_H2,               # mass flow hydrogen
    ):

    # I'll hard code the values and we will correct later
    fuel_name = 'Hydrogen'
    eta_turb = 0.85
    eta_comp = 0.85

    # --- State B0: Before the reheat
    TB0 = 500 + 273.15 # K
    PB0 = 3e7          # Pa
    HB0 = PropsSI('H', 'P', PB0, 'T', TB0, fuel_name)

    # --- State B1: Hydrogen exiting reheat ---
    TB1 = 100+273.15
    PB1 = 3e7
    SB1 = PropsSI('S','T',TB1,'P',PB1,fuel_name)
    HB1 = PropsSI('H','T',TB1,'P',PB1,fuel_name) # This gives a 5879213.584799376 J which is pretty good, for the reheat

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
    PB3_N2 = 0.78 * PB3
    HB3_O2 = PropsSI('H', 'P', PB3_O2, 'T', TB3, 'Oxygen')
    HB3_N2 = PropsSI('H', 'P', PB3_N2, 'T', TB3, 'Nitrogen')
    SB3_O2 = PropsSI('S', 'P', PB3_O2, 'T', TB3, 'Oxygen')
    SB3_N2 = PropsSI('S', 'P', PB3_N2, 'T', TB3, 'Nitrogen')
    HB3_air = 0.21 * HB3_O2 + 0.78 * HB3_N2

    # --- State B4: Air After compression ---
    # No mass flow here, It's found through the lambda controller at state 5 :)
    PB4 = PB2
    PB4_O2 = 0.21 * PB4
    PB4_N2 = 0.78 * PB4
    HB4s_O2 = PropsSI('H', 'P', PB4_O2, 'S', SB3_O2, 'Oxygen')
    HB4s_N2 = PropsSI('H', 'P', PB4_N2, 'S', SB3_N2, 'Nitrogen')
    HB4_O2 = HB3_O2 + (HB4s_O2 - HB3_O2) / eta_comp
    HB4_N2 = HB3_N2 + (HB4s_N2 - HB3_N2) / eta_comp
    TB4_O2 = PropsSI('T', 'P', PB4_O2, 'H', HB4_O2, 'Oxygen')
    TB4_N2 = PropsSI('T', 'P', PB4_N2, 'H', HB4_N2, 'Nitrogen')
    HB4_air = 0.21 * HB4_O2 + 0.78 * HB4_N2


    # --- State B5: After Combustion --- we cooking here
    lambda_ = 1.3   # air/fuel
    M_H2 = 0.002016  # kg/mol
    M_O2 = 0.032
    M_N2 = 0.0280134 * 2
    n_H2 = m_H2 / M_H2
    n_O2_stoich = 0.5 * n_H2
    n_O2_supplied = lambda_ * n_O2_stoich
    n_O2_excess = n_O2_supplied - n_O2_stoich
    n_H2O = n_H2  # watah produced
    n_N2 = n_O2_supplied * (78.0 / 21.0)

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

    def energy_balance(T_guess):
        h_H2O = PropsSI('H', 'P', PB4, 'T', T_guess, 'Water')
        h_O2  = PropsSI('H', 'P', PB4, 'T', T_guess, 'Oxygen')
        h_N2  = PropsSI('H', 'P', PB4, 'T', T_guess, 'Nitrogen')
        
        m_H2O = n_H2O * 0.018015      # kg/s)
        m_O2_excess = n_O2_excess * M_O2
        m_N2_out = n_N2 * M_N2
        
        # Sensible enthalpy of products:
        H_products = (m_H2O * h_H2O) + (m_O2_excess * h_O2) + (m_N2_out * h_N2)
        
        # Final Balance
        return H_products - (H_reactants + H_rxn)

    T_B5_guess = 1800
    TB5 = fsolve(energy_balance, T_B5_guess)[0]
    print(f"Temperature of the products after combustion: {TB5:.1f} K")
    # Partial pressures again
    n_total = n_H2O + n_O2_excess + n_N2
    X_H2O = n_H2O / n_total
    X_O2  = n_O2_excess / n_total
    X_N2  = n_N2 / n_total

    PB5_H20 = X_H2O * PB4
    PB5_O2  = X_O2  * PB4
    PB5_02  = X_N2  * PB4

    # Enthalpy
    HB5_H20 = PropsSI('H','P',PB5_H20,'T',TB5,'Water')
    SB5_H20 = PropsSI('S','P',PB5_H20,'T',TB5,'Water')
    HB5_O2  = PropsSI('H','P',PB5_O2,'T',TB5,'Oxygen')
    SB5_O2  = PropsSI('S','P',PB5_O2,'T',TB5,'Oxygen')
    HB5_N2  = PropsSI('H','P',PB5_02,'T',TB5,'Nitrogen')
    SB5_N2  = PropsSI('S','P',PB5_02,'T',TB5,'Nitrogen')

    # --- State B6: After Expansion --- 
    PB6 = 101325 # 1atm basically
    # Ideal isentropic expansion (still use entropy)
    HB6s_H2O = PropsSI('H', 'P', PB6, 'S', SB5_H20, 'Water')
    HB6s_O2  = PropsSI('H', 'P', PB6, 'S', SB5_O2,  'Oxygen')
    HB6s_N2  = PropsSI('H', 'P', PB6, 'S', SB5_N2,  'Nitrogen')

    # Real expansion with efficiency
    HB6_H2O = HB5_H20 + eta_turb * (HB6s_H2O - HB5_H20)
    HB6_O2  = HB5_O2  + eta_turb * (HB6s_O2  - HB5_O2)
    HB6_N2  = HB5_N2  + eta_turb * (HB6s_N2  - HB5_N2)

    # And then idk mother nature will manage the rest


    # --- Work! Work! Work! ---
    m_H2O = n_H2O * 0.018015      # kg/s)
    m_O2_excess = n_O2_excess * M_O2
    m_N2_out = n_N2 * M_N2
    m_air = m_O2_supplied + m_N2

    #print(f"HB4_air {HB4_air}")
    #print(f"HB3_air {HB3_air}")

    W_comp = m_air * (HB4_air - HB3_air)
    #print(f"Compression work: {W_comp:.2f}")
    W_turb = m_H2O*(HB5_H20 - HB6_H2O) + m_O2_excess*(HB5_O2 - HB6_O2) + m_N2_out*(HB5_N2-HB6_N2)
    #print(f"Turbine work out: {W_turb:.2f}")

    W_net = W_turb - W_comp

    return W_net



 # Testing Unit
if __name__ == "__main__":

    m_H2 = 1.0087472201630838

    # Run the Brayton cycle function
    Dog = brayton_cycle_h2(
        m_H2
    )

    print(Dog)