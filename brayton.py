import math

# 1) Polynomials for cp in kJ/(kmol.K) and the converter to J/(kg.K)
def cp_air_kj_per_kmolK(T):
    """
    Returns cp of air in kJ/(kmol.K) at temperature T [K].
    Example polynomial constants below.
    """
    a = 28.11
    b = 0.1967e-2
    c = 0.4802e-5
    d = -1.966e-9   
    return a + b*T + c*(T**2) + d*(T**3)

def cp_air_j_per_kgK(T):
    """
    Returns cp of air in J/(kg.K) at temperature T [K],
    using the polynomial above.
    """
    # Molar mass of air [kg/kmol]
    M_air = 28.97
    cp_kj_per_kmolK = cp_air_kj_per_kmolK(T)
    cp_j_per_kgK = (cp_kj_per_kmolK * 1000.0) / M_air
    return cp_j_per_kgK

def brayton_cycle_h2(
    m_H2,            # [kg/s] mass flow of hydrogen
    T1,              # [K] inlet temperature to compressor
    P1,              # [Pa] inlet pressure to compressor (1 bar = 1e5 Pa)
    pr,              # [-] overall pressure ratio (compressor outlet / inlet)
    gamma_guess,     # [-] an approximate gamma (used for the isentropic T2_ideal)
    LHV_H2,          # [J/kg] lower heating value of hydrogen
    m_air,           # [kg/s] mass flow of air
    eta_comp=1.0,    # [-] isentropic efficiency of compressor (ideal=1.0)
    eta_turb=1.0     # [-] isentropic efficiency of turbine (ideal=1.0)
):
    """
    Returns:
      W_net    : Net shaft power [W]
      T2       : Compressor outlet temperature [K]
      T3       : Turbine inlet temperature (post-combustion) [K]
      T4       : Turbine outlet temperature [K]

    This version uses temperature-dependent cp (from the polynomial) to compute
    compressor and turbine work more accurately. However, it still uses a
    simple 'textbook' approach to isentropic temperature ratios with an
    approximate gamma_guess. That means T2_ideal and T4_ideal come from
    T2_ideal = T1*(pr)^((gamma_guess-1)/gamma_guess), etc.
    Then we use the average of cp(T1) and cp(T2) to estimate W_comp, etc.
    """
    # ------------------------------------------------------------------------
    # 1) Compressor: isentropic temperature T2_ideal from T1, then actual T2
    # ------------------------------------------------------------------------
    P2 = pr * P1

    # Ideal isentropic temperature rise (still using a constant gamma approximation):
    T2_ideal = T1 * (pr)**((gamma_guess - 1.0) / gamma_guess)

    # Actual T2 factoring in isentropic efficiency
    #   (T2 - T1) = (T2_ideal - T1)/eta_comp
    T2 = T1 + (T2_ideal - T1)/eta_comp

    # Approximate the compressor work via average cp over [T1, T2]
    cp_comp_in  = cp_air_j_per_kgK(T1)
    cp_comp_out = cp_air_j_per_kgK(T2)
    cp_comp_avg = 0.5*(cp_comp_in + cp_comp_out)
    W_comp = m_air * cp_comp_avg * (T2 - T1)

    # ------------------------------------------------------------------------
    # 2) Combustor: add heat from H2, find T3
    #    Q_dot_in = m_H2 * LHV_H2  [J/s]
    #    (m_air + m_H2)*cp_avg*(T3 - T2) = Q_dot_in
    #
    # For a first approximation, we'll evaluate cp at some average temperature
    # guess or just at T2.  A more rigorous approach would iterate or do a
    # polynomial enthalpy solve. We'll keep it simple here.
    # ------------------------------------------------------------------------
    Q_dot_in = m_H2 * LHV_H2
    m_total = m_air + m_H2

    # cp for the mixture? For simplicity, just use air's cp at T2 or an average
    cp_combust_in  = cp_air_j_per_kgK(T2)  # could do something more elaborate
    # We'll guess T3 by ignoring the T-dependence in the combustor, except at T2
    T3 = T2 + Q_dot_in/( m_total * cp_combust_in )

    # ------------------------------------------------------------------------
    # 3) Turbine expansion from (T3, P2) down to P1
    #    Again we do an isentropic T4_ideal from the same gamma_guess,
    #    then apply turbine efficiency, and use an average cp for the enthalpy.
    # ------------------------------------------------------------------------
    expansion_ratio = P1 / P2
    T4_ideal = T3 * (expansion_ratio)**((gamma_guess - 1.0) / gamma_guess)
    T4 = T3 - eta_turb*(T3 - T4_ideal)

    # Approx turbine work with average cp over [T3, T4]
    cp_turb_in  = cp_air_j_per_kgK(T3)
    cp_turb_out = cp_air_j_per_kgK(T4)
    cp_turb_avg = 0.5*(cp_turb_in + cp_turb_out)
    W_turb = m_total * cp_turb_avg * (T3 - T4)

    # ------------------------------------------------------------------------
    # 4) Net work [W]
    # ------------------------------------------------------------------------
    W_net = W_turb - W_comp

    return W_net, T2, T3, T4

if __name__ == "__main__":
    # ------------------------------
    # Example usage
    # ------------------------------

    # Mass flow of H2 [kg/s]
    #m_H2 = 0.05  # 0.05 kg/s of hydrogen
    m_H2 = 0.1117 * 8.75

    # Inlet conditions
    T_in_C = 200.58          # degC
    T1 = T_in_C + 273.15     # K
    P1_bar = 1.0
    P1 = P1_bar * 1e5        # Pa

    # Assumed Brayton cycle parameters
    pr          = 10.0        # Overall pressure ratio
    gamma_guess = 1.4         # Just an approximate gamma for the isentropic formula
    LHV_H2      = 120e6       # J/kg (120 MJ/kg for H2)
    m_air       = 60         # [kg/s] air flow
    eta_comp    = 0.88        # compressor isentropic efficiency
    eta_turb    = 0.90        # turbine isentropic efficiency

    # Run the Brayton cycle function
    W_net, T2, T3, T4 = brayton_cycle_h2(
        m_H2, T1, P1, pr,
        gamma_guess,
        LHV_H2, m_air,
        eta_comp=eta_comp,
        eta_turb=eta_turb
    )

    # Print results
    print(f"Compressor outlet temperature T2 = {T2:.2f} K")
    print(f"Turbine inlet temperature    T3 = {T3:.2f} K")
    print(f"Turbine outlet temperature   T4 = {T4:.2f} K")
    print(f"Net shaft power (approx)     W_net = {(W_net / 1e6 ):.2f} MW")
