# =============== BRAYTON CYCLE FUNCTION ======================
def brayton_cycle_h2(
    m_H2,            # [kg/s] mass flow of hydrogen
    T1,              # [K] inlet temperature to compressor
    P1,              # [Pa] inlet pressure to compressor (1 bar = 1e5 Pa)
    pr,              # [-] overall pressure ratio (compressor outlet / inlet)
    cp_air,          # [J/(kg·K)] specific heat of air
    gamma_air,       # [-] ratio of specific heats for air
    LHV_H2,          # [J/kg] lower heating value of hydrogen
    m_air,           # [kg/s] mass flow of air
    eta_comp=1.0,    # [-] isentropic efficiency of compressor (ideal=1.0)
    eta_turb=1.0     # [-] isentropic efficiency of turbine (ideal=1.0)
):
    """
    Returns:
      W_net : Net shaft power from Brayton cycle [W]
      T2    : Compressor outlet temperature [K]
      T3    : Turbine inlet temperature (post-combustion) [K]
      T4    : Turbine outlet temperature [K]
    """
    # ------------------------------------------------------------------------
    # 1) Compressor
    # ------------------------------------------------------------------------
    P2 = pr * P1
    T2_ideal = T1 * (pr)**((gamma_air - 1.0) / gamma_air)
    T2 = T1 + (T2_ideal - T1)/eta_comp  # corrected for compressor efficiency
    W_comp = m_air * cp_air * (T2 - T1) # compressor work (W)

    # ------------------------------------------------------------------------
    # 2) Combustion
    # ------------------------------------------------------------------------
    Q_dot_in = m_H2 * LHV_H2  # chemical energy input from H2 (J/s)
    m_total = m_air + m_H2
    # approximate cp of mixture as cp_air (idealization)
    cp_mix = cp_air
    T3 = T2 + Q_dot_in / (m_total * cp_mix)

    # ------------------------------------------------------------------------
    # 3) Turbine expansion (P2 -> P1)
    # ------------------------------------------------------------------------
    expansion_ratio = P1 / P2
    T4_ideal = T3 * (expansion_ratio)**((gamma_air - 1.0) / gamma_air)
    T4 = T3 - eta_turb * (T3 - T4_ideal)
    W_turb = m_total * cp_mix * (T3 - T4)  # turbine work (W)

    # ------------------------------------------------------------------------
    # 4) Net Brayton cycle work
    # ------------------------------------------------------------------------
    W_net = W_turb - W_comp
    return W_net, T2, T3, T4