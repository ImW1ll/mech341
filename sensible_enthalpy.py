from CoolProp.CoolProp import PropsSI

def sensible_enthalpy_PT(P, T, fluid='Water', T_ref=(25 + 273.15), tol=0.05):
    P_ref = 101325

    # Ref enthalpy given in the project description
    h_ref = PropsSI('H', 'T', T_ref, 'P', P_ref, fluid)

    try:
        T_sat = PropsSI('T', 'P', P, 'Q', 0, fluid)
    except:
        T_sat = None

    # Now, I use TP in the function, if its a mix, coolprops goes kaboom, this fixed that
    try:
        h_actual = PropsSI('H', 'P', P, 'T', T, fluid)
    except ValueError:
        h_actual = None

    if (T_sat is None) or (T < (T_sat - tol)):
        if h_actual is None:
            h_actual = PropsSI('H', 'P', P, 'T', T, fluid)
        return h_actual - h_ref

    if T > (T_sat + tol):
        if h_actual is None:
            h_actual = PropsSI('H', 'P', P, 'T', T, fluid)
        h_sat_liq = PropsSI('H', 'P', P, 'Q', 0, fluid)
        h_sat_vap = PropsSI('H', 'P', P, 'Q', 1, fluid)
        sensible_liquid = h_sat_liq - h_ref
        superheat_portion = h_actual - h_sat_vap
        return sensible_liquid + superheat_portion

    if T <= T_sat:
        h_sat_liq = PropsSI('H', 'P', P, 'Q', 0, fluid)
        return h_sat_liq - h_ref
    else:
        h_sat_vap = PropsSI('H', 'P', P, 'Q', 1, fluid)
        h_sat_liq = PropsSI('H', 'P', P, 'Q', 0, fluid)
        sensible_liquid = h_sat_liq - h_ref
        superheat_portion = 0.0
        return sensible_liquid + superheat_portion
