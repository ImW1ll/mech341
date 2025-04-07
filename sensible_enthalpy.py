from CoolProp.CoolProp import PropsSI

def sensible_enthalpy_PT(P, T, fluid='Water', T_ref=(25 + 273.15), tol=0.05):
    P_ref = 101325

    # === Special case: Manual N2-O2 mixture ===
    if fluid == 'Nitrogen[0.79]&Oxygen[0.21]':
        try:

            mass_frac_O2 = 0.232
            mass_frac_N2 = 0.768
            P_O2 = 0.21 * P
            P_N2 = 0.79 * P

            # Current enthalpies
            h_O2 = PropsSI('H', 'P', P_O2, 'T', T, 'Oxygen')
            h_N2 = PropsSI('H', 'P', P_N2, 'T', T, 'Nitrogen')

            h_actual = mass_frac_O2 * h_O2 + mass_frac_N2 * h_N2

            # Reference enthalpies
            P_O2_ref = 0.21 * P_ref
            P_N2_ref = 0.79 * P_ref
            h_O2_ref = PropsSI('H', 'P', P_O2_ref, 'T', T_ref, 'Oxygen')
            h_N2_ref = PropsSI('H', 'P', P_N2_ref, 'T', T_ref, 'Nitrogen')

            h_ref = mass_frac_O2 * h_O2_ref + mass_frac_N2 * h_N2_ref

            return h_actual - h_ref

        except Exception as e:
            raise ValueError(f"Failed to compute sensible enthalpy for N2/O2 mixture: {e}")

    # === Combustion go boom ===
    if fluid == 'CombustionProducts':
        try:
            # Mole fractions
            X_H2O = 0.19004524886877827
            X_O2 = 0.09502262443438914
            X_N2 = 0.7149321266968326

            M_H2O = 0.018015
            M_O2 = 0.032
            M_N2 = 0.0280134

            # Partial P (pushing P?)
            P_H2O = X_H2O * P
            P_O2 = X_O2 * P
            P_N2 = X_N2 * P

            # Current enthalpies
            h_H2O = PropsSI('H', 'P', P_H2O, 'T', T, 'Water')
            h_O2 = PropsSI('H', 'P', P_O2, 'T', T, 'Oxygen')
            h_N2 = PropsSI('H', 'P', P_N2, 'T', T, 'Nitrogen')

            h_actual = (h_H2O * X_H2O * M_H2O + h_O2 * X_O2 * M_O2 + h_N2 * X_N2 * M_N2) / (X_H2O * M_H2O + X_O2 * M_O2 + X_N2 * M_N2)

            # Reference enthalpies
            P_H2O_ref = X_H2O * P_ref
            P_O2_ref = X_O2 * P_ref
            P_N2_ref = X_N2 * P_ref
            h_H2O_ref = PropsSI('H', 'P', P_H2O_ref, 'T', T_ref, 'Water')
            h_O2_ref = PropsSI('H', 'P', P_O2_ref, 'T', T_ref, 'Oxygen')
            h_N2_ref = PropsSI('H', 'P', P_N2_ref, 'T', T_ref, 'Nitrogen')

            h_ref = (h_H2O_ref * X_H2O * M_H2O + h_O2_ref * X_O2 * M_O2 + h_N2_ref * X_N2 * M_N2) / (X_H2O * M_H2O + X_O2 * M_O2 + X_N2 * M_N2)

            return h_actual - h_ref

        except Exception as e:
            raise ValueError(f"Failed to compute sensible enthalpy for combustion products: {e}")

    # === Normal Shit===
    h_ref = PropsSI('H', 'T', T_ref, 'P', P_ref, fluid)

    # Try saturation temperature
    try:
        T_sat = PropsSI('T', 'P', P, 'Q', 0, fluid)
    except:
        T_sat = None

    # Try getting current enthalpy
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
