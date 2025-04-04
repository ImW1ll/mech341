# I'll basically write a function that takes the q_in values, and simply computes what the water mass flow needs to be and what the state of the liquid will be afterwards
# I'm trying to force sat vapour so that it flies off in the steam tower.

from CoolProp.CoolProp import PropsSI

def towerCooling(heat_in):
    # River state 
    TR = 15 + 273.15  # chilly day I guess
    PR = 101325 # 1 atm
    
    HR = PropsSI('H', 'T', TR, 'P', PR, 'Water')  # H of that river water
    
    QC = 1  # we want it to float in the air
    HC = PropsSI('H', 'Q', QC, 'P', PR, 'Water') 
    TC = PropsSI('T', 'Q', QC, 'P', PR, 'Water') 
    
    # Mass flow rate required
    m_flow = heat_in / (HC - HR)  # kg/s
    
    return m_flow, TC


if __name__ == "__main__":
    water, temp = towerCooling(10e6)
    print(f"Water flow is : {water:.2f}")
    print(f"Water temp is : {temp:.2f}")