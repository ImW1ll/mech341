extra = 15.160e6 # j/kg of alu
m_alumina = 17.02 # kg/s
m_h2 = 1.01 # kg/s
m_al = 9

# We know temp of water and alu is at 0, so now we need to remove the enthalpy to heat up H2 and Alumina to 500
from metalAndShit import *
from CoolProp.CoolProp import PropsSI

highAlumina = enthalpy_alumina_highT(500+273.15, 3e7)
lowAlumina = enthalpy_alumina_lowT(25+273.15, 3e7)
delAlumina = (highAlumina - lowAlumina) * m_alumina # energy to heat up from stap to press

# Now H2
highHydrogen = PropsSI('H', 'T', (500+273.15), 'P', 3e7, 'Hydrogen')
lowHydrogen = PropsSI('H', 'T', (25+273.15), 'P', 3e7, 'Hydrogen')
delHydrogen = ( highHydrogen - lowHydrogen) * m_h2

print(((extra*9) - delAlumina - delHydrogen) /(m_al*1e6))