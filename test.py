import sys
print(sys.path)
from Coolprop.Plots import PropertyPlot
plot = PropertyPlot('HEOS::Water', 'TS')
plot.calc_isolines()
plot.show()