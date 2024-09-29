# Imports e Configurações ---------------------------------------------------------------------------------------------

from classeCarta import *

#Plotagem

j = sp.I
Z = 60 + 60*j
chart = CartaSmith(Z0 = 60, grid_color = '#bdbdbd', figureSize = 7)
chart.plotS((Z, True))