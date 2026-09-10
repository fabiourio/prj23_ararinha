'''
This script shows how one can generate an airfoil with Strechinsky's parameters
and then export the airfoil coordinates to a text file.
'''

# IMPORTS
import airfoil_mod as am
import matplotlib.pyplot as plt

# INPUTS

# Strechinsky's coefficients
rle = 0.016174085507412374
xtmax = 0.299251037175028
tmax = 0.12006008605082567
phi = 0.28079208414754636
theta = 0.12545631889286601
epsilon = -0.07568231343983031
xcmax = 0.435612821228977
cmax = 0.019980465726939244
ctmax = 0.018475754870048083
cte = -7.699598486022169e-05
tte = 0.0025357931536794526

# Output file
filename = 'airfoil.dat'

# EXECUTION
xx = am.cos_bunching(0.0, 1.0, 61)

airfoil = am.strefoil(xx, rle, xtmax, tmax, phi,
                      theta, epsilon, xcmax, cmax, ctmax,
                      cte, tte, plot=True)

am.export_airfoil(airfoil, filename=filename)