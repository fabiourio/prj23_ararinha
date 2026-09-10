'''
This script shows how one can generate an airfoil with NACA parameters
and then export the airfoil coordinates to a text file.
'''

# IMPORTS
import airfoil_mod as am
import matplotlib.pyplot as plt

# INPUTS

# 4- or 5-digit NACA
naca_string = '23012'

# Output file
filename = 'airfoil.dat'

# EXECUTION
xx = am.cos_bunching(0.0, 1.0, 61)

airfoil = am.nacafoil(naca_string, xx, plot=True)

am.export_airfoil(airfoil, filename=filename)