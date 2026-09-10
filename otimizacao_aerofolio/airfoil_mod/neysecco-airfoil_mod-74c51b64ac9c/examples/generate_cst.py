'''
This script shows how one can generate an airfoil with CST parameters
and then export the airfoil coordinates to a text file.
'''

# IMPORTS
import airfoil_mod as am
import matplotlib.pyplot as plt

# INPUTS

# CST coefficients
Au = [0.12566921, 0.14573835, 0.15191353, 0.21353185, 0.17871782, 0.20866563]
Al = [-0.13318821, -0.11891598, -0.22044409, -0.12810375, -0.08111686,  0.05141437]
tte = 0.0

# Output file
filename = 'airfoil.dat'

# EXECUTION
xx = am.cos_bunching(0.0, 1.0, 61)

airfoil = am.cstfoil(Au, Al, xx, N1=0.5, N2=1.0, tte=tte, plot=True)

am.export_airfoil(airfoil, filename=filename)