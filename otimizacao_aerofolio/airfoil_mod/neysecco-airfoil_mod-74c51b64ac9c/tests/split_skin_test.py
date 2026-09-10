'''
This script reads coordinates from a dat file taken from AirfoilTools website,
and then splits the lower and upper skin coordinates
'''

# IMPORTS
import airfoil_mod as am
import matplotlib.pyplot as plt

# INPUTS
file_path = 'airfoil_data/naca23012.dat'

# EXECUTION

xx,yy = am.read_coordinates(file_path)

xl,yl,xu,yu = am.split_skins(xx, yy)

fig = plt.figure()
plt.plot(xx,yy,'o',label='reference')
plt.plot(xl,yl,label='lower skin')
plt.plot(xu,yu,label='upper skin')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best')
plt.axis('equal')
plt.tight_layout()
plt.show()