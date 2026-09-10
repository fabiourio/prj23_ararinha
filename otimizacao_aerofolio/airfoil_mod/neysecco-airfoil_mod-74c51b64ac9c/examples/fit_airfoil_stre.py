'''
This script receives a set of airfoil coordinates following
Selig's format and finds the corresponding CST values.
'''

# IMPORTS
import airfoil_mod as am
import matplotlib.pyplot as plt

#========================================
# INPUTS

# File with airfoil coordinates
airfoil_file = 'airfoil_data/naca2412.dat'
#airfoil_file = 'airfoil_data/naca64a006.dat'

# Additional characteristics
force_sharp_te = False
force_symmetric = False

#========================================
# EXECUTION

xx,yy = am.read_coordinates(airfoil_file)
xl,yl,xu,yu = am.split_skins(xx,yy)

rle, xtmax, tmax, phi, theta, epsilon, xcmax, cmax, ctmax, cte, tte = am.fit_stre(xu,yu,xl,yl,
                                                                                  force_sharp_te=force_sharp_te,
                                                                                  force_symmetric=force_symmetric)

xx_sample = am.cos_bunching(0.0, 1.0, 121)
airfoil = am.strefoil(xx_sample, rle, xtmax, tmax, phi,
                      theta, epsilon, xcmax, cmax, ctmax,
                      cte, tte, plot=False)

xf = airfoil['x_coord']
yf = airfoil['y_coord']

fig = plt.figure()
plt.plot(xx,yy,'o',label='samples')
plt.plot(xf,yf,label='Strechinsky')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best')
plt.axis('equal')
plt.tight_layout()
plt.show()