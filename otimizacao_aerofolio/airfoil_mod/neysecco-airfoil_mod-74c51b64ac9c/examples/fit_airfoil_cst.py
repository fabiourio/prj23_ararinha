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
airfoil_file = 'airfoil_data/rae2822.dat'
#airfoil_file = 'airfoil_data/naca64a006.dat'

# Number of CST coefficients on each skin
NN = 6

# Additional characteristics
force_sharp_te = True
force_symmetric = False

#========================================
# EXECUTION

xx,yy = am.read_coordinates(airfoil_file)
xl,yl,xu,yu = am.split_skins(xx,yy)

Au, Al, tte = am.fit_cst(xu,yu,xl,yl,NN,
                         force_sharp_te=force_sharp_te,
                         force_symmetric=force_symmetric)

xx_sample = am.cos_bunching(0.0, 1.0, 121)
airfoil = am.cstfoil(Au, Al, xx_sample, N1=0.5, N2=1.0, tte=tte, plot=False)

xf = airfoil['x_coord']
yf = airfoil['y_coord']

fig = plt.figure()
plt.plot(xx,yy,'o',label='samples')
plt.plot(xf,yf,label='CST')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best')
plt.axis('equal')
plt.tight_layout()
plt.show()