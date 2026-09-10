'''
This script is an example of how to use the EulerBlock interface
to analyze an airfoil.

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2022
'''

# Modify path to include the Euler solver
import sys
with open('../eulerblock_path.txt') as f:
    exec(f.read())
sys.path.append(eulerblock_path)
from eulerblock import euler_mod as eb
import airfoil_mod as am

import numpy as np
import matplotlib.pyplot as plt

# Flow conditions
alpha = 1.405*np.pi/180.0 # Angle of attack [rad]
mach = 0.75 # mach number
gamma = 1.4 # ratio of specific heats of the fluid

# Solver parameters
order = 2 # 1 for 1st; 2 for 2nd
iter = 20000 # Use 50000 for first order
dt = 0.001 # Global time step (only used if use_local_dt==0) .Use 0.002 for first order and 0.001 for second order.
CFL = 0.2 # Courant number used to compute local time steps (only used if use_local_dt==1)
use_local_dt = 1 # Time stepping-strategy. 0 for global dt; 1 for local dt based on CFL.
res_NK = 1e-4
res_tol = 1e-6
reinitialize = 0
adj_funcs = ['cl_jlow', 'cd_jlow', 'cm_jlow']

# Discretization (higher grid levels have higher number of elements)
grid_level = 1.0

# Reference mesh properties
Nchord0 = 30 # Number of elements along the chord
NJ0 = 48 # Number of elements along the normal direction
s00 = 0.5e-2 # Reference mesh spacing

Nchord = int((Nchord0*grid_level) + 1)
NJ = int((NJ0*grid_level) + 1)
s0 = NJ0/(NJ-1)*s00

# Airfoil CST values - ClarkY
#Al = [-0.1294, -0.0036, -0.0666]
#Au = [0.206, 0.2728, 0.2292]

# Airfoil CST values - NACA1411
Al = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]
Au = [ 0.16146332,  0.18349204,  0.14126241,  0.18194397]

# Airfoil CST values - NACA0012
#Al = [-0.16941, -0.15145, -0.13904, -0.13988]
#Au = [ 0.16941,  0.15145,  0.13904,  0.13988]

# Airfoil CST values - RAE2822 (remember to lower alpha)
#Au = [0.12569744, 0.14564925, 0.1520378,  0.21346977, 0.17867926, 0.20872853]
#Al = [-0.13334359, -0.11827498, -0.22179445, -0.12644608, -0.08233977,  0.05193736]

results = eb.run_cst(Al, Au, Nchord, alpha, mach,
                     gamma, order,
                     iter, dt, CFL, use_local_dt, res_NK, res_tol,
                     reinitialize, adj_funcs, plot=True, NJ=NJ, s0=s0)

CL = results['CL']
CD = results['CD']
CM = results['CM']
maxt = results['maxt']
xmaxt = results['xmaxt']
mint = results['mint']
maxc = results['maxc']
xmaxc = results['xmaxc']
xx = results['distrib']['xx']
Cp = results['distrib']['Cp']
MachD = results['distrib']['Mach']
grads = results['grads']

print('Results')
print('max thickness=',maxt)
print('max thickness x/c=',xmaxt)
print('max camber=',maxc)
print('max camber x/c=',xmaxc)
print('CL=',CL)
print('CD=',CD)
print('CM=',CM)
print('grads=',grads)

# Refine coordinates for xfoil export
x = (1-np.cos(np.linspace(0, 1, 81)*np.pi))/2
airfoil = am.cstfoil(Au, Al, x)
xf = airfoil['x_coord']
yf = airfoil['y_coord']
maxt = airfoil['max_thickness']
xmaxt = airfoil['x_max_thickness']
mint = airfoil['min_thickness']
maxc = airfoil['max_camber']
xmaxc = airfoil['x_max_camber']
am.export_airfoil(airfoil, title='AIRFOIL', filename='airfoil.dat')

# Example of how to plot Cp curve and airfoil coordinates
fig = plt.figure()
plt.subplot(311)
plt.plot(xx,Cp)
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.subplot(312)
plt.plot(xx,MachD)
plt.ylabel('Mach')
plt.subplot(313)
plt.plot(xf,yf)
plt.axis('equal')
plt.ylabel('y/c')
plt.xlabel('x/c')
plt.tight_layout()
plt.show()