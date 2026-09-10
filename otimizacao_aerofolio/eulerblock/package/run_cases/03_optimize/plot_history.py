'''
This script can be used to monitor the optimization convergence.

While executing the optimization, run this script in a separate
terminal to generate the optimization history.

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

import pickle
import matplotlib.pyplot as plt
import numpy as np

with open('opt_results.pickle','rb') as fid:
    hist_dict = pickle.load(fid)
    xxhist = hist_dict['xxhist']
    CDhist = hist_dict['CDhist']
    CLhist = hist_dict['CLhist']
    maxthist = hist_dict['maxthist']
    minthist = hist_dict['minthist']
    mach = hist_dict['mach']
    Al = hist_dict['Al']
    Au = hist_dict['Au']

xxhist = np.array(xxhist)
CDhist = np.array(CDhist)
CLhist = np.array(CLhist)
maxthist = np.array(maxthist)
minthist = np.array(minthist)

fig = plt.figure()
plt.subplot(411)
plt.plot(xxhist)
plt.ylabel(r'$x$')
plt.subplot(412)
plt.plot(CDhist)
plt.ylabel(r'$C_D$')
plt.subplot(413)
plt.plot(CLhist)
plt.ylabel(r'$C_L$')
plt.subplot(414)
plt.plot(maxthist)
plt.plot(minthist)
plt.ylabel(r'$t/c$')
plt.xlabel('Function evaluations')
plt.tight_layout()

# The number of airfoil coefficients is half of xx
# after discounting the alpha
Nvar = (len(xxhist[0,:])-1)//2

Al = xxhist[-1,:Nvar]
Au = xxhist[-1,Nvar:2*Nvar]
alpha = xxhist[-1,-1]
xa = (1-np.cos(np.linspace(0, 1, 81)*np.pi))/2

airfoil = am.cstfoil(Au, Al, xa)
x = airfoil['x_coord']
y = airfoil['y_coord']
maxt = airfoil['max_thickness']
xmaxt = airfoil['x_max_thickness']
mint = airfoil['min_thickness']
maxc = airfoil['max_camber']
xmaxc = airfoil['x_max_camber']

# Export airfoil coordinates
am.export_airfoil(airfoil)

# Read the most recent wall.dat file stored in the folder
CL,CD,CM,xxD,Cp,MachD = eb.postprocess(x,y,alpha,mach,plot=True)

print('Results for the last airfoil in history')
print('max thickness=',maxt)
print('max thickness x/c=',xmaxt)
print('max camber=',maxc)
print('max camber x/c=',xmaxc)
print('CL=',CL)
print('CD=',CD)
print('CM=',CM)

plt.show()
