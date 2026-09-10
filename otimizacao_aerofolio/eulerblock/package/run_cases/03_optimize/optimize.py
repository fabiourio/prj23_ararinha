'''
This script sets up an optimization problem to design airfoils
analyzed with Euler Equations.

The optimization problem definition is:
min CD
w.r.t alpha, CST coefficients
s.t. t >= tref
     CL = CLref

The user may modify the parameters under the "INPUTS" section
of this script.

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2024
'''

# Modify path to include the Euler solver
import sys
with open('../eulerblock_path.txt') as f:
    exec(f.read())
sys.path.append(eulerblock_path)
from eulerblock import euler_mod as eb

import numpy as np
from scipy.optimize import minimize, Bounds
import pickle
import time

#===============================================
# INPUTS

# Flow conditions
mach = 0.75
CLref = 0.5

# Minimum thickness allowed
tref = 0.11
tlower = 0.01

# Euler solver parameters
order = 2 # 1 for 1st; 2 for 2nd
iter = 20000 # max number of iterations (Use 50000 for first order)
dt = 0.001 # Global time step (only used if use_local_dt==0). Use 0.002 for first order and 0.001 for second order. Reduce if you are getting NaNs or unstable solutions
CFL = 0.2 # Courant number used to compute local time steps (only used if use_local_dt==1). Reduce if you are getting NaNs or unstable solutions
use_local_dt = 1 # Time stepping-strategy. 0 for global dt; 1 for local dt based on CFL.
res_NK = 1e-5 # Tolerance to switch from Runge-Kutta to Newton-Krylov solver. This value should be greater than res_tol.
res_tol = 1e-8 # Tolerance on residual MSE to stop iterations. I do not recommend using values greater than 1e-6 (e.g. 1e-5), otherwise the numerical noise will break finite-difference estimations

# Guess for initial angle of attack (this will be changed during the optimization)
#alpha = 2.06*np.pi/180.0 # CL=0.5 for first order
alpha = 1.405*np.pi/180.0 # CL=0.5 for second order

# Airfoil discretization
Nchord = 31
cref = 1.0 # reference chord

# Reference airfoil CST variables - NACA 1411
Al = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]
Au = [ 0.16146332,  0.18349204,  0.14126241,  0.18194397]

# Bounds for design variables
Al_lower = [-1.00, -1.00, -1.00, -1.00]
Al_upper = [-0.05, -0.05, -0.05, -0.05]

Au_lower = [ 0.05,  0.05,  0.05,  0.05]
Au_upper = [ 1.00,  1.00,  1.00,  1.00]

# Number of airfoil CST parameters (design variables of the problem)
Nvar = len(Al) # This value is used for the upper and lower skin separately (the total number of DVs is 2*Nvar)

# Angle of attack bounds during optimizations
alpha_min = -3*np.pi/180
alpha_max =  7*np.pi/180

# Select if we should use reinitialization for the first CFD analysis.
# reinitialize = [0] does not use reinitialization for the first CFD analysis.
# reinitialize = [1] uses reinitialization for the first CFD analysis.
# We set it as a list so it remains accessible within the the check_point function.
reinitialize = [0]

# Option to reuse results from previous optimization.
# This may be very useful to continue a previous optimization that was interrupted.
# Make a copy of the opt_result.pickle file from a previous optimization, then add its name
# to the variable below.
# Use None to avoid this reinitialization.
reload_opt = None
#reload_opt = 'opt_results_copy.pickle'

#===============================================
# EXECUTION

# Initial airfoil design variables
# The structure of the DV array is:
# xx = [lower skin CST coefficients (Al); upper skin CST coefficients (Au); alpha]
xx0 = np.hstack([Al, Au, alpha])

# Initialize optimization history
if reload_opt is None:
    xxhist = []
    CDhist = []
    CLhist = []
    maxthist = []
    minthist = []
    gradshist = []

else:

    with open(reload_opt,'rb') as fid:
        hist_dict = pickle.load(fid)
        xxhist = hist_dict['xxhist']
        CDhist = hist_dict['CDhist']
        CLhist = hist_dict['CLhist']
        maxthist = hist_dict['maxthist']
        minthist = hist_dict['minthist']
        gradshist = hist_dict['gradshist']

# Define function to check if the desired point was already executed
def check_point(xx, tol=1e-11, force_run=False):

    # Loop over the history to find previous design point within tolerance
    ii_closest = -1
    for ii,xx_curr in enumerate(xxhist):
        delta = np.sqrt(np.sum((xx_curr - xx)**2))

        if delta < tol:
            print('Found previous solution')
            print('ii:',ii)
            print('delta:',delta)
            ii_closest = ii
            break

    # Get design variables
    Al = xx[:Nvar]
    Au = xx[Nvar:2*Nvar]
    alpha = xx[-1]
    
    print('Evaluating the following point')
    print('Al:',Al)
    print('Au:',Au)
    print('alpha [deg]:',alpha*180/np.pi)

    # Choose whether to run new case or use history
    if (ii_closest == -1) or force_run: # i_closest keeps the default value if no match is found

        # Call Euler solver with adjoint for each function
        results = eb.run_cst(Al, Au, Nchord, alpha, mach,
                             gamma=1.4, order=order,
                             iter=iter, dt=dt, CFL=CFL, use_local_dt=use_local_dt, res_NK=res_NK,
                             res_tol=res_tol, reinitialize=reinitialize[0], plot=False,
                             adj_funcs=['cl_jlow','cd_jlow'])
        CL = results['CL']
        CD = results['CD']
        maxt = results['maxt']
        mint = results['mint']
        grads1 = results['grads']

        # Gather gradients w.r.t. xx
        grads = {}
        grads['CL'] = np.hstack([grads1['cl_jlow']['dAl'],
                                 grads1['cl_jlow']['dAu'],
                                 [grads1['cl_jlow']['alpha']]])

        grads['CD'] = np.hstack([grads1['cd_jlow']['dAl'],
                                 grads1['cd_jlow']['dAu'],
                                 [grads1['cd_jlow']['alpha']]])

        grads['maxt'] = np.hstack([grads1['maxt']['dAl'],
                                   grads1['maxt']['dAu'],
                                   [0.0]])
        grads['mint'] = np.hstack([grads1['mint']['dAl'],
                                   grads1['mint']['dAu'],
                                   [0.0]])

        xxhist.append(xx.copy())
        CDhist.append(CD)
        CLhist.append(CL)
        maxthist.append(maxt)
        minthist.append(mint)
        gradshist.append(grads)
        with open('opt_results.pickle','wb') as fid:
            pickle.dump({'xxhist':xxhist,
                         'CDhist':CDhist,
                         'CLhist':CLhist,
                         'maxthist':maxthist,
                         'minthist':minthist,
                         'mach':mach,
                         'gradshist':gradshist,
                         'Al':Al,
                         'Au':Au}, fid)


    else:
        
        # Gather values from previous solution
        CL = CLhist[ii_closest]
        CD = CDhist[ii_closest]
        maxt = maxthist[ii_closest]
        mint = minthist[ii_closest]
        grads = gradshist[ii_closest]

    # Print and store history
    print('evaluation:',len(CDhist))
    print('CL:',CL)
    print('CD:',CD)
    print('maxt:',maxt)
    print('mint:',mint)
    print('Al:',Al)
    print('Au:',Au)
    print('alpha [deg]:',alpha*180/np.pi)
    print('grads:',grads)
    print('')

    # Check if we need to avoid solution reinitialization if NaNs show up
    if np.isnan(CL) or np.isnan(CD):
        reinitialize[0] = 0
    else:
        reinitialize[0] = 1

    # Return design functions
    return CL, CD, maxt, mint, grads


# Define objective function
def objfun(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    return CD

def objfungrad(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    return grads['CD']

# Define constraints
def ineqconfun(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    g1 = (maxt - tref)
    g2 = (tlower - mint)

    return g1, g2

def ineqconfungrad(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    return grads['maxt'], -grads['mint']

# Define constraints
def eqconfun(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    h1 = (CL - CLref)

    return h1

def eqconfungrad(xx):

    CL, CD, maxt, mint, grads = check_point(xx)

    return grads['CL']

# Run initial case
objfun(xx0)
reinitialize[0] = 1

# Create list of constraints
con1 = {'type': 'ineq',
        'fun': ineqconfun,
        'jac': ineqconfungrad}
con2 = {'type': 'eq',
        'fun': eqconfun,
        'jac': eqconfungrad}
cons = [con1, con2]

# Set DV bounds
#lb = [-Amax]*Nvar + [Amin]*Nvar + [alpha_min]
#ub = [-Amin]*Nvar + [Amax]*Nvar + [alpha_max]
lb = Al_lower + Au_lower + [alpha_min]
ub = Al_upper + Au_upper + [alpha_max]
bounds = Bounds(lb, ub, keep_feasible=True)

# Set optimizer options.
# Adjusted eps (finite difference size) according to expected DV bounds.
options={'maxiter': 100, 'ftol': 1e-05, 'iprint': 1, 'disp': True}#, 'finite_diff_rel_step': None}

# Run optimizer
start_time = time.time()
result = minimize(objfun, x0=xx0, jac=objfungrad,
                  constraints=cons, bounds=bounds,
                  method='slsqp',
                  options=options)
end_time = time.time()

# Print results
print('')
print('=========================================')
print('OPTIMIZATION RESULTS')
print(result)
print('Optimization time: %d seconds'%(end_time-start_time))
print('=========================================')
print('')
xopt = result.x

# Run optimum point to make sure it is the last item
# of the history file
CL, CD, maxt, mint, grads = check_point(xopt, force_run=True)

print('')
print('=========================================')
print('OPTIMUM POINT')
print('CL:',CL)
print('CD:',CD)
print('maxt:',maxt)
print('xopt:',xopt)
print('=========================================')
print('')
