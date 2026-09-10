'''
This module is the main interface to call EulerBlock
from Python.

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2022
Rev: Sept 2024 - added adjoint method
'''

from . import gridgen2d
from . import gridwarp
import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import airfoil_mod as am

#=======================================

def run_case(X, Y, alpha, mach, bcs, xref=0.25, yref=0.0, cref=1.0,
             gamma=1.4, order=2,
             iter=10000, dt=0.001, CFL=0.2, use_local_dt=1, res_NK=1e-4, res_tol=1e-5,
             reinitialize=0, adj_funcs=[],
             plot=False, use_f2py=False):

    '''
    This script takes a general grid defined by matrices X and Y and executes the
    Euler solver.

    bcs: list -> 4 integers indicating boundary conditions at ilow,
                 jlow, ihigh, and jhigh boundaries (1:farfield, 2:wall, 3:periodic)
    '''

    #-------------------------------
    # WRITING GRID FILE

    # Plot the grid if the user wants to visulize it
    if plot:
        fig = gridgen2d.plot_grid(X,Y,show=True)
        fig.savefig('grid.png',dpi=300)

    # Save the grid file that will be read by the Fortran code
    gridgen2d.export_plot3d(X,Y,'grid.xyz',1,2)
    
    #-------------------------------
    # FLOW SOLUTION

    if use_f2py:

        from .euler_f2py import run_case_f2py

        wallData, grads = run_case_f2py(X, Y, alpha, mach, bcs, xref, yref, cref,
                                        gamma, order,
                                        iter, dt, CFL, use_local_dt, res_NK, res_tol, reinitialize,
                                        plot, adjoint_funcs=adj_funcs)

    else:

        # Write input file
        write_settings(alpha, mach, xref, yref, cref, gamma, order,
                       iter, dt, CFL, use_local_dt, res_NK, res_tol, bcs, reinitialize, adj_funcs)

        # Call the Fortran code
        subprocess.run(os.path.join(os.path.dirname(os.path.realpath(__file__)),'eulerblock.exe'))

        # Post-processing
        wallData = read_walls()

        # Adjoint not yet implemented in Fortran
        grads = read_derivatives(adj_funcs)

    return wallData, grads

#=======================================

def run_airfoil(x, y, alpha, mach, xref=0.25, yref=0.0, cref=1.0,
                gamma=1.4, order=2,
                iter=10000, dt=0.001, CFL=0.2, use_local_dt=1, res_NK=1e-4, res_tol=1e-5,
                reinitialize=0, adj_funcs=[],
                plot=False, NJ=49, s0=0.5e-2, use_f2py=False,
                dx=None, dy=None):
    '''
    This script generates an O-mesh around an airfoil and calls
    the EulerBlock block to solve the flow around it.
    The x and y coordinates should start and finish at the trailing edge,
    looping first around the lower skin.
    The trailing edge point should NOT be repeated.

    gamma: float -> Ratio of specific heats of the fluid
    order: integer -> 1 for 1st order; 2 for 2nd order solution
    iter: integer -> maximum number of timestep iterations
    dt: float -> timestep value. Decrease this value if you are getting
                 NaNs during the iterations. This is only used if use_local_dt==0.
    CFL: float -> Courant number used to compute dt for local time-stepping.
                  This is only used if use_local_dt==1.
    use_local_dt: integer -> 0 for global dt; 1 for local dt based on CFL
    res_NK: float -> Residual threshold to change from RK5 to NK
    res_tol: float -> Residual MSE value to stop iterations
    reinitialize: logical -> Reloads the flowfield stored in solution.vtk file
    plot: logical -> Draw mesh and solution plots
    NJ: integer -> Number of mesh layers to be extruded to the farfield
    s0: float -> Initial layer size
    dx and dy: arrays if the same size as x and y representing surface
               perturbations for mesh warping. We use mesh warping because
               it is easier to compute derivatives.
    '''
    
    # Join x and y coordinates in a single array.
    # This is how gridgen2d reads the airfoil mesh.
    # We also ensure that the airfoil is closed (first and last point matches).
    if (x[0]==x[-1]) and (y[0]==y[-1]):
        r0 = np.vstack([x,
                        y]).T.ravel()
    else:
        r0 = np.vstack([np.hstack([x,x[0]]),
                        np.hstack([y,y[0]])]).T.ravel()
    
    bc1 = 'continuous'
    bc2 = 'continuous'

    # Non-dimensionalized mesh extension (farfield position)
    extension = 40

    options = {
        'epsE0':2.0, # Dissipation parameter (change this if mesh is crossing itself)
    }

    # Generate the mesh
    X,Y = gridgen2d.march_layer(r0,bc1,bc2,s0,NJ,extension,**options)

    # Transpose to follow the right-hand rule of i and j coordinates
    X = X.T
    Y = Y.T

    if (dx is not None) and (dy is not None):
        X, Y = gridwarp.warp_mesh(x, y, dx, dy, X, Y)

    # Boundary conditions
    bcs = [3, 2, 3, 1]

    # Call the Euler Solver
    wallData, grads = run_case(X, Y, alpha, mach, bcs, xref, yref, cref,
                        gamma, order,
                        iter, dt, CFL, use_local_dt, res_NK, res_tol, reinitialize, adj_funcs,
                        plot, use_f2py)

    # Gather outputs
    if dx is not None:
        CL, CD, CM, xx, Cp, Mach = postprocess(x+dx, y+dy, alpha, mach, wallData, plot=plot)
    else:
        CL, CD, CM, xx, Cp, Mach = postprocess(x, y, alpha, mach, wallData, plot=plot)

    results = {'alpha':alpha,
               'Mach':Mach,
               'CL':CL,
               'CD':CD,
               'CM':CM,
               'distrib':{'xx':xx,
                          'Cp':Cp,
                          'Mach':Mach},
               'grads':grads,
               'X':X,
               'Y':Y,
               'x':x,
               'y':y}

    return results

#=======================================

def run_naca(naca_string, Nchord, alpha, mach,
             gamma=1.4, order=2,
             iter=10000, dt=0.001, CFL=0.2, use_local_dt=1, res_NK=1e-4, res_tol=1e-5,
             reinitialize=0, adj_funcs=[],
             plot=False, NJ=49, s0=0.5e-2, use_f2py=False):
    '''
    This is a specialized version of the run_case function
    that employs the 4 or 5-digit NACA parametrization over a unit chord
    airfoil.

    Surface mesh desrivatives are not implemented yet.

    Nchord: number of points along the upper/lower skin of the
    airfoil (the total number of points will be 2*Nchord-2)
    '''

    xref = 0.25
    yref = 0.0
    cref = 1.0

    x = (1-np.cos(np.linspace(0, 1, Nchord)*np.pi))/2

    airfoil = am.nacafoil(naca_string, x)

    xf = airfoil['x_coord']
    yf = airfoil['y_coord']
    maxt = airfoil['max_thickness']
    xmaxt = airfoil['x_max_thickness']
    mint = airfoil['min_thickness']
    maxc = airfoil['max_camber']
    xmaxc = airfoil['x_max_camber']

    results = run_airfoil(xf, yf, alpha, mach, xref, yref, cref,
                          gamma, order, iter, dt, CFL, use_local_dt, res_NK, res_tol,
                          reinitialize, adj_funcs, plot, NJ, s0, use_f2py)
    
    results['maxt'] = maxt
    results['xmaxt'] = xmaxt
    results['mint'] = mint
    results['maxc'] = maxc
    results['xmaxc'] = xmaxc

    return results

#=======================================

def run_cst(Al, Au, Nchord, alpha, mach,
            gamma=1.4, order=2,
            iter=10000, dt=0.001, CFL=0.2, use_local_dt=1, res_NK=1e-4, res_tol=1e-5,
            reinitialize=0, adj_funcs=[],
            plot=False, NJ=49, s0=0.5e-2, use_f2py=False):
    '''
    This is a specialized version of the run_case function
    that employs the CST parametrization over a unit chord
    airfoil.
    This function is the same as run_cst_warp, but here we
    compute derivatives directly over the Al and Au coefficients
    rather than their deformations (dAl, dAu).
    The baseline mesh is fixed and based on an
    airfoil of constant coefficients of 0.1.
    '''

    # Make sure inputs are arrays
    Al = np.array(Al)
    Au = np.array(Au)

    # Reference airfoil CST variables
    Alref = Al.copy() #np.ones_like(Al)*(-0.1)
    Auref = Au.copy() #np.ones_like(Au)*0.1

    results = run_cst_warp(Alref, Auref, Nchord, alpha, mach,
                           gamma=gamma, order=order,
                           iter=iter, dt=dt, CFL=CFL, use_local_dt=use_local_dt, res_NK=res_NK, res_tol=res_tol,
                           reinitialize=reinitialize, adj_funcs=adj_funcs,
                           plot=plot, NJ=NJ, s0=s0, use_f2py=use_f2py,
                           dAl=Al-Alref, dAu=Au-Auref)

    return results

#=======================================

def run_cst_warp(Al, Au, Nchord, alpha, mach,
                 gamma=1.4, order=2,
                 iter=10000, dt=0.001, CFL=0.2, use_local_dt=1, res_NK=1e-4, res_tol=1e-5,
                 reinitialize=0, adj_funcs=[],
                 plot=False, NJ=49, s0=0.5e-2, use_f2py=False,
                 dAl=None, dAu=None):
    '''
    This is a specialized version of the run_case function
    that employs the CST parametrization over a unit chord
    airfoil.
    Mesh warping from the original airfoil is also enabled
    through the dAl and dAu perturbations. This facilitates
    the computation of derivatives.
    NOTE: Unless the warping procedure is failing, it is 
    highly recommended to use run_cst rather than this function.

    Nchord: number of points along the upper/lower skin of the
    airfoil (the total number of points will be 2*Nchord-2)
    '''

    xref = 0.25
    yref = 0.0
    cref = 1.0

    x = (1-np.cos(np.linspace(0, 1, Nchord)*np.pi))/2

    airfoil = am.cstfoil(Au, Al, x)

    xf = airfoil['x_coord']
    yf = airfoil['y_coord']
    maxt = airfoil['max_thickness']
    xmaxt = airfoil['x_max_thickness']
    mint = airfoil['min_thickness']
    maxc = airfoil['max_camber']
    xmaxc = airfoil['x_max_camber']

    # Check if the user requested mesh-warping
    if (dAl is not None) and (dAu is not None):

        airfoil = am.cstfoil(Au+dAu, Al+dAl, x)

        xfw = airfoil['x_coord']
        yfw = airfoil['y_coord']
        maxt = airfoil['max_thickness']
        xmaxt = airfoil['x_max_thickness']
        mint = airfoil['min_thickness']
        maxc = airfoil['max_camber']
        xmaxc = airfoil['x_max_camber']

        dxf = xfw - xf
        dyf = yfw - yf

    else:

        dxf = None
        dyf = None

    results = run_airfoil(xf, yf, alpha, mach, xref, yref, cref,
                          gamma, order, iter, dt, CFL, use_local_dt, res_NK, res_tol,
                          reinitialize, adj_funcs, plot, NJ, s0, use_f2py, dxf, dyf)

    # Propagate derivatives of the mesh warping algorithm
    if (dAl is not None) and (dAu is not None) and (len(adj_funcs) > 0):
        X = results['X']
        Y = results['Y']
        results['grads'] = gridwarp.cst_deriv(Au, Al, dAu, dAl, x, X, Y, adj_funcs, results['grads'])

    results['maxt'] = maxt
    results['xmaxt'] = xmaxt
    results['mint'] = mint
    results['maxc'] = maxc
    results['xmaxc'] = xmaxc

    return results

#=======================================

def write_settings(alpha, mach, xref, yref, cref, gamma, order,
                   iter, dt, CFL, use_local_dt, res_NK, res_tol, bcs, reinitialize=0, adj_funcs=[]):
    '''
    This function writes the settings.txt file, which
    is read by EulerBlock to set up the simulation.
    '''

    # Convert the adjoint function identifier to its corresponding number
    func_dict = {None:0,
                 'cx_ilow':1,
                 'cy_ilow':2,
                 'cm_ilow':3,
                 'cl_ilow':4,
                 'cd_ilow':5,
                 'cx_jlow':6,
                 'cy_jlow':7,
                 'cm_jlow':8,
                 'cl_jlow':9,
                 'cd_jlow':10,
                 'cx_ihigh':11,
                 'cy_ihigh':12,
                 'cm_ihigh':13,
                 'cl_ihigh':14,
                 'cd_ihigh':15,
                 'cx_jhigh':16,
                 'cy_jhigh':17,
                 'cm_jhigh':18,
                 'cl_jhigh':19,
                 'cd_jhigh':20}
    
    #if adj_funcs is None:
    #    adj_func_ids = []
    #else:
    adj_func_ids = [func_dict[func_name] for func_name in adj_funcs]

    lines = []
    lines.append('# alpha (rad)\n')
    lines.append('%g\n'%alpha)
    lines.append('# mach\n')
    lines.append('%g\n'%mach)
    lines.append('# reference x position\n')
    lines.append('%g\n'%xref)
    lines.append('# reference y position\n')
    lines.append('%g\n'%yref)
    lines.append('# reference chord\n')
    lines.append('%g\n'%cref)
    lines.append('# gamma\n')
    lines.append('%g\n'%gamma)
    lines.append('# order (1 for 1st; 2 for 2nd)\n')
    lines.append('%d\n'%order)
    lines.append('# time-stepping strategy (0 for global dt; 1 for local dt based on CFL)\n')
    lines.append('%d\n'%use_local_dt)
    lines.append('# maximum number of time steps\n')
    lines.append('%d\n'%iter)
    lines.append('# dt (only used for global time-steeping)\n')
    lines.append('%g\n'%dt)
    lines.append('# CFL (only used for local time-steeping)\n')
    lines.append('%g\n'%CFL)
    lines.append('# residuals MSE tolerance to switch to NK solver (should be lower than the next residual value)\n')
    lines.append('%g\n'%res_NK)
    lines.append('# residuals MSE tolerance to finish the execution\n')
    lines.append('%g\n'%res_tol)
    lines.append('# BC types for ilow, jlow, ihigh, and jhigh boundaries (1:farfield, 2:wall, 3:periodic)\n')
    lines.append('%d %d %d %d\n'%tuple(bcs))
    lines.append('# Reinitialize solution with solution.vtk? (0-No 1-Yes)\n')
    lines.append('%d\n'%reinitialize)
    lines.append('# First line: Quantity of functions to compute derivatives (use 0 for no derivatives). Next lines: Function flags (1-cx_ilow 2-cy_ilow 3-cm_ilow 4-cl_ilow 05-cd_ilow 06-cx_jlow 07-cy_jlow 08-cm_jlow 09-cl_jlow 10-cd_jlow 11-cx_ihigh 12-cy_ihigh 13-cm_ihigh 14-cl_ihigh 15-cd_ihigh 16-cx_jhigh 17-cy_jhigh 18-cm_jhigh 19-cl_jhigh 20-cd_jhigh)\n')
    lines.append('%d\n'%len(adj_func_ids))
    for adj_func_id in adj_func_ids:
        lines.append('%d\n'%adj_func_id)

    with open('settings.txt','w') as fid:
        fid.writelines(lines)

#=======================================

def read_walls():
    '''
    This script reads the wall.dat file generated by the Fortran code
    with the code results.
    '''

    # Initialize dictionary of results
    wallData = {'ilow':None,
                'jlow':None,
                'ihigh':None,
                'jhigh':None}

    # Read the output file
    with open('wall.dat','r') as fid:
        lines = fid.readlines()

    # Set reading state flags
    store_cp = False
    store_F = False

    # Initialize data fields
    xx = []
    yy = []
    Cp = []
    M = []

    # Loop over each line
    for line in lines:

        # Identity beginning of wall section.
        # If this happens, we get ready to read forces
        # in the next line and reset list of data
        if 'ilow' in line:
            curr_wall = 'ilow'
            store_F = True

        elif 'jlow' in line:
            curr_wall = 'jlow'
            store_F = True

        elif 'ihigh' in line:
            curr_wall = 'ihigh'
            store_F = True

        elif 'jhigh' in line:
            curr_wall = 'jhigh'
            store_F = True
            

        # Skip other types of comments
        elif '#' in line:
            pass

        # Find the integrated forces and prepare to read wall
        # distributions from the next lines
        elif store_F:
            Cx, Cy, Cm, Cl, Cd = map(float,line.split())
            store_F = False
            store_cp = True

        elif store_cp:
            if len(line) > 3:
                xi, yi, Cpi, Mi = map(float,line.split())
                xx.append(xi)
                yy.append(yi)
                Cp.append(Cpi)
                M.append(Mi)

            # The Cp curve ends when we find a blank line.
            # Then we store the data in the appropriate dictionary
            else:
                store_cp = False

                wallData[curr_wall] = {'x':xx,
                                       'y':yy,
                                       'Cp':Cp,
                                       'Mach':M,
                                       'Cx':Cx,
                                       'Cy':Cy,
                                       'Cm':Cm,
                                       'Cl':Cl,
                                       'Cd':Cd,}

                # Reinitialize data fields
                xx = []
                yy = []
                Cp = []
                M = []
    
    return wallData

#=======================================

def read_derivatives(adj_func):
    '''
    This function reads the file generated by the Fortran code
    with derivatives computed by the adjoint method.
    '''

    if adj_func is None:
        return {}

    # Read the output file
    with open('derivatives.dat','r') as fid:
        lines = fid.readlines()

    # Initialize dictionary to hold derivatives
    grad = {}

    found_nodes = False

    # Loop over every line to get results
    for lid,line in enumerate(lines):

        if 'Grid size' in line:
            ni, nj = lines[lid+1].split()
            ni = int(ni)
            nj = int(nj)

        # Check if we found a new function
        elif 'Derivatives with respect to the following function (f)' in line:

            # Store the new function name. The [:-2] is to avoid the \n
            curr_func = lines[lid+1][:-2]

            # Initialize dictionary for this new function
            grad[curr_func] = {'X': np.zeros((ni,nj)),
                               'Y': np.zeros((ni,nj))}

        elif 'df/dMach_inf' in line:
            grad[curr_func]['mach'] = float(lines[lid+1])

        elif 'df/dalpha_inf' in line:
            grad[curr_func]['alpha'] = float(lines[lid+1])

        elif 'i, j, X, Y, df/dX, df/dY' in line:
            found_nodes = True

        elif len(line) < 4:
            found_nodes = False

        elif found_nodes:

            broken_line = line.split(',')

            if len(broken_line) < 6: # We reached the end of the file
                break

            ii, jj, xx, yy, xb, yb = line.split(',')

            ii = int(ii)-1
            jj = int(jj)-1
            xb = float(xb)
            yb = float(yb)

            grad[curr_func]['X'][ii,jj] = xb
            grad[curr_func]['Y'][ii,jj] = yb

    return grad

#=======================================

def postprocess(x, y, alpha, mach, wallData=None, plot=False):
    '''
    x, y are airfoil coordinates (same inputs as run_case)
    '''

    # Read the output file
    if wallData is None:
        wallData = read_walls()
    
    xx = wallData['jlow']['x']
    Cp = wallData['jlow']['Cp']
    Mach = wallData['jlow']['Mach']
    Cx = wallData['jlow']['Cx']
    Cy = wallData['jlow']['Cy']
    CM = wallData['jlow']['Cm']
    CL = wallData['jlow']['Cl']
    CD = wallData['jlow']['Cd']

    if plot:
        fig = plt.figure()
        plt.subplot(311)
        plt.plot(xx,Cp)
        plt.ylabel('Cp')
        plt.title(r'$\alpha=%.2f^o \quad M=%g \quad C_L=%.4f \quad C_D=%.5f \quad C_M=%.4f$'%(alpha*180.0/np.pi, mach, CL, CD, CM))
        plt.gca().invert_yaxis()
        plt.subplot(312)
        plt.plot(xx,Mach)
        plt.ylabel('Mach')
        plt.subplot(313)
        plt.plot(x,y)
        plt.ylabel('y')
        plt.xlabel('x')
        plt.axis('equal')
        plt.show()

    return CL, CD, CM, xx, Cp, Mach

#=======================================