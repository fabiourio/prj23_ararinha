'''
This module contains functions to generate airfoil geometries

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2022
'''

import numpy as np
import matplotlib.pyplot as plt

#=======================================

def nacafoil(naca_string, x, sharp_te=False, plot=False):
    '''
    This function generates coordinates for a NACA 4- or 5-digit airfoil.
    We used the modified version where the last coefficient of the
    thickness equation is adjusted to give a sharp trailing edge.

    INPUTS
    naca_string: string -> string with 4 naca digits (e.g. '2412')
    x: real(n) -> array of positions going from 0 to cref
    '''

    # Get reference chord
    cref = x[-1]

    # Normalize x
    xn = x/cref

    # Thickness distribution coefficients
    tt0 =  0.2969
    tt1 = -0.1260
    tt2 = -0.3516
    tt3 =  0.2843
    if sharp_te:
        tt4 = -0.1036
    else:
        tt4 = -0.1015

    if len(naca_string) == 4:

        # Get airfoil information from digits
        mm = int(naca_string[0])/100.0
        pp = int(naca_string[1])/10.0
        tc = int(naca_string[2:])/100.0

        # Thickness distribution
        yt = 5*tc*(tt0*np.sqrt(xn) + tt1*xn + tt2*xn**2 + tt3*xn**3 + tt4*xn**4)

        # Camber distribution
        if pp > 0:

            yc1 = mm/pp**2*(2*pp*xn - xn**2)
            yc2 = mm/(1-pp)**2*(1 - 2*pp + 2*pp*xn - xn**2)
            yc = (xn <= pp)*yc1 + (xn > pp)*yc2
    
            dyc1dx = mm/pp**2*(2*pp - 2*xn)
            dyc2dx = mm/(1-pp)**2*(2*pp - 2*xn)
            dycdx = (xn <= pp)*dyc1dx + (xn > pp)*dyc2dx

        else:

            yc = np.zeros_like(xn)
            dycdx = np.zeros_like(xn)

    elif len(naca_string) == 5:

        # Get airfoil information from digits
        LL = int(naca_string[0])
        pp = int(naca_string[1])/20.0
        ss = int(naca_string[2])
        tc = int(naca_string[3:])/100.0

        # Thickness distribution
        yt = 5*tc*(tt0*np.sqrt(xn) + tt1*xn + tt2*xn**2 + tt3*xn**3 + tt4*xn**4)

        # Camber distribution
        if pp > 0:

            if ss == 0:

                rr = np.interp(pp, [0.05, 0.10, 0.15, 0.20, 0.25],
                                   [0.058, 0.126, 0.2025, 0.2900, 0.391])

                k1 = np.interp(pp, [0.05, 0.10, 0.15, 0.20, 0.25],
                                   [361.40, 51.640, 15.957, 6.643, 3.230])

                yc1 = k1/6*(xn**3 - 3*rr*xn**2 + rr**2*(3-rr)*xn)
                yc2 = k1/6*rr**3*(1-xn)
    
                dyc1dx = k1/6*(3*xn**2 - 6*rr*xn + rr**2*(3-rr))
                dyc2dx = -k1/6*rr**3
                

            elif ss == 1:

                rr = np.interp(pp, [0.10, 0.15, 0.20, 0.25],
                                   [0.130, 0.217, 0.318, 0.441])

                k1 = np.interp(pp, [0.10, 0.15, 0.20, 0.25],
                                   [51.990, 15.793, 6.520, 3.191])

                kr = np.interp(pp, [0.10, 0.15, 0.20, 0.25],
                                   [0.000764, 0.00677, 0.0303, 0.1355])

                yc1 = k1/6*((xn-rr)**3 - kr*(1-rr)**3*xn - rr**3*xn + rr**3)
                yc2 = k1/6*(kr*(xn-rr)**3 - kr*(1-rr)**3*xn - rr**3*xn + rr**3)
    
                dyc1dx = k1/6*(3*(xn-rr)**2 - kr*(1-rr)**3 - rr**3)
                dyc2dx = k1/6*(3*kr*(xn-rr)**2 - kr*(1-rr)**3 - rr**3)

            yc = (xn <= rr)*yc1 + (xn > rr)*yc2
            dycdx = (xn <= rr)*dyc1dx + (xn > rr)*dyc2dx

            # Apply scaling factor according to the design lift coefficient (first digit).
            # The table is for CL_des = 0.3, what corresponds to LL=2.
            # If LL=0, then CL_des = 0.05. We linearly scale within this range.
            # Airfoiltools limit values to CL=0.05, what corresponds to yc_factor=1.0/6.0.
            yc_factor = max(1.0/6.0,LL/2.0)   #(1.0 + 2.5*LL)/6.0
            yc = yc_factor*yc
            dycdx = yc_factor*dycdx

        else:

            yc = np.zeros_like(xn)
            dycdx = np.zeros_like(xn)

    # Get skin coordinates
    xl,yl,xu,yu = tt_cc_to_skins(xn,yt,yc,dycdx,t_method='perpendicular')
    xf, yf = join_skins(xl,yl,xu,yu)
    xf = xf*cref
    yf = yf*cref

    # Get other properties
    maxt = 2*np.max(yt)
    xmaxt = x[np.argmax(yt)]/cref
    mint = 2*np.min(yt)
    maxc = np.max(yc)
    xmaxc = x[np.argmax(yc)]/cref

    # Build output
    airfoil = build_airfoil_dict(xf, yf, maxt, xmaxt, mint, maxc, xmaxc, xl, yl, xu, yu)

    # Plot airfoil if requested by the user
    if plot:
        plot_airfoil(airfoil)

    return airfoil

#=======================================

def cstfoil(Au, Al, x, N1=0.5, N2=1.0, tte=0.0, plot=False):
    '''
    This function generates the airfoil ordinates at the locations
    described by x.
    Based on the CST description given by:
    Inverse Airfoil Design Utilizing CST Parameterization.
    Kevin A. Lane and David D. Marshall.
    48th AIAA Aerospace Sciences Meeting Including the New Horizons Forum and Aerospace Exposition.
    4 - 7 January 2010, Orlando, Florida.
    AIAA 2010-1228

    INPUTS
    Au: real(nu) -> weights of the upper surface.
    Al: real(nl) -> weights of the lower surface (can be a different quantity of the upper surface).
    x: real(n) -> array of positions going from 0 to cref
    tte: real -> trailing edge thickness

    OUTPUTS
    maxt: maximum thickness
    xmaxt: chordwise position maximum thickness
    mint: minimum thickness (should be >= zero for a feasible airfoil)
    maxc: maximum camber
    xmaxc: chordwise position maximum camber
    '''

    from scipy.special import comb

    # Get the length of the profile
    cref = max(x)

    # Make sure x-coordinates are normalized
    psi = x/cref

    # Compute class function
    C = psi**N1*(1.0-psi)**N2

    # Get order of Bernstein polynomials, which is
    # one less than the number of coefficients
    Nu = len(Au)-1
    Nl = len(Al)-1

    # Compute shape functions
    Su = 0.0
    Sl = 0.0
    for ii in range(Nu+1):
        K = comb(Nu,ii)
        Su = Su + Au[ii]*K*psi**ii*(1.0-psi)**(Nu-ii)
    for ii in range(Nl+1):
        K = comb(Nl,ii)
        Sl = Sl + Al[ii]*K*psi**ii*(1.0-psi)**(Nl-ii)

    # Determine coordinates
    csiu = C*Su + psi*0.5*tte/cref
    csil = C*Sl - psi*0.5*tte/cref

    # Compute airfoil properties
    tt = csiu-csil # Thickness distribution
    cond1 = x > 0.05
    cond2 = x < 0.97
    chosen_id = np.where(cond1 & cond2)
    tt = tt[chosen_id]
    imaxt = np.argmax(np.real(tt))
    maxt = tt[imaxt]
    xmaxt = psi[imaxt]
    mint = min(tt)

    # Use KS function for maximum thickness.
    # We use this instead of the max function so that it has
    # a continuous variation for gradient-based optimization.
    rhoKS = 500
    maxt = maxt + 1/rhoKS*np.log(np.sum(np.exp(rhoKS*(tt-maxt))))
    mint = mint - 1/rhoKS*np.log(np.sum(np.exp(rhoKS*(mint-tt))))

    cc = 0.5*(csiu+csil) # Camber distribution
    imaxc = np.argmax(np.abs(np.real(cc)))
    maxc = cc[imaxc]
    xmaxc = psi[imaxc]

    # Denormalize results
    xu = psi*cref
    xl = psi*cref
    yu = csiu*cref
    yl = csil*cref

    xf, yf = join_skins(xl,yl,xu,yu)
        
    # Build output
    airfoil = build_airfoil_dict(xf, yf, maxt, xmaxt, mint, maxc, xmaxc, xl, yl, xu, yu)

    # Plot airfoil if requested by the user
    if plot:
        plot_airfoil(airfoil)

    return airfoil

#=======================================

def fit_cst(xu,yu,xl,yl,N,options=None,
            force_sharp_te=False, force_symmetric=False):
    '''
    This function fits CST parameters for given airfoil coordinates.

    INPUTS
    xu: real(nu) -> array of chordwise coordinates for upper skin samples
    yu: real(nu) -> upper skin coordinates
    xu: real(nl) -> array of chordwise coordinates for lower skin samples
    yl: real(nl) -> lower skin coordinates
    N: integer -> number of desired CST parameters (N for upper and N for lower skin)

    OUTPUTS
    Au: real(N) -> upper skin CST parameters
    Al: real(N) -> lower skin CST parameters
    '''

    from scipy.optimize import minimize

    print()
    print('Fitting CST coefficients')

    # Define objetive function that computes MSE error between
    # computed and provided coordinates.
    def mse_error(aa):

        # The design variables for this optimization are:
        # Al - CST coefficients of the lower skin
        # dAu - (Au-Al), Difference between CST coefficients of upper and lower skin
        # tte - trailing edge thickness
        #
        # We use dAu instead of Au so that we constraint all Au values
        # to be greater that their corresponding Al values by just
        # setting dAu >= 0 in the optimization bounds.

        Al = aa[:N]
        dAu = aa[N:-1]
        tte = aa[-1]

        if force_symmetric:
            Au = -Al
        else:
            Au = Al + dAu

        # Generate samples at the original lower skin location
        airfoil = cstfoil(Au, Al, xl, N1=0.5, N2=1.0, tte=tte, plot=False)
        yl_trial = airfoil['y_lower']

        # Generate samples at the original upper skin location
        airfoil = cstfoil(Au, Al, xu, N1=0.5, N2=1.0, tte=tte, plot=False)
        yu_trial = airfoil['y_upper']
        
        mse = np.sum((yl-yl_trial)**2) + np.sum((yu-yu_trial)**2)

        return mse

    # Define gradient of MSE for optimization
    mse_error_grad = complex_step_gradient(mse_error)

    # Initial guess for optimization
    Al0 = -np.ones(N)
    dAu0 = 2*np.ones(N)
    tte0 = yu[-1]-yl[-1]
    aa0 = np.hstack([Al0, dAu0, tte0])

    # Design variable bounds
    min_val = 1e-10
    if force_symmetric:
        dAu_max = 1e-10
    else:
        dAu_max = None
    
    if force_sharp_te:
        tte_max = min_val
    else:
        tte_max = None

    Al_bounds = [(None, None)]*N
    dAu_bounds = [(0.0, dAu_max)]*N
    tte_bounds = [(0.0, tte_max)]
    bounds = Al_bounds + dAu_bounds + tte_bounds

    # Optimization options
    if options is None:
        options = {'ftol':1e-12}

    results = minimize(mse_error, aa0, method='SLSQP', jac=mse_error_grad,
                       bounds=bounds, options=options)

    print('Optimization log:')
    print(results)

    Al = results.x[:N]
    dAu = results.x[N:-1]
    tte = results.x[-1]

    if force_symmetric:
        Au = -Al
    else:
        Au = Al + dAu

    print()
    print('Final coefficients')
    print('Au:',Au)
    print('Al:',Al)
    print('tte:',tte)
    print('Coordinates MSE:',results.fun)
    print()

    return Au, Al, tte

#=======================================

def strefoil(xx, rle, xtmax, tmax, phi,
             theta, epsilon, xcmax, cmax, ctmax,
             cte=0, tte=0, plot=False, t_method='vertical'):
    '''
    This function generates airfoils using the Streshinsky's parametrization.
    Please see:
    STRESHINSKY, J. R., OVCHARENKO, V. V.,
    Aerodynamic Design Transonic Wing Using CFD and Optimization Methods
    ICAS-94.2.1.4
    https://www.icas.org/icas_archive/ICAS1994/ICAS-94-2.1.4.pdf

    The linear systems are further detailed in:
    SECCO, N. R.,
    Training Artificial Neural Networks to Predict Aerodynamic
    Coefficients of Airliner Wing-Fuselage Configurations. Master Thesis.
    Instituto Tecnologico de Aeronautica, 2014.
    http://www.bdita.bibl.ita.br/tesesdigitais/lista_resumo.php?num_tese=66533

    Note that there is a typo in the phi definition. We flipped its sign here.

    xx: array of floats -> monotonic sequence ranging from 0 to 1 where coordinates are sampled
    rle: leading edge radius
    xtmax: maximum thickness position
    tmax: maximum thickness
    phi: trailing edge angle measured between upper and lower skins
    theta: leading edge camber angle
    epsilon: trailing edge camber angle
    xcmax: maximum camber position
    cmax: maximum camber
    ctmax: camber at maximum thickness position
    cte: trailing edge camber (distance from chordline)
    tte: trailing edge thickness
    '''

    # Assemble the thickness linear system
    # Note that there is a typo in the phi definition.
    # We flipped its sign in the nt definition.
    Mt = np.array([[1, 0, 0, 0, 0],
                   [2, 2, 2, 2, 2],
                   [0.5/np.sqrt(xtmax), 1, 2*xtmax, 3*xtmax**2, 4*xtmax**3],
                   [2*np.sqrt(xtmax), 2*xtmax, 2*xtmax**2, 2*xtmax**3, 2*xtmax**4],
                   [1, 2, 4, 6, 8]])
    nt = np.array([np.sqrt(2*rle), tte, 0, tmax, -phi])

    # Assemble the camber linear system
    Mc = np.array([[1, 0, 0, 0, 0, 0],
                   [1, 2, 3, 4, 5, 6],
                   [1, 2*xcmax, 3*xcmax**2, 4*xcmax**3, 5*xcmax**4, 6*xcmax**5],
                   [xcmax, xcmax**2, xcmax**3, xcmax**4, xcmax**5, xcmax**6],
                   [xtmax, xtmax**2, xtmax**3, xtmax**4, xtmax**5, xtmax**6],
                   [1, 1, 1, 1, 1, 1]])
    nc = np.array([theta, epsilon, 0, cmax, ctmax, cte])

    # Solve linear systems to obtain polynomials
    aa = np.linalg.solve(Mt, nt)
    bb = np.linalg.solve(Mc, nc)

    # Thickness and camber distributions
    tt = aa[0]*np.sqrt(xx) + aa[1]*xx + aa[2]*xx**2 + aa[3]*xx**3 + aa[4]*xx**4
    cc = bb[0]*xx + bb[1]*xx**2 + bb[2]*xx**3 + bb[3]*xx**4 + bb[4]*xx**5 + bb[5]*xx**6
    dcdx = bb[0] + 2*bb[1]*xx + 3*bb[2]*xx**2 + 4*bb[3]*xx**3 + 5*bb[4]*xx**4 + 6*bb[5]*xx**5

    # Get skin coordinates
    xl,yl,xu,yu = tt_cc_to_skins(xx,tt,cc,dcdx,t_method)
    xf, yf = join_skins(xl,yl,xu,yu)

    # Get other properties
    imaxt = np.argmax(tt)
    maxt = 2.0*tt[imaxt]
    xmaxt = xx[imaxt]
    mint = np.min(tt)
    imaxc = np.argmax(np.abs(cc))
    maxc = 2.0*cc[imaxc]
    xmaxc = xx[imaxc]

    # Build output
    airfoil = build_airfoil_dict(xf, yf, maxt, xmaxt, mint, maxc, xmaxc, xl, yl, xu, yu)

    # Plot airfoil if requested by the user
    if plot:
        plot_airfoil(airfoil)

    return airfoil

#=======================================

def fit_stre(xu,yu,xl,yl,opt_options=None,
             force_sharp_te=False, force_symmetric=False):
    '''
    This function fits Strechinsky's parameters for given airfoil coordinates.

    INPUTS
    xu: real(nu) -> array of chordwise coordinates for upper skin samples
    yu: real(nu) -> upper skin coordinates
    xu: real(nl) -> array of chordwise coordinates for lower skin samples
    yl: real(nl) -> lower skin coordinates

    OUTPUTS
    Au: real(N) -> upper skin CST parameters
    Al: real(N) -> lower skin CST parameters
    '''

    from scipy.optimize import minimize

    print()
    print("Fitting Strechinsky's coefficients")

    # Define objetive function that computes MSE error between
    # computed and provided coordinates.
    def mse_error(params):
        
        # Split design variables
        [rle,
         xtmax,
         tmax,
         phi,
         theta,
         epsilon,
         xcmax,
         cmax,
         ctmax,
         cte,
         tte] = params

        # Generate samples at the original lower skin location
        # We use the 'vertical' t_method to ensure that coordinated
        # are sampled at the right location
        airfoil = strefoil(xl, rle, xtmax, tmax, phi,
                           theta, epsilon, xcmax, cmax, ctmax,
                           cte, tte, plot=False, t_method='vertical')
        yl_trial = airfoil['y_lower']

        # Generate samples at the original upper skin location
        airfoil = strefoil(xu, rle, xtmax, tmax, phi,
                           theta, epsilon, xcmax, cmax, ctmax,
                           cte, tte, plot=False, t_method='vertical')
        yu_trial = airfoil['y_upper']
        
        mse = np.sum((yl-yl_trial)**2) + np.sum((yu-yu_trial)**2)

        return mse
    
    # Define gradient of MSE for optimization
    mse_error_grad = complex_step_gradient(mse_error)

    # Initial guess for optimization
    params0 = [0.005, # rle
               0.4, # xtmax
               0.12, # tmax
               0.10, # phi
               0.0, # theta
               0.0, # epsilon
               0.5, # xcmax
               0.0, # cmax
               0.0, # ctmax
               0.0, # cte
               0.0] # tte

    # Design variable bounds
    min_val = 1e-10
    if force_sharp_te:
        min_tte = 0.0
        max_tte =  min_val
    else:
        min_tte = 0.0
        max_tte = 0.1

    if force_symmetric:
        min_theta = -min_val
        max_theta =  min_val
        min_epsilon = -min_val
        max_epsilon =  min_val
        min_xcmax = 0.5-min_val
        max_xcmax = 0.5+min_val
        min_cmax = -min_val
        max_cmax =  min_val
        min_ctmax = -min_val
        max_ctmax =  min_val
        min_cte = -min_val
        max_cte =  min_val
    else:
        min_theta = -np.pi/4
        max_theta =  np.pi/4
        min_epsilon = -np.pi/4
        max_epsilon =  np.pi/4
        min_xcmax = 0.1
        max_xcmax = 0.9
        min_cmax = -0.1
        max_cmax =  0.1
        min_ctmax = -0.1
        max_ctmax =  0.1
        min_cte = -0.1
        max_cte =  0.1

    bounds =[(0.0, 0.1),   # rle
             (0.1, 0.9),    # xtmax
             (0.0, 0.4),   # tmax
             (0.0, np.pi/4),   # phi
             (min_theta, max_theta),  # theta
             (min_epsilon, max_epsilon),  # epsilon
             (min_xcmax, max_xcmax),    # xcmax
             (min_cmax, max_cmax),   # cmax
             (min_ctmax, max_ctmax),   # ctmax
             (min_cte, max_cte),  # cte
             (min_tte, max_tte)]   # tte

    # Optimization options
    if opt_options is None:
        opt_options = {'ftol':1e-12}

    results = minimize(mse_error, params0, method='SLSQP', jac=mse_error_grad,
                       bounds=bounds, options=opt_options, tol=1e-12)

    # I also tested the differential evolution, but I got the same answer at additional cost.
    #from scipy.optimize import differential_evolution
    #results = differential_evolution(mse_error, bounds, maxiter=1000, popsize=15, tol=1e-6, mutation=(0.5, 1), recombination=0.7, disp=False, polish=True, atol=0)

    print('Optimization log:')
    print(results)

    [rle,
     xtmax,
     tmax,
     phi,
     theta,
     epsilon,
     xcmax,
     cmax,
     ctmax,
     cte,
     tte] = results.x

    print()
    print('Final coefficients')
    print('rle:',rle)
    print('xtmax:',xtmax)
    print('tmax:',tmax)
    print('phi:',phi)
    print('theta:',theta)
    print('epsilon:',epsilon)
    print('xcmax:',xcmax)
    print('cmax:',cmax)
    print('ctmax:',ctmax)
    print('cte:',cte)
    print('tte:',tte)
    print('Coordinates MSE:',results.fun)
    print()

    return rle, xtmax, tmax, phi, theta, epsilon, xcmax, cmax, ctmax, cte, tte

#=======================================

def export_airfoil(airfoil, title='AIRFOIL', filename='airfoil.dat', reverse=False, close_te=True):
    '''
    This function generates a text file with airfoil coordinates
    according to the Xfoil input format (starts at trailing edge and loops around
    the airfoil from the lower skin first).
    The reverse keyword changes the order of points to loop over the upper skin first.
    If title is None, the file will have only coordinates.
    '''

    # Gather complete airfoil coordinates
    xf = airfoil['x_coord']
    yf = airfoil['y_coord']

    # Initialize list containing lines to be written in the file
    lines = []

    # Add airfoil name
    if title is not None:
        lines.append(title+'\n')

    # Add each x,y coordinate
    for ii in range(len(xf)):
        lines.append('%.13f %.13f\n'%(xf[ii], yf[ii]))

    # Repeat first point to close the airfoil
    if close_te:
        lines.append('%.13f %.13f\n'%(xf[0], yf[0]))

    # Write file
    with open(filename, 'w') as fid:
        if not reverse:
            fid.writelines(lines)
        else:
            fid.writelines(lines[0:1] + lines[-1:0:-1])

    print()
    print('Exported airfoil coordinates to:',filename)

#=======================================

def cos_bunching(xmin, xmax, NN):
    '''
    Simple cosenoidal bunching (more refined at the beginning and the end)
    between xmin and xmax.
    '''

    tt = np.linspace(0, np.pi, NN)
    xx = 0.5*(1-np.cos(tt))

    xx = xmin + xx*(xmax - xmin)

    return xx

#=======================================

def get_bunching(NN, le_sep=None):
    '''
    This functions generates samples of x coordinates between 0 and 1
    to use as input for the airfoil parametrization functions.
    The usar can force a point in le_sep to facilitate splits for mesh generation.
    '''

    if le_sep is None:
        xx = cos_bunching(0,1,NN)
    else:
        xx = np.hstack([cos_bunching(0, le_sep, int(20*NN*le_sep)), cos_bunching(le_sep, 1.0, int(NN*(1-5*le_sep)))[1:]])

    return xx

#=======================================

def plot_airfoil(airfoil):

    # Gather complete airfoil coordinates
    xf = airfoil['x_coord']
    yf = airfoil['y_coord']

    fig = plt.figure()
    plt.plot(xf,yf,'-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()

#=======================================

def read_coordinates(file_path, scale=1.0):
    '''
    This function reads airfoil coordinates from a text file.
    We assume that this text file has two columns (x and y coordinates).
    '''

    xx = []
    yy = []

    print()
    print('Reading airfoil coordinates from file:', file_path)

    with open(file_path, 'r') as file:

        for line in file:

            try:
                x_curr, y_curr = map(float, line.split())
            except:
                print('Failed to read line:')
                print(line[:-1])
                print('Check if coordinates are correctly loaded.')
                continue

            # Apply scalings if necessary
            x_curr = x_curr*scale
            y_curr = y_curr*scale

            xx.append(x_curr)
            yy.append(y_curr)
    
    print('Number of coordinates:',len(xx))
    print()

    xx = np.array(xx)
    yy = np.array(yy)

    return xx, yy

#=======================================

def build_airfoil_dict(xf, yf, maxt, xmaxt, mint, maxc, xmaxc, xl, yl, xu, yu):
    '''
    This function builds the airfoil dictionary which is the
    main type of output from airfoil_mod.
    '''

    airfoil = {'x_coord': xf,
               'y_coord': yf,
               'max_thickness': maxt,
               'x_max_thickness': xmaxt,
               'min_thickness': mint,
               'max_camber': maxc,
               'x_max_camber': xmaxc,
               'x_lower': xl,
               'y_lower': yl,
               'x_upper': xu,
               'y_upper': yu}
    
    return airfoil

#=======================================

def split_skins(xx, yy):
    '''
    This functions splits the upper and lower skin coordinates from
    the coordinate vector given in Selig format.
    '''

    # We will assign the first coordinates to the lower skin
    # until the changes in x from consecutive points change directions,
    # what we expect to be the leading edge.
    # We will check if this first set of points actually belongs to the lower
    # skin later on and make the necessary adjustments.

    first_delta_x = xx[1]-xx[0]

    if first_delta_x < 0: # Coordinates starting at trailing edge
        edge_index = np.where(xx[1:]-xx[:-1] > 0.0)[0][0]
    else: # Coordinates starting at leading edge
        edge_index = np.where(xx[1:]-xx[:-1] < 0.0)[0][0]

    # Initialize vectors for both skins
    xl = xx[:edge_index+1]
    yl = yy[:edge_index+1]
    xu = xx[edge_index:]
    yu = yy[edge_index:]

    # Revert coordinates so that they always go from leading to trailing edge
    if first_delta_x < 0:
        xl = xl[::-1]
        yl = yl[::-1]
    else:
        xu = xu[::-1]
        yu = yu[::-1]

    # Check average values of each skin to identify the upper and lower ones
    yl_avg = np.average(yl)
    yu_avg = np.average(yu)

    if yl_avg > yu_avg:
        xu,xl = xl,xu
        yu,yl = yl,yu

    return xl, yl, xu, yu

#=======================================

def join_skins(xl,yl,xu,yu):
    '''
    This function merges coordinates in a single array
    that begins at the trailing edge and loops around
    the leading edge, starting with the lower skin.
    '''

    xf = np.hstack([xl[::-1], xu[1:]])
    yf = np.hstack([yl[::-1], yu[1:]])

    return xf, yf

#=======================================

def tt_cc_to_skins(xx,tt,cc,dcdx,t_method='perpendicular'):
    '''
    This function computes the upper and lower skin coordinates
    based on the thickness and camber distribution along the airfoil.
    Two thickness methods are available:
    'perpendicular': thickness is applied perpendicular to the camber line
    'vertical': thickness is applied perpendicular to the chord line (x-direction)
    '''

    if t_method == 'perpendicular':

        theta = np.arctan(dcdx)
        xu = xx - tt*np.sin(theta)
        xl = xx + tt*np.sin(theta)
        yu = cc + tt*np.cos(theta)
        yl = cc - tt*np.cos(theta)

    elif t_method == 'vertical':
        
        xu = xx.copy()
        xl = xx.copy()
        yu = cc + tt
        yl = cc - tt

    return xl,yl,xu,yu

#=======================================

def complex_step_gradient(func,hh=1e-20):
    '''
    This function returns another function
    that gives the gradients of func using complex-step.
    '''

    def grad_func(xx):

        grad = np.zeros_like(xx)

        for ii in range(len(xx)):

            xx_hh = np.array(xx, dtype=complex)
            xx_hh[ii] = xx[ii] + 1j*hh

            f_hh = func(xx_hh)

            grad[ii] = np.imag(f_hh)/hh

        return grad
    
    return grad_func