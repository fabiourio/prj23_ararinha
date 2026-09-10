'''
2D Grid generator
Written by:
Ney Rafael Secco
ney@ita.br
Jan 2016
Rev: Aug 2022

If you use this code, please cite:
SÊCCO, N. R.; KENWAY, G. K. W. ; HE, P. ; MADER, C. ; MARTINS, J. R. R. A.
Efficient Mesh Generation and Deformation for Aerodynamic Shape Optimization.
AIAA Journal, v. 59, p. 1151-1168, 2021. DOI: http://dx.doi.org/10.2514/1.J059491. 

TO DO

-set up scaling factor based on distance
-blend the angle-based dissipation coefficient
'''

def march_layer(r0, bc1, bc2, s0, NJ, extension, \
                epsE0 = 1.0, theta = 0.0, \
                alphaP0 = 0.25, num_smoothing_passes = 0, \
                nuArea = 0.16, num_area_passes = 0, \
                sigmaSplay = 0.3, \
                metricCorrection = False, \
                ratioGuess = 2):

    # This function will march N steps

    # IMPORTS
    from numpy import zeros, copy

    print('2D Grid generator')
    print('Written by:')
    print('Ney Rafael Secco')
    print('ney@ita.br')
    print('Jan 2016')
    print('Rev: Aug 2022')
    print('')

    # Initialize 2D array that contains all the grid points
    R = zeros((NJ,len(r0)))
    R[0,:] = r0
    r_prevprev = r0

    # Initialize step size and total marched distance
    s = s0
    s_tot = 0

    # Find the characteristic radius of the mesh
    radius = findRadius(r0)

    # Find the desired marching distance
    s_max = radius*(extension-1.0)

    # Compute the growth ratio necessary to match this distance
    growth = findRatio(s_max, s0, NJ, ratioGuess)

    # Print growth ratio
    print('Growth ratio: ',growth)
    
    '''
    ==========================================
    Now let's define the functions that do the heavy work.
    These functions are the actual hyperbolic generation
    calculations.
    I define them here so that the several optional parameters
    are readily available to them.
    ==========================================
    '''

    def sub_iteration(r_prev, r_prevprev, s, s_tot, layer_index):

        # IMPORTS
        from numpy.linalg import solve

        # Smooth coordinates
        r_prev = smoothing(r_prev,layer_index+2)

        # Compute area factor
        A = area_factor(r_prev,s)

        # Generate matrices of the linear system
        K,f = compute_matrices(r_prev,r_prevprev,A,layer_index,s_tot)

        # Solve the linear system
        dr = solve(K,f)

        # Update r
        r_next = r_prev + dr

        # RETURNS
        return r_next

    def area_factor(r0,s):

        # This function computes the area factor distribution

        # IMPORT
        from numpy import hstack, reshape
        from numpy.linalg import norm

        # Extrapolate the end points
        r0minus1 = 2*r0[0:2] - r0[2:4]
        r0plus1 = 2*r0[-2:] - r0[-4:-2]
        r0_extrap = hstack([r0minus1, r0, r0plus1])

        # Reshape so we have 2D array
        R0_extrap = reshape(r0_extrap,(-1,2)).T

        # Compute the distance of each node to its neighbors
        d = 0.5*(norm(R0_extrap[:,1:-1] - R0_extrap[:,:-2],axis=0) + norm(R0_extrap[:,2:] - R0_extrap[:,1:-1],axis=0))

        # Adjust the splay
        d[0] = d[0]*(1-sigmaSplay)
        d[-1] = d[-1]*(1-sigmaSplay)

        # Multiply distances by the step size to get the areas
        A = d*s

        # Do the requested number of averagings
        for index in range(num_area_passes):
            # Store previous values
            Aplus = A[1]
            Aminus = A[-2]
            # Do the averaging
            A[1:-1] = (1-nuArea)*A[1:-1] + nuArea/2*(A[:-2] + A[2:])
            A[0] = (1-nuArea)*A[0] + nuArea*Aplus
            A[-1] = (1-nuArea)*A[-1] + nuArea*Aminus

        # RETURNS
        return A

    def metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index):

        # This function applies the metric correction algorithm to update x0_eta and y0_eta

        # Compute distance vectors (Eq. 6.9)
        r_plus = r0_next - r0_curr
        r_minus = r0_prev - r0_curr

        # Compute corrected zeta derivatives (Eq. 7.2)
        r_zeta = 0.25*(abs(r_plus) + abs(r_minus))*(r_plus/abs(r_plus) - r_minus/abs(r_minus))
        xp_zeta = r_zeta[0]
        yp_zeta = r_zeta[1]

        # Find the adjusted eta derivatives (Eq. 7.1)
        gamma_p = xp_zeta*xp_zeta + yp_zeta*yp_zeta
        xp_eta = -dA_curr/gamma_p*yp_zeta
        yp_eta = dA_curr/gamma_p*xp_zeta

        # Compute metric correction factor
        nu_mc = 2**(-layer_index)

        # Apply metric correction (Eq. 7.3)
        x0_eta = (1-nu_mc)*x0_eta + nu_mc*xp_eta
        y0_eta = (1-nu_mc)*y0_eta + nu_mc*yp_eta

        # RETURNS
        return x0_eta, y0_eta

    def compute_matrices(r0,rm1,dA,layer_index,s_tot):

        # This function computes the derivatives r_zeta and r_eta for a group
        # of coordinates given in r.
        # It shoud be noted that r is a 1d array that contains x and y
        # coordinates of all points in a layer:
        # r0 = [x1, y1, x2, y2, x3, y3, ... ]

        # IMPORTS
        from numpy import zeros, sqrt, eye, array, arctan2, pi
        from numpy.linalg import norm

        # Find number of nodes
        NI = int(len(r0)/2)

        # Initialize arrays
        K = zeros((2*NI,2*NI))
        f = zeros(2*NI)

        # Define assembly routines
        def central_matrices(prev_index,curr_index,next_index):

            # Using central differencing for zeta = 2:NI-1
            x0_zeta = 0.5*(r0[2*(next_index)] - r0[2*(prev_index)])
            y0_zeta = 0.5*(r0[2*(next_index)+1] - r0[2*(prev_index)+1])

            # Computing gamma
            gamma = x0_zeta*x0_zeta + y0_zeta*y0_zeta

            # Use forward differencing for eta
            x0_eta = -dA[curr_index]/gamma*y0_zeta
            y0_eta = dA[curr_index]/gamma*x0_zeta

            if metricCorrection:
                r0_prev = r0[2*(prev_index):2*(prev_index)+2]
                r0_curr = r0[2*(curr_index):2*(curr_index)+2]
                r0_next = r0[2*(next_index):2*(next_index)+2]
                dA_curr = dA[curr_index]
                x0_eta, y0_eta = metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index)

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[2*(next_index):2*(next_index)+2]-rm1[2*(index):2*(index)+2]) + norm(rm1[2*(prev_index):2*(prev_index)+2]-rm1[2*(index):2*(index)+2])
            dden = norm(r0[2*(next_index):2*(next_index)+2]-r0[2*(index):2*(index)+2]) + norm(r0[2*(prev_index):2*(prev_index)+2]-r0[2*(index):2*(index)+2])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = give_angle(r0[2*(prev_index):2*(prev_index)+2], \
                               r0[2*(curr_index):2*(curr_index)+2], \
                               r0[2*(next_index):2*(next_index)+2])

            # Sharp convex corner detection
            if angle > 240*pi/180: # Corner detected

                # Populate matrix with Eq 8.3
                K[2*curr_index:2*curr_index+2,2*next_index:2*next_index+2] = -eye(2)
                K[2*curr_index:2*curr_index+2,2*curr_index:2*curr_index+2] = 2*eye(2)
                K[2*curr_index:2*curr_index+2,2*prev_index:2*prev_index+2] = -eye(2)
                f[2*curr_index:2*curr_index+2] = array([0,0])

            else:

                # Compute other auxiliary variables
                alpha = x0_zeta*x0_eta - y0_zeta*y0_eta
                beta = x0_zeta*y0_eta + x0_eta*y0_zeta

                # Compute C = Binv*A
                C = array([[alpha, beta],[beta, -alpha]])/gamma

                # Compute smoothing coefficients
                epsE, epsI = dissipation_coefficients(layer_index, gamma, x0_eta, y0_eta, d_sensor, angle)

                # Compute RHS components
                Binvg = array([x0_eta, y0_eta])
                De = epsE*(r0[2*(prev_index):2*(prev_index)+2] - 2*r0[2*(index):2*(index)+2] + r0[2*(next_index):2*(next_index)+2])

                # Compute block matrices
                L_block = 0.5*(1+theta)*C - epsI*eye(2)
                M_block = (1 + 2*epsI)*eye(2)
                N_block = -0.5*(1+theta)*C - epsI*eye(2)
                f_block = Binvg + De

                # Populate matrix
                K[2*(curr_index):2*(curr_index)+2,2*(next_index):2*(next_index)+2] = L_block
                K[2*(curr_index):2*(curr_index)+2,2*(curr_index):2*(curr_index)+2] = M_block
                K[2*(curr_index):2*(curr_index)+2,2*(prev_index):2*(prev_index)+2] = N_block
                f[2*(curr_index):2*(curr_index)+2] = f_block

        def forward_matrices(curr_index,next_index,nextnext_index):

            # Using central differencing for zeta = 2:NI-1
            x0_zeta = 0.5*(-3*r0[2*(curr_index)] + 4*r0[2*(next_index)] - r0[2*(nextnext_index)])
            y0_zeta = 0.5*(-3*r0[2*(curr_index)+1] + 4*r0[2*(next_index)+1] - r0[2*(nextnext_index)+1])

            # Computing gamma
            gamma = x0_zeta*x0_zeta + y0_zeta*y0_zeta

            # Use forward differencing for eta
            x0_eta = -dA[curr_index]/gamma*y0_zeta
            y0_eta = dA[curr_index]/gamma*x0_zeta

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[2*(next_index):2*(next_index)+2]-rm1[2*(index):2*(index)+2]) + norm(rm1[2*(nextnext_index):2*(nextnext_index)+2]-rm1[2*(index):2*(index)+2])
            dden = norm(r0[2*(next_index):2*(next_index)+2]-r0[2*(index):2*(index)+2]) + norm(r0[2*(nextnext_index):2*(nextnext_index)+2]-r0[2*(index):2*(index)+2])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

            # Compute other auxiliary variables
            alpha = x0_zeta*x0_eta - y0_zeta*y0_eta
            beta = x0_zeta*y0_eta + x0_eta*y0_zeta

            # Compute C = Binv*A
            C = array([[alpha, beta],[beta, -alpha]])/gamma

            # Compute smoothing coefficients
            epsE, epsI = dissipation_coefficients(layer_index, gamma, x0_eta, y0_eta, d_sensor, angle)

            # Compute RHS components
            Binvg = array([x0_eta, y0_eta])
            De = epsE*(r0[2*(nextnext_index):2*(nextnext_index)+2] - 2*r0[2*(next_index):2*(next_index)+2] + r0[2*(curr_index):2*(curr_index)+2])

            # Compute block matrices
            L_block = -0.5*(1+theta)*C - epsI*eye(2)
            M_block = 2*(1+theta)*C + 2*epsI*eye(2)
            N_block = -1.5*(1+theta)*C + (1-epsI)*eye(2)
            f_block = Binvg + De

            # Populate matrix
            K[2*(curr_index):2*(curr_index)+2,2*(nextnext_index):2*(nextnext_index)+2] = L_block
            K[2*(curr_index):2*(curr_index)+2,2*(next_index):2*(next_index)+2] = M_block
            K[2*(curr_index):2*(curr_index)+2,2*(curr_index):2*(curr_index)+2] = N_block
            f[2*(curr_index):2*(curr_index)+2] = f_block

        def backward_matrices(curr_index,prev_index,prevprev_index):

            # Using central differencing for zeta = 2:NI-1
            x0_zeta = 0.5*(3*r0[2*(curr_index)] - 4*r0[2*(prev_index)] + r0[2*(prevprev_index)])
            y0_zeta = 0.5*(3*r0[2*(curr_index)+1] - 4*r0[2*(prev_index)+1] + r0[2*(prevprev_index)+1])

            # Computing gamma
            gamma = x0_zeta*x0_zeta + y0_zeta*y0_zeta

            # Use forward differencing for eta
            x0_eta = -dA[curr_index]/gamma*y0_zeta
            y0_eta = dA[curr_index]/gamma*x0_zeta

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[2*(prev_index):2*(prev_index)+2]-rm1[2*(index):2*(index)+2]) + norm(rm1[2*(prevprev_index):2*(prevprev_index)+2]-rm1[2*(index):2*(index)+2])
            dden = norm(r0[2*(prev_index):2*(prev_index)+2]-r0[2*(index):2*(index)+2]) + norm(r0[2*(prevprev_index):2*(prevprev_index)+2]-r0[2*(index):2*(index)+2])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

            # Compute other auxiliary variables
            alpha = x0_zeta*x0_eta - y0_zeta*y0_eta
            beta = x0_zeta*y0_eta + x0_eta*y0_zeta

            # Compute C = Binv*A
            C = array([[alpha, beta],[beta, -alpha]])/gamma

            # Compute smoothing coefficients
            epsE, epsI = dissipation_coefficients(layer_index, gamma, x0_eta, y0_eta, d_sensor, angle)

            # Compute RHS components
            Binvg = array([x0_eta, y0_eta])
            De = epsE*(r0[2*(prevprev_index):2*(prevprev_index)+2] - 2*r0[2*(prev_index):2*(prev_index)+2] + r0[2*(curr_index):2*(curr_index)+2])

            # Compute block matrices
            L_block = 0.5*(1+theta)*C - epsI*eye(2)
            M_block = -2*(1+theta)*C + 2*epsI*eye(2)
            N_block = 1.5*(1+theta)*C + (1-epsI)*eye(2)
            f_block = Binvg + De

            # Populate matrix
            K[2*(curr_index):2*(curr_index)+2,2*(prevprev_index):2*(prevprev_index)+2] = L_block
            K[2*(curr_index):2*(curr_index)+2,2*(prev_index):2*(prev_index)+2] = M_block
            K[2*(curr_index):2*(curr_index)+2,2*(curr_index):2*(curr_index)+2] = N_block
            f[2*(curr_index):2*(curr_index)+2] = f_block

        # Now loop over each node

        for index in [0]:

            if bc1 is 'free':

                # Name indexes
                curr_index = index
                next_index = index+1
                nextnext_index = index+2

                # Call assembly routine
                forward_matrices(curr_index,next_index,nextnext_index)

            elif bc1 is 'continuous':

                # Name indexes
                prev_index = NI-2
                curr_index = index
                next_index = index+1

                # Call assembly routine
                central_matrices(prev_index,curr_index,next_index)

            elif bc1 is 'splay':

                # Get coordinates
                r_curr = r0[2*(index):2*(index)+2]
                r_next = r0[2*(index+1):2*(index+1)+2]

                # Get vector that connects r_next to r_curr
                d_vec = r_curr - r_next
                d_vec_rot = rotMatrix(-pi/2).dot(d_vec)

                # Populate matrix
                K[2*index:2*index+2,2*index:2*index+2] = array([[d_vec[0], d_vec[1]],[d_vec_rot[0], d_vec_rot[1]]])
                f[2*index:2*index+2] = array([0, dA[index]*(1-sigmaSplay)])
                #f[2*index:2*index+2] = array([0, s*(1-sigmaSplay)*norm(d_vec_rot)])

            elif bc1 is 'constX':

                # Populate matrix
                K[2*index:2*index+2,2*(index+1):2*(index+1)+2] = [[0, 0],[0,-1]]
                K[2*index:2*index+2,2*index:2*index+2] = eye(2)
                f[2*index:2*index+2] =  [0, 0]

            elif bc1 is 'constY':

                # Populate matrix
                K[2*index:2*index+2,2*(index+1):2*(index+1)+2] = [[-1,0],[0, 0]]
                K[2*index:2*index+2,2*index:2*index+2] = eye(2)
                f[2*index:2*index+2] =  [0, 0]

        for index in range(1,NI-1):

            # Name indexes
            prev_index = index-1
            curr_index = index
            next_index = index+1

            # Call assembly routine
            central_matrices(prev_index,curr_index,next_index)

        for index in [NI-1]:

            if bc2 is 'free':

                # Name indexes
                curr_index = index
                prev_index = index-1
                prevprev_index = index-2

                # Call assembly routine
                backward_matrices(curr_index,prev_index,prevprev_index)

            elif bc2 is 'continuous':

                # Populate matrix (use same displacements of first node)
                K[2*index:2*index+2,2*index:2*index+2] = eye(2)
                K[2*index:2*index+2,:2] = -eye(2)
                f[2*index:2*index+2] = [0, 0]

            elif bc2 is 'splay':

                # Get coordinates
                r_curr = r0[2*(index):2*(index)+2]
                r_prev = r0[2*(index-1):2*(index-1)+2]

                # Get vector that connects r_next to r_curr
                d_vec = r_curr - r_prev
                d_vec_rot = rotMatrix(pi/2).dot(d_vec)

                # Populate matrix
                K[2*index:2*index+2,2*index:2*index+2] = array([[d_vec[0], d_vec[1]],[d_vec_rot[0], d_vec_rot[1]]])
                f[2*index:2*index+2] = array([0, dA[index]*(1-sigmaSplay)])
                #f[2*index:2*index+2] = array([0, s*(1-sigmaSplay)*norm(d_vec_rot)])

            elif bc2 is 'constX':

                # Populate matrix
                K[2*index:2*index+2,2*(index-1):2*(index-1)+2] = [[0, 0],[0,-1]]
                K[2*index:2*index+2,2*index:2*index+2] = eye(2)
                f[2*index:2*index+2] =  [0, 0]

            elif bc2 is 'constY':

                # Populate matrix
                K[2*index:2*index+2,2*(index-1):2*(index-1)+2] = [[-1, 0],[0, 0]]
                K[2*index:2*index+2,2*index:2*index+2] = eye(2)
                f[2*index:2*index+2] =  [0, 0]

        # RETURNS
        return K,f

    def smoothing(r,eta):

        # This function does the grid smoothing

        # IMPORTS
        from numpy import copy, zeros
        from numpy.linalg import norm

        # Find number of nodes
        NI = int(len(r)/2)

        # Loop over the desired number of smoothing passes
        for index_pass in range(num_smoothing_passes):

            # Initialize array of smoothed coordinates
            r_smooth = zeros(2*NI)

            # Copy the edge nodes
            r_smooth[:2] = r[:2]
            r_smooth[-2:] = r[-2:]

            # Smooth every node
            for index in range(1,NI-1):

                # Get coordinates
                r_curr = r[2*(index):2*(index)+2]
                r_next = r[2*(index+1):2*(index+1)+2]
                r_prev = r[2*(index-1):2*(index-1)+2]

                # Compute distances
                lp = norm(r_next - r_curr)
                lm = norm(r_curr - r_prev)

                # Compute alpha'
                alphaP = min(alphaP0, alphaP0*(eta-2)/NJ)

                # Compute smoothed coordinates
                r_smooth[2*index:2*index+2] = (1-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)

            # Copy coordinates to allow next pass
            r = copy(r_smooth)

        # RETURNS
        return r

    def dissipation_coefficients(layer_index, gamma, x0_eta, y0_eta, d_sensor, angle):

        # IMPORTS
        from numpy import sqrt, cos, pi

        # Compute N (Eq. 6.3)
        N = sqrt((x0_eta*x0_eta + y0_eta*y0_eta)/gamma)

        # Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layer_index+2
        ltrans = int(3/4*NJ)
        if l <= ltrans:
            Sl = sqrt((l-1)/(ltrans-1))
        else:
            Sl = sqrt((ltrans-1)/(NJ-1))

        # Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max([d_sensor**(2/Sl), 0.1])

        # Compute a (Eq 6.12 adjusted for entire angle (angle=2*alpha))
        if angle >= pi: # Convex corner
            a = 1.0
        else:
            a = 1.0/(1.0 - cos(angle/2)*cos(angle/2))

        # Compute auxiliary variable R (Eq. 6.4)
        R = Sl*dbar*a

        # Compute the dissipation coefficients
        epsE = epsE0*R*N
        epsI = 2*epsE

        # RETURNS
        return epsE, epsI

    #===========================================================

    '''
    The marching function actually begins here
    '''

    # MARCH!!!
    for layer_index in range(NJ-1):
        # Get the previous coordinates
        r_prev = R[layer_index,:]

        # Use it as a guess to the next iteration
        r_next = copy(r_prev)

        # Compute the total marched distance
        s_tot = s_tot + s

        # Subiter
        for index_sub in range(1):
            r_next = sub_iteration(r_prev, r_prevprev, s, s_tot, layer_index)

        # Store grid points
        R[layer_index+1,:] = r_next

        # Update step size
        s = s*growth

        # Update r_prevprev
        r_prevprev = R[layer_index,:]

    # Convert to X and Y
    X = R[:,::2]
    Y = R[:,1::2]

    # RETURNS
    return X,Y

'''
==============================================
MORE AUXILIARY FUNCTIONS
==============================================
'''

#=============================================
#=============================================

def give_angle(r0,r1,r2):

    '''
    This function gives the angle between the vectors joining
    r0 to r1 and r1 to r2. We assume that the body is to the right
    of the vector, while the mesh is propagated to the left

          r0
          |
     body | mesh
          |
          V angle
          r1--------->r2
              body

    angles > pi indicate convex corners, while angles < pi
    indicate concave corners
    '''

    # IMPORTS
    from numpy import arctan2, pi

    dx1 = r1 - r0
    dx2 = r2 - r1

    dot = dx1.dot(dx2) #dot product
    det = dx1[0]*dx2[1] - dx1[1]*dx2[0]      # determinant
    angle = pi - arctan2(det, dot)

    return angle

#=============================================
#=============================================

def findRatio(s_max, s0, NJ, ratioGuess):

    '''
    This function returns the geometrical progression ratio that satisfies
    the farfield distance and the number of cells. Newton search is used
    INPUTS
    s_max: distance that should be reached
    s0: cell edge length at the wall
    NJ: number of cells used to reach farfield
    '''

    # Extra parameters
    nIters = 400 # Maximum number of iterations for Newton search
    q0 = ratioGuess # Initial ratio

    # Initialize ratio
    q = q0

    # Newton search loop
    for i in range(nIters):
       # Residual function
       R = s0*(1-q**(NJ-1)) - s_max*(1-q)

       # Residual derivative
       Rdot = -(NJ-1)*s0*q**(NJ-2) + s_max

       # Update ratio with Newton search
       q = q - R/Rdot

    # Check if we got a reasonable value
    if (q <= 1) or (q >= q0):
       print('Ratio may be too large...')
       print('Increase number of cells or reduce extension')
       from sys import exit
       exit()

    # RETURNS
    return q

#=============================================
#=============================================

def findRadius(r):

    '''
    This function find the largest radius of the bounding box
    that encompasses the given mesh. This length will be used
    to calculate the maching distance
    '''

    # IMPORTS
    from numpy import max, min

    # Split coordinates
    x = r[::2]
    y = r[1::2]

    # Find bounds
    minX = min(x)
    maxX = max(x)
    minY = min(y)
    maxY = max(y)

    # Find longest radius (we give only half of the largest side to be considered as radius)
    radius = max([maxX-minX, maxY-minY])/2

    # RETURNS
    return radius

#=============================================
#=============================================

def rotMatrix(theta):

    '''
    This function gives the 2D rotation matrix for a counter clock-wise
    rotation of theta radians
    '''

    # IMPORTS
    from numpy import array, sin, cos

    # RETURNS
    return array([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])

#=============================================
#=============================================

def plot_grid(X,Y,show=False):

    # This function plots the grid

    # IMPORTS
    import matplotlib.pyplot as plt

    # Initialize figure
    fig = plt.figure()

    # Plot zeta lines
    plt.plot(X,Y,'k')

    # Plot eta lines
    plt.plot(X.T,Y.T,'k')

    # Plot the initial surface
    plt.plot(X[0,:],Y[0,:],'r',linewidth=2)

    # Add labels
    #plt.xlabel('x')
    #plt.ylabel('y')
    plt.tick_params(axis='x',
                    which='both',
                    top=False,
                    bottom=False,
                    labelbottom=False)
    plt.tick_params(axis='y',
                    which='both',
                    left=False,
                    right=False,
                    labelleft=False)

    # Adjust axis
    plt.axis('equal')
    #plt.xlim([-0.5, 1.5])
    #plt.ylim([-0.05, 1.35])

    if show:
        # Show
        plt.show()

    # Return the figure handle
    return fig

#=============================================
#=============================================

def export_plot3d(X,Y,filename,zSpan=1,zNodes=2):

    '''
    This function exports a 3D mesh in plot3d format.
    The user specifies the span and number of nodes in the z direction.
    '''

    # IMPORTS
    from .plot3d_writer import Grid, export_plot3d
    from numpy import array, copy, ones, linspace

    # Initialize grid object
    myGrid = Grid()

    # Expand the coordinate matrices
    X3d = array([copy(X) for dummy in range(zNodes)])
    Y3d = array([copy(Y) for dummy in range(zNodes)])

    # Create the span-wise coordinates
    Z3d = array([z*ones(X.shape) for z in linspace(0,zSpan,zNodes)])

    # Add block to the grid
    myGrid.add_block(X3d, Y3d, Z3d)

    # Export grid
    export_plot3d(myGrid, filename)

#=============================================
#=============================================

def export_pickle(X,Y,filename):

    import pickle

    with open(filename,'wb') as fid:
        pickle.dump([X,Y],fid)
