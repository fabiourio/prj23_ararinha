'''
This module contains grid warping subroutines
'''

import airfoil_mod as am
import numpy as np

#=====================================

def warp_mesh(x, y, dx, dy, X, Y):

    Xbase = np.vstack([x,y]).T
    Ubase = np.vstack([dx,dy]).T
    Xsample = np.vstack([X.flatten(), Y.flatten()]).T

    Usample = idwInterp(Xbase, Ubase, Xsample, p=2.5, distTol=1e-9)

    dX = Usample[:,0].reshape(X.shape)
    dY = Usample[:,1].reshape(Y.shape)

    return X+dX, Y+dY

#=====================================

def cst_deriv(Au, Al, dAu, dAl, x, X, Y, adj_funcs, grads):

    def cst_plus_warp(dAu, dAl):
        '''
        Define function that concatenates surface generation and mesh warping
        to propagate derivatives.
        '''

        airfoil = am.cstfoil(Au, Al, x)

        xf = airfoil['x_coord']
        yf = airfoil['y_coord']
        maxt = airfoil['max_thickness']
        xmaxt = airfoil['x_max_thickness']
        mint = airfoil['min_thickness']
        maxc = airfoil['max_camber']
        xmaxc = airfoil['x_max_camber']

        airfoil = am.cstfoil(Au+dAu, Al+dAl, x)

        xfw = airfoil['x_coord']
        yfw = airfoil['y_coord']
        maxtw = airfoil['max_thickness']
        xmaxtw = airfoil['x_max_thickness']
        mintw = airfoil['min_thickness']
        maxcw = airfoil['max_camber']
        xmaxcw = airfoil['x_max_camber']

        dxf = xfw-xf
        dyf = yfw-yf

        Xbase = np.vstack([xf,yf]).T
        Ubase = np.vstack([dxf,dyf]).T
        Xsample = np.vstack([X.flatten(), Y.flatten()]).T

        Usample = idwInterp(Xbase, Ubase, Xsample, p=2.5, distTol=1e-9, use_complex=True)

        dX = Usample[:,0].reshape(X.shape)
        dY = Usample[:,1].reshape(Y.shape)

        return X+dX, Y+dY, maxtw, xmaxtw, mintw, maxcw, xmaxcw

    # Apply complex step
    step = 1e-30

    #gradAl = []
    dXdAl = []
    dYdAl = []
    dfdAl = []
    dmaxtdAl = []
    #dxmaxtdAl = [] # Not implemented because function is discontinuous
    dmintdAl = []
    #dmaxcdAl = []
    #dxmaxcdAl = []
    

    for ii in range(len(Al)):

        delta = np.array(np.zeros_like(Al), dtype=complex)

        delta[ii] = 1j*step

        Xi, Yi, maxti, xmaxti, minti, maxci, xmaxci = cst_plus_warp(dAu, dAl+delta)

        Xd_CS = np.imag(Xi)/step
        Yd_CS = np.imag(Yi)/step
        maxt_CS = np.imag(maxti)/step
        #xmaxt_CS = np.imag(xmaxti)/step
        mint_CS = np.imag(minti)/step
        #maxc_CS = np.imag(maxci)/step
        #xmaxc_CS = np.imag(xmaxci)/step

        #dfdAlii = np.sum(Xd_CS*grads[adj_func]['X'] + Yd_CS*grads[adj_func]['Y'])

        #dfdAl.append(dfdAlii)

        dXdAl.append(Xd_CS)
        dYdAl.append(Yd_CS)
        dmaxtdAl.append(maxt_CS)
        #dxmaxtdAl.append(xmaxt_CS)
        dmintdAl.append(mint_CS)
        #dmaxcdAl.append(maxc_CS)
        #dxmaxcdAl.append(xmaxc_CS)

    dXdAu = []
    dYdAu = []
    dmaxtdAu = []
    #dxmaxtdAu = []
    dmintdAu = []
    #dmaxcdAu = []
    #dxmaxcdAu = []

    for ii in range(len(Au)):

        delta = np.array(np.zeros_like(Au), dtype=complex)

        delta[ii] = 1j*step

        Xi, Yi, maxti, xmaxti, minti, maxci, xmaxci = cst_plus_warp(dAu+delta, dAl)

        Xd_CS = np.imag(Xi)/step
        Yd_CS = np.imag(Yi)/step
        maxt_CS = np.imag(maxti)/step
        #xmaxt_CS = np.imag(xmaxti)/step
        mint_CS = np.imag(minti)/step
        #maxc_CS = np.imag(maxci)/step
        #xmaxc_CS = np.imag(xmaxci)/step

        #dfdAuii = np.sum(Xd_CS*grads[adj_func]['X'] + Yd_CS*grads[adj_func]['Y'])

        #dfdAu.append(dfdAuii)

        dXdAu.append(Xd_CS)
        dYdAu.append(Yd_CS)
        dmaxtdAu.append(maxt_CS)
        #dxmaxtdAu.append(xmaxt_CS)
        dmintdAu.append(mint_CS)
        #dmaxcdAu.append(maxc_CS)
        #dxmaxcdAu.append(xmaxc_CS)

    # Gradients of the CFD function
    for adj_func in grads.keys():
    
        # Chain rule to take derivatives from mesh nodes to airfoil perturbations
        dfdAl = []
        for ii in range(len(Al)):
            dfdAl_ii = np.sum(dXdAl[ii]*grads[adj_func]['X'] + dYdAl[ii]*grads[adj_func]['Y'])
            dfdAl.append(dfdAl_ii)

        dfdAu = []
        for ii in range(len(Au)):
            dfdAu_ii = np.sum(dXdAu[ii]*grads[adj_func]['X'] + dYdAu[ii]*grads[adj_func]['Y'])
            dfdAu.append(dfdAu_ii)

        grads[adj_func]['dAl'] = np.array(dfdAl)
        grads[adj_func]['dAu'] = np.array(dfdAu)

    # Gradients of geometric parameters
    grad_geo = {'maxt':{'dAl':np.array(dmaxtdAl),
                        'dAu':np.array(dmaxtdAu)},
                #'xmaxt':{'dAl':np.array(dxmaxtdAl),
                #         'dAu':np.array(dxmaxtdAu)},
                'mint':{'dAl':np.array(dmintdAl),
                        'dAu':np.array(dmintdAu)},
                #'maxc':{'dAl':np.array(dmaxcdAl),
                #        'dAu':np.array(dmaxcdAu)},
                #'xmaxc':{'dAl':np.array(dxmaxcdAl),
                #         'dAu':np.array(dxmaxcdAu)}
                }
    grads.update(grad_geo)

    return grads


#====================================

def idwInterp(Xbase, Ubase, Xsample, p=2, distTol=1e-6, use_complex=False):

    '''
    INPUTS:

    Xbase: float[npts,ndim] -> Nodal coordinates of the baseline data.
           npts is the number of points and ndim is the dimension of the space.

    Ubase: float[npts,nfuncs] -> Function values at the baseline nodes.
           npts is the number of points and nfuncs is the number of field
           variables that will be interpolated. Each field variable is
           independently interpolated.

    Xsample: float[nsamples,ndim] -> Nodal coordinates of the sampling points.
             nsample is the number of sampling points and ndim is the
             dimension of the space.

    p: float -> IDW order. Higher order gives more emphasis to nearby points.

    distTol: float -> Tolerance to avoid singularities in the IDW interpolation.
             Nodes whose distance are smaller than distTol from a baseline node
             will get the full value from the baseline node.

    OUTPUTS:

    Usample: float[nsamples,nfuncs] -> Function values at the baseline nodes.
             nsample is the number of sampling points and nfuncs is the
             number of field variables that will be interpolated.
    '''

    # Get problem size
    npts = Xbase.shape[0]
    ndim = Xbase.shape[1]
    nfuncs = Ubase.shape[1]
    nsamples = Xsample.shape[0]

    ### BUILD WEIGHT MATRIX
    # This matrix W is of shape [nsample,npts]
    # Every row contains the weights used to interpolate the baseline nodes
    # to get the value for a sampling node.

    # Initialize weight matrix
    W = np.zeros((nsamples, npts))

    if use_complex:
        W = np.array(W, dtype=complex)

    # Loop for every sampling node
    for isample in range(nsamples):

        # Get coordinate of the sampling node
        Xref = Xsample[isample,:]

        # Compute distances of the sampling node to the baseline nodes
        # We also apply the IDW order here.
        dist = np.sqrt(np.sum((Xbase-Xref)**2, axis=1))

        # Find the closest node and minimum distance
        closestNodeID = np.argmin(dist)
        minDist = dist[closestNodeID]

        # If the sampling node is too close to a baseline node, we forget
        # about the averaging and just take the baseline value
        if minDist < distTol:

            W[isample, closestNodeID] = 1.0

        # Otherwise we do the standard IDW averaging
        else:

            # Compute the interpolation weights using the IDW order
            W[isample,:] = dist**(-p)

    # Create array with the sum of the weights
    Wsum = np.sum(W, axis=1)

    # Transform sum array into a column vector
    Wsum = np.array([Wsum]).T

    ### PERFORM IDW INTERPOLATION FOR EVERY FIELD VARIABLE

    Usample = W.dot(Ubase)/Wsum

    ### RETURNS
    return Usample