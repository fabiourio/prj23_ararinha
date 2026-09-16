# GENERAL IMPORTS
import numpy as np
from .constants import gravity
from .balance import tank_properties

#========================================

def moment_of_inertia(airplane, fuel_frac=1.0, payload_frac=1.0):

    # Unpack dictionary
    #W0 = airplane['W0']
    W_payload = airplane['inputs']['W_payload']
    xcg_payload = airplane['inputs']['xcg_payload']
    W_crew = airplane['inputs']['W_crew']
    xcg_crew = airplane['inputs']['xcg_crew']
    W_empty = airplane['thrust_matching']['W_empty']
    xcg_empty = airplane['empty_weight']['xcg_empty']
    W_fuel = airplane['thrust_matching']['W_fuel']
    c_tank_c_w = airplane['inputs']['c_tank_c_w']
    x_tank_c_w = airplane['inputs']['x_tank_c_w']
    b_tank_b_w_start = airplane['inputs']['b_tank_b_w_start']
    b_tank_b_w_end = airplane['inputs']['b_tank_b_w_end']
    sweep_w = airplane['inputs']['sweep_w']
    b_w = airplane['geometry']['b_w']
    xr_w = airplane['inputs']['xr_w']
    cr_w = airplane['geometry']['cr_w']
    ct_w = airplane['geometry']['ct_w']
    tcr_w = airplane['inputs']['tcr_w']
    tct_w = airplane['inputs']['tct_w']
    rho_fuel = airplane['inputs']['rho_fuel']

    # CG position
    V_maxfuel, W_maxfuel, xcg_fuel, ycg_fuel = tank_properties(cr_w, ct_w, tcr_w, tct_w, b_w, sweep_w, xr_w,
                                                               x_tank_c_w, c_tank_c_w, b_tank_b_w_start, b_tank_b_w_end,
                                                               rho_fuel, gravity)
    
    xcg = (W_empty*xcg_empty + fuel_frac*W_fuel*xcg_fuel + payload_frac*W_payload*xcg_payload + W_crew*xcg_crew)/(W_empty + fuel_frac*W_fuel + payload_frac*W_payload + W_crew)
    
    Ixxt = 0.0
    Iyyt = 0.0
    Izzt = 0.0
    Ixyt = 0.0
    Ixzt = 0.0
    Iyzt = 0.0
    
    # Reference point
    xref = xcg
    yref = 0.0
    zref = 0.0
    
    # Right wing + fuel (SIMPLIFICATION: fuel stored in the wing)    
    xr = airplane['inputs']['xr_w']
    yr = 0.0
    zr = airplane['inputs']['zr_w']
    cr = airplane['geometry']['cr_w']
    tcr = airplane['inputs']['tcr_w']
    xt = airplane['geometry']['xt_w']
    yt = airplane['geometry']['yt_w']
    zt = airplane['geometry']['zt_w']
    ct = airplane['geometry']['ct_w']
    tct = airplane['inputs']['tct_w']
    mass = (airplane['empty_weight']['W_w']+fuel_frac*W_fuel)*0.5/gravity
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = surf_moi(xr,yr,zr,cr,tcr,
                                       xt,yt,zt,ct,tct,
                                       xref,yref,zref,
                                       mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Left wing      
    yt = -airplane['geometry']['yt_w']
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = surf_moi(xr,yr,zr,cr,tcr,
                                       xt,yt,zt,ct,tct,
                                       xref,yref,zref,
                                       mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Right HT                
    xr = airplane['geometry']['xr_h']
    yr = 0.0
    zr = airplane['inputs']['zr_h']
    cr = airplane['geometry']['cr_h']
    tcr = airplane['inputs']['tcr_h']
    xt = airplane['geometry']['xt_h']
    yt = airplane['geometry']['yt_h']
    zt = airplane['geometry']['zt_h']
    ct = airplane['geometry']['ct_h']
    tct = airplane['inputs']['tct_h']
    mass = airplane['empty_weight']['W_h']*0.5/gravity
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = surf_moi(xr,yr,zr,cr,tcr,
                                       xt,yt,zt,ct,tct,
                                       xref,yref,zref,
                                       mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Left HT
    yt = -airplane['geometry']['yt_h']
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = surf_moi(xr,yr,zr,cr,tcr,
                                       xt,yt,zt,ct,tct,
                                       xref,yref,zref,
                                       mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # VT
    xr = airplane['geometry']['xr_v']
    yr = 0.0
    zr = airplane['inputs']['zr_v']
    cr = airplane['geometry']['cr_v']
    tcr = airplane['inputs']['tcr_v']
    xt = airplane['geometry']['xt_v']
    yt = 0.0
    zt = airplane['geometry']['zt_v']
    ct = airplane['geometry']['ct_v']
    tct = airplane['inputs']['tct_v']
    mass = airplane['empty_weight']['W_v']/gravity
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = surf_moi(xr,yr,zr,cr,tcr,
                                       xt,yt,zt,ct,tct,
                                       xref,yref,zref,
                                       mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Fuselage + allelse + payload + crew
    xn = 0.0
    yn = 0.0
    zn = 0.0
    Lcyl = airplane['inputs']['L_f']
    Dcyl = airplane['inputs']['D_f']
    mass = (airplane['empty_weight']['W_f']+airplane['empty_weight']['W_allelse'] + \
            payload_frac*W_payload + W_crew)/gravity
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = cyl_moi(xn,yn,zn,
                                      Lcyl,Dcyl,
                                      xref,yref,zref,
                                      mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Right engine
    xn = airplane['inputs']['x_n']
    yn = airplane['inputs']['y_n']
    zn = airplane['inputs']['z_n']
    Lcyl = airplane['inputs']['L_n']
    Dcyl = airplane['inputs']['D_n']
    mass = airplane['empty_weight']['W_eng']*0.5/gravity
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = cyl_moi(xn,yn,zn,
                                      Lcyl,Dcyl,
                                      xref,yref,zref,
                                      mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Left engine
    yn = -airplane['inputs']['y_n']
    
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz = cyl_moi(xn,yn,zn,
                                      Lcyl,Dcyl,
                                      xref,yref,zref,
                                      mass)
    
    Ixxt += Ixx
    Iyyt += Iyy
    Izzt += Izz
    Ixyt += Ixy
    Ixzt += Ixz
    Iyzt += Iyz
    
    # Nose landing gear
    if airplane['inputs']['x_nlg'] is not None:
        
        x = airplane['inputs']['x_nlg']
        y = 0.0
        z = airplane['inputs']['z_lg']
        mass = airplane['empty_weight']['W_nlg']/gravity
    
        Ixx,Iyy,Izz,Ixy,Ixz,Iyz = point_moi(x,y,z,
                                            xref,yref,zref,
                                            mass)
    
        Ixxt += Ixx
        Iyyt += Iyy
        Izzt += Izz
        Ixyt += Ixy
        Ixzt += Ixz
        Iyzt += Iyz
    
    # Right landing gear
    if airplane['inputs']['x_mlg'] is not None:
        
        x = airplane['inputs']['x_mlg']
        y = airplane['inputs']['y_mlg']
        z = airplane['inputs']['z_lg']
        mass = airplane['empty_weight']['W_mlg']*0.5/gravity
    
        Ixx,Iyy,Izz,Ixy,Ixz,Iyz = point_moi(x,y,z,
                                            xref,yref,zref,
                                            mass)
    
        Ixxt += Ixx
        Iyyt += Iyy
        Izzt += Izz
        Ixyt += Ixy
        Ixzt += Ixz
        Iyzt += Iyz
    
        # Left landing gear
        y = -airplane['inputs']['y_mlg']
    
        Ixx,Iyy,Izz,Ixy,Ixz,Iyz = point_moi(x,y,z,
                                            xref,yref,zref,
                                            mass)
    
        Ixxt += Ixx
        Iyyt += Iyy
        Izzt += Izz
        Ixyt += Ixy
        Ixzt += Ixz
        Iyzt += Iyz

    airplane['moment_of_inertia'] = {}
    airplane['moment_of_inertia']['Ixx'] = Ixxt
    airplane['moment_of_inertia']['Iyy'] = Iyyt
    airplane['moment_of_inertia']['Izz'] = Izzt
    airplane['moment_of_inertia']['Ixy'] = Ixyt
    airplane['moment_of_inertia']['Ixz'] = Ixzt
    airplane['moment_of_inertia']['Iyz'] = Iyzt
    
    return None

#----------------------------------------

def surf_moi(xr,yr,zr,cr,tcr,
             xt,yt,zt,ct,tct,
             xref,yref,zref,
             mass,
             Nc=100,Nb=200):
    '''
    This function computes the moment of inertia
    of a half-surface.

    Returns
    -------
    None.

    '''
    
    # Determine leading and trailing edge points
    xrLE = xr
    yrLE = yr
    zrLE = zr
    xrTE = xr+cr
    yrTE = yr
    zrTE = zr
    
    xtLE = xt
    ytLE = yt
    ztLE = zt
    xtTE = xt+ct
    ytTE = yt
    ztTE = zt
    
    # Inverse matrix for bilinear interpolation
    Ainv = np.linalg.inv(np.array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,0,1,0],
                                   [1,1,1,1]]))
    
    # Create a mesh for properties varying along the wing.
    # We do this by creating a linear interpolation of properties.
    # uu=0 at the LE and uu=1 at TE
    # vv=0 at the root and vv=1 at the tip
    uu = np.linspace(0,1,Nc)
    vv = np.linspace(0,1,Nb)
    
    UU, VV = np.meshgrid(uu, vv)
    
    # X interpolation
    iw = Ainv.dot([xrLE,xrTE,xtLE,xtTE])
    XX = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Y interpolation
    iw = Ainv.dot([yrLE,yrTE,ytLE,ytTE])
    YY = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Z interpolation
    iw = Ainv.dot([zrLE,zrTE,ztLE,ztTE])
    ZZ = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # absolute thickness interpolation
    iw = Ainv.dot([tcr*cr,tcr*cr,tct*ct,tct*ct])
    TT = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Compute areas of each quadrilateral using cross product
    # of diagonals
    dX1 = XX[1:,1:]  - XX[:-1,:-1]
    dX2 = XX[1:,:-1] - XX[:-1,1:]
    dY1 = YY[1:,1:]  - YY[:-1,:-1]
    dY2 = YY[1:,:-1] - YY[:-1,1:]
    dZ1 = ZZ[1:,1:]  - ZZ[:-1,:-1]
    dZ2 = ZZ[1:,:-1] - ZZ[:-1,1:]
    
    # Cross product and area
    dXv = dY1*dZ2 - dY2*dZ1
    dYv = dZ1*dX2 - dX1*dZ2
    dZv = dX1*dY2 - dY1*dX2
    
    AA = 0.5*np.sqrt(dXv**2 + dYv**2 + dZv**2)
    
    # Get average thickness and centroid of each quadrilateral
    Xc = 0.25*(XX[:-1,:-1] + XX[1:,:-1] + XX[1:,1:] + XX[:-1,1:])
    Yc = 0.25*(YY[:-1,:-1] + YY[1:,:-1] + YY[1:,1:] + YY[:-1,1:])
    Zc = 0.25*(ZZ[:-1,:-1] + ZZ[1:,:-1] + ZZ[1:,1:] + ZZ[:-1,1:])
    Tc = 0.25*(TT[:-1,:-1] + TT[1:,:-1] + TT[1:,1:] + TT[:-1,1:])
    
    # Compute volume of each element and of the wing
    Vc = AA*Tc
    vol = np.sum(Vc)
    
    # Wing density
    rho = mass/vol
    
    # Moments of inertia
    Ixx = np.sum(rho*Vc*((Yc-yref)**2 + (Zc-zref)**2))
    Iyy = np.sum(rho*Vc*((Xc-xref)**2 + (Zc-zref)**2))
    Izz = np.sum(rho*Vc*((Xc-xref)**2 + (Yc-yref)**2))
    Ixy = np.sum(rho*Vc*(Xc-xref)*(Yc-yref))
    Ixz = np.sum(rho*Vc*(Xc-xref)*(Zc-zref))
    Iyz = np.sum(rho*Vc*(Yc-yref)*(Zc-zref))
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz

#----------------------------------------

def cyl_moi(xn,yn,zn,
            Lcyl,Dcyl,
            xref,yref,zref,
            mass,
            Nr=200,Nt=200,Nl=10):
    '''
    This function computes the moment of inertia
    of a half-surface.

    Returns
    -------
    None.

    '''
    
    # Compute volume of the cylinder
    vol = np.pi*Dcyl**2/4*Lcyl
    rho = mass/vol
    
    # Compute length of each disk
    tdisk = Lcyl/Nl
    
    # First we get the distribution of one disk
    # perpendicular to the x axis
            
    # Create a mesh for properties varying along the disk.
    # We do this by creating a linear interpolation of properties.
    # uu=0 at the center and uu=1 at the rim
    # vv=0 at for theta=0 and vv=1 at for theta=2*pi, starting at the left
    uu = np.linspace(0,Dcyl/2,Nr)
    vv = np.linspace(0,2*np.pi,Nt)
    
    UU, VV = np.meshgrid(uu, vv)
    
    YY = UU*np.cos(VV) + yn
    ZZ = UU*np.sin(VV) + zn
    XX = np.zeros_like(ZZ)
    
    # Compute areas of each quadrilateral using cross product
    # of diagonals
    dX1 = XX[1:,1:]  - XX[:-1,:-1]
    dX2 = XX[1:,:-1] - XX[:-1,1:]
    dY1 = YY[1:,1:]  - YY[:-1,:-1]
    dY2 = YY[1:,:-1] - YY[:-1,1:]
    dZ1 = ZZ[1:,1:]  - ZZ[:-1,:-1]
    dZ2 = ZZ[1:,:-1] - ZZ[:-1,1:]
    
    # Cross product and area
    dXv = dY1*dZ2 - dY2*dZ1
    dYv = dZ1*dX2 - dX1*dZ2
    dZv = dX1*dY2 - dY1*dX2
    
    AA = 0.5*np.sqrt(dXv**2 + dYv**2 + dZv**2)
    
    # Get average thickness and centroid of each quadrilateral
    Xc = 0.25*(XX[:-1,:-1] + XX[1:,:-1] + XX[1:,1:] + XX[:-1,1:])
    Yc = 0.25*(YY[:-1,:-1] + YY[1:,:-1] + YY[1:,1:] + YY[:-1,1:])
    Zc = 0.25*(ZZ[:-1,:-1] + ZZ[1:,:-1] + ZZ[1:,1:] + ZZ[:-1,1:])
    
    # Compute volume of each element and of the disk
    Vc = AA*tdisk
    
    # Moments of inertia of each disk
    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    Ixy = 0.0
    Ixz = 0.0
    Iyz = 0.0
    for ii in range(Nl):
        
        # Ajust X coordinate
        Xc[:] = xn + tdisk/2 + ii*tdisk
        
        Ixx = Ixx + np.sum(rho*Vc*((Yc-yref)**2 + (Zc-zref)**2))
        Iyy = Iyy + np.sum(rho*Vc*((Xc-xref)**2 + (Zc-zref)**2))
        Izz = Izz + np.sum(rho*Vc*((Xc-xref)**2 + (Yc-yref)**2))
        Ixy = Ixy + np.sum(rho*Vc*(Xc-xref)*(Yc-yref))
        Ixz = Ixz + np.sum(rho*Vc*(Xc-xref)*(Zc-zref))
        Iyz = Iyz + np.sum(rho*Vc*(Yc-yref)*(Zc-zref))
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz

#----------------------------------------
    
def point_moi(x,y,z,
              xref,yref,zref,
              mass):
    '''
    This function computes the moment of inertia
    of a point mass.

    Returns
    -------
    None.

    '''
    
    # Moments of inertia
    Ixx = mass*((y-yref)**2 + (z-zref)**2)
    Iyy = mass*((x-xref)**2 + (z-zref)**2)
    Izz = mass*((x-xref)**2 + (y-yref)**2)
    Ixy = mass*(x-xref)*(y-yref)
    Ixz = mass*(x-xref)*(z-zref)
    Iyz = mass*(y-yref)*(z-zref)
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz

#----------------------------------------