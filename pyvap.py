# pyvap : simple model of particle evaporation
# author: Adam Birdsall

from scipy.constants import pi, k
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from math import pow

def dn(t, y, Ma, rho, cinf, po, Dg, T, has_water=False, xh2o=None):
    '''Construct differential equations to describe change in n.
    
    Parameters:
    -----------
    t : float
    Time
    y: ndarray
    1D-array of number of molecules of each component.
    Ma : ndarray
    1D-array of molecular weights, kg molec^-1.
    rho : ndarray
    1D-array of number of densities, kg m^-3.
    cinf : ndarray
    1D-array of concentration of substance at infinite distance, molec m^-3.
    po : ndarray
    1D-array of pure compound vapor pressures, Pa.
    Dg : ndarray
    1D-array of gas-phase diffusivities, m^2 s^-1.
    T : float
    Temperature, K.
    has_water : boolean
    Include water in calculation.
    xh2o : float
    Fixed water mole fraction. Only used if has_water is True.
    
    Outputs
    -------
    dn : ndarray
    1D-array of dn for all components.
    '''
    ytot = y.sum()
    v = y*Ma/rho # array of partial volumes, m^3    
    if has_water:
        # known mole fraction of water
        ntot = 1/(1-xh2o) * ytot
        nh2o = xh2o * ntot
        vh2o = nh2o*Ma_h2o/rho_h2o
        vtot = v.sum()+vh2o # particle volume, m^3
        # print("vh2o: {:.1e}".format(vh2o))
    else:
        ntot = ytot
        vtot = v.sum()
        
    r = (3*vtot/(4*pi))**(1/3) # radius, m^3
        
    x = y/ntot # mole fractions, where ntot includes h2o if present
    # assume ideality in vapor pressure calculation
    cs = x*po/(k*T) # gas-phase concentration at surface, molec m^-3
    # array of differential changes in number of molecules
    dn = 4*pi*r*Dg*(cinf-cs)
    return dn

def evaporate(components, ninit, T, tsteps, dt, has_water=False, xh2o=None):
    '''calculate evaporation of multicomponent particle.
    
    num : int
    total number of integrated time points, including t=0'''
    
    # extract data from components and build data arrays
    Ma = np.empty(len(components))
    rho = np.empty(len(components))
    cinf = np.empty(len(components))
    po = np.empty(len(components))
    Dg = np.empty(len(components))
    for i, component in enumerate(components):
        Ma[i] = component['Ma']
        rho[i] = component['rho']
        cinf[i] = component['cinf']
        po[i] = calcp0(component['p0_a'], component['p0_b'], T)
        Dg[i] = component['Dg']
       
    # set up ode
    output = np.empty((int(num), len(components)))
    output[0, :] = ninit

    r = ode(dn)
    r.set_integrator('lsoda', with_jacobian=False,)
    r.set_initial_value(ninit, t=0)
    r.set_f_params(Ma, rho, cinf, po, Dg, T, has_water, xh2o)

    # integrate and save output
    entry = 0
    # use `entry` condition to avoid rounding errors causing
    # possible problems with `r.t` condition
    while r.successful() and entry < num-1:
        entry = int(round(r.t/dt))+1
        nextstep = r.integrate(r.t + dt)
        output[entry, :] = nextstep

    return output

def calcv(components, ns, has_water=False, xh2o=None):
    '''calculate volume of particle over time for given components'''
    vtot = np.zeros_like(ns.shape[0])
    for i, c in enumerate(components):
        v = ns[:,i]*c['Ma']/c['rho']
        vtot = vtot + v
    if has_water:
        # water not explicitly tracked, instead fixed
        # xh2o for fixed ambient RH (ah2o)
        nh2o = xh2o/(1-xh2o)*ns.sum(axis=1)
        vh2o = nh2o*Ma_h2o/rho_h2o
        vtot = vtot + vh2o
    return vtot

def calcr(components, ns, has_water=False, xh2o=None):
    '''given array of n values in time and list of components, calculate radius'''
    vtot = calcv(components, ns, has_water, xh2o)
    r = (3*vtot/(4*pi))**(1/3)
    return r

def calcp0(a, b, temp):
    '''given regression line parameters, calculate vapor pressure at given temperature.'''
    log_p0 = a + b*(1000./temp)
    p0 = pow(10, log_p0)
    return p0

def plotevap(components, ns, tsteps, dt, labels):
    '''convenience function for quick plot of evaporation'''
    fig, ax = plt.subplots()
    x = np.arange(start=0, stop=(tsteps+1)*dt, step=dt)/3600
    for i, n in enumerate(ns.transpose()):
        ax.plot(x, n, label=labels[i])
    ax.set_ylabel("quantity compound / molec")
    ax.set_xlabel("time / h")
    ax.legend(loc='lower right')

    r = calcr(components, ns)
    ax2 = ax.twinx()
    ax2.plot(x, r, 'k--', label='radius (right axis)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax2.set_ylabel("particle radius / m")
    ax2.legend()
    plt.show()
