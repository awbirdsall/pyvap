# pyvap : simple model of particle evaporation
# author: Adam Birdsall

from scipy.constants import pi, k
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from math import pow

def dn(t, y, Ma, rho, cinf, po, Dg, T):
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
    
    Outputs
    -------
    result : ndarray
    1D-array describing dn for all components.
    '''
    v = y*Ma/rho # array of partial volumes, m^3
    vtot = v.sum() # particle volume, m^3
    r = (3*vtot/(4*pi))**(1/3) # radius, m^3
    # print("calced r: {:.3e}".format(r))
    x = y/y.sum() # mole fractions
    # assume ideality in vapor pressure calculation
    cs = x*po/(k*T) # gas-phase concentration at surface, molec m^-3
    result = 4*pi*r*Dg*(cinf-cs)
    return result

def evaporate(components, ninit, T, tsteps, dt):
    '''calculate evaporation of multicomponent particle.'''
    
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
    output = np.empty((tsteps+1, len(components)))
    output[0, :] = ninit

    r = ode(dn)
    r.set_integrator('lsoda', with_jacobian=False,)
    r.set_initial_value(ninit, t=0)
    r.set_f_params(Ma, rho, cinf, po, Dg, T)

    # integrate and save output
    while r.successful() and r.t < tsteps*(dt-1):
        entry = int(round(r.t/dt))+1
        nextstep = r.integrate(r.t + dt)
        # print("{}: dn at {:.1e} s: {}".format(entry, r.t, dn(r.t, r.y, Ma, rho, cinf, po, Dg, T)))
        output[entry, :] = nextstep

    return output

def calcv(components, ns):
    '''calculate volume of particle over time for given components'''
    vtot = np.empty(ns.shape[0])
    for i, c in enumerate(components):
        v = ns[:,i]*c['Ma']/c['rho']
        vtot = vtot + v
    return vtot

def calcr(components, ns):
    '''given array of n values in time and list of components, calculate radius'''
    vtot = calcv(components, ns)
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
