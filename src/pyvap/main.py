# pyvap : simple model of particle evaporation
# author: Adam Birdsall

from scipy.constants import pi, k
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from math import pow, e

Ma_h2o = 0.018/6.02e23
rho_h2o = 1.000e3

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

def evaporate(components, ninit, T, num, dt, has_water=False, xh2o=None):
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

def plot_evap(x, molec_data, r_data, series_labels, xlabel):
    fig, ax = plt.subplots()
    if series_labels is not None:
        for i in np.arange(molec_data.shape[1]):
            ax.plot(x, molec_data[:,i], label=series_labels[i])
        ax.legend(loc='lower right', ncol=3)
    else:
        ax.plot(x, molec_data[:,i])

    ax.set_ylabel("quantity compound / molec")
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, np.max(molec_data[0,:])*1.1)

    ax2 = ax.twinx()
    ax2.plot(x[:-1], r_data[:-1], 'k--', label='radius (right axis)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax2.set_ylabel("particle radius / m")
    ax2.set_ylim(0, r_data[0]*1.1)
    ax2.legend()
    
    return fig, (ax, ax2)

def efold_time(t, y):
    '''Calculate "e-fold time" given timeseries and corresponding y values.
    
    Return -1 if e-fold time not reached in given series.
    
    Assume that "e-fold time" is sensical values for data series (i.e.,
    monotonic decay) and that value at t=inf is 0.'''
    efold_a = np.where(y<=1./e*y[0])[0]
    if efold_a.size > 0:
        efold_time = t[efold_a[0]]
    else:
        efold_time = -1
    return efold_time


def analyze_evap(cmpds, comp, complabels, r, t, num, temp, makefig=False, has_water=False, xh2o=None):
    '''all-in-one function to run kinetics model and plot and return output'''
    output_dict = dict()
    
    # calc initial total num molecules
    avg_rho = np.average([x['rho'] for x in cmpds],
                         weights=comp) # kg m^-3
    total_mass = 4./3.*pi * r**3 * avg_rho # kg
    avg_molec_mass = np.average([x['Ma'] for x in cmpds],
                                weights=comp) # kg molec^-1
    total_molec = total_mass/avg_molec_mass # molec
    
    ncomp = comp/comp.sum() # make sure composition is normalized to 1.
    molec_init = ncomp*total_molec

    # set up and integrate ODE
    t_a, dt_evap = np.linspace(start=0, stop=t, num=num, retstep=True)
    evap_a = evaporate(cmpds, ninit=molec_init,
                       T=temp, num=num, dt=dt_evap, has_water=has_water, xh2o=xh2o)
    output_dict.update({'t_a': t_a, 'evap_a': evap_a})
    
    # back out radius timeseries
    
    r_a = calcr(cmpds, evap_a, has_water, xh2o)
    output_dict.update({'r_a': r_a})
    
    # plot
    if makefig:
        xlabel = "time / h"
        # diag_extra = "tstep={:.2e} s, temp={} K".format(tstep, temp)
        fig, (ax, ax2) = plot_evap(x=output_dict['t_a']/3600,
                            molec_data=evap_a,
                            r_data=r_a,
                            series_labels=complabels,
                            xlabel=xlabel)
        output_dict.update({'evap_fig': (fig, (ax, ax2))})
    else:
        output_dict.update({'evap_fig': None})
    
    # e-folding times, converted from seconds to hours
    efold_dict = {l: efold_time(output_dict['t_a']/3600, evap_a[:, i]) for i, l in enumerate(complabels)}
    output_dict.update({'efold_dict': efold_dict})
        
    return output_dict