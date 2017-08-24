# pyvap : simple model of particle evaporation
# author: Adam Birdsall

from __future__ import division
from math import e
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, k, R, N_A
from scipy.integrate import ode

MA_H2O = 0.018/6.02e23 # mass water, kg molec^-1
RHO_H2O = 1.000e3 # density water, kg m^-3

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
        vh2o = nh2o * MA_H2O / RHO_H2O
        vtot = v.sum()+vh2o # particle volume, m^3
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
    '''Calculate evaporation of multicomponent particle.

    Parameters
    ----------
    components : list
    List of components of particles, each represented as a dict of parameters.

    ninit : numpy.ndarray
    1D array of starting number of molecules of each component in `components`.

    T : float
    Temperature of evaporation, K.

    num : int
    Total number of integrated time points, including t=0.

    dt : float
    Integration timestep, s.

    has_water : boolean (optional, default False)
    Toggle whether the presence of water is considered in calculating
    evaporation (separate from list of `components`).

    xh2o : float (optional, default None)
    Fixed mole fraction of water to include in particle. Only considered if
    `has_water` is True.

    Returns
    -------
    output : numpy.ndarray
    2D array of results from integration: number of molecules of each component
    remaining at each timestep. First index is along `num` timesteps and second
    index is along `len(components)` components.
    '''

    # extract data from components and build data arrays
    Ma = np.empty(len(components))
    rho = np.empty(len(components))
    cinf = np.empty(len(components))
    po = np.empty(len(components))
    Dg = np.empty(len(components))
    for i, component in enumerate(components):
        # use Ma (molecular weight, kg molec^-1) if provided
        # otherwise use M (molar weight, kg  mol^-1).
        if 'Ma' in component:
            Ma[i] = component['Ma']
        else:
            Ma[i] = component['M']/N_A
        rho[i] = component['rho']
        cinf[i] = component['cinf']
        # use p298 and delh if available. otherwise use p0_a and p0_b
        if ('p298' in component) and ('delh' in component):
            p0a, p0b = convert_p0_enth_a_b(component['p298'],
                                           component['delh'], 298.15)
            po[i] = calcp0(p0a, p0b, T)
        else:
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
    '''Calculate volume of particle over time for given components.
    
    Parameters
    ----------
    components : list
    List of dicts for each component.
    
    ns : numpy.array
    2D numpy array of molar amounts of material. First index is entry number
    (e.g., each timestep), second index is index of component in `components`.
    
    has_water : Boolean (optional, default False)
    Whether implicit water is added in addition to `components`.
    
    xh2o : float (optional, default None)
    Fixed mole fraction of water added to particle in calculating value. Only
    considered if has_water is True.
    
    Returns
    -------
    vtot : numpy.array
    1D array of total volumes calculated for each row in `ns`, with possible
    addition of water, in m^3.
    '''
    vtot = np.zeros_like(ns.shape[0])
    for i, c in enumerate(components):
        # convert number of molecules in ns, using molecular/molar mass Ma/M
        # (units kg (molec or mol)^-1) and density rho (kg m^-3), to volume in
        # m^3
        if 'Ma' in c:
            v = ns[:, i] * c['Ma'] / c['rho']
        else:
            v = ns[:, i] * c['M'] / c['rho'] / N_A
        vtot = vtot + v
    if has_water:
        # water not explicitly tracked, instead fixed
        # xh2o for fixed ambient RH (ah2o)
        nh2o = xh2o/(1-xh2o)*ns.sum(axis=1)
        vh2o = nh2o * MA_H2O / RHO_H2O
        vtot = vtot + vh2o
    return vtot

def calcr(components, ns, has_water=False, xh2o=None):
    '''Given array of n values in time and list of components, calculate radii.

    Parameters
    ----------
    components : list
    List of dicts for each component.
    
    ns : numpy.array
    2D numpy array of molar amounts of material. First index is entry number
    (e.g., each timestep), second index is index of component in `components`.
    
    has_water : Boolean (optional, default False)
    Whether implicit water is added in addition to `components`.
    
    xh2o : float (optional, default None)
    Fixed mole fraction of water added to particle in calculating value. Only
    considered if has_water is True.

    Returns
    -------
    r : numpy.array
    Array of radii, in m, for each row of components given in `ns`.

    '''
    vtot = calcv(components, ns, has_water, xh2o)
    r = (3*vtot/(4*pi))**(1/3)
    return r

def convert_p0_enth_a_b(p0, del_enth, t0):
    '''Convert p0 and delta enthalpy to vapor pressure temp dependence params.

    Parameters
    ----------
    p0 : float or ndarray
    Vapor pressure at reference temperature, Pa.
    del_enth : float or ndarray
    Enthalpy of vaporization (or sublimation), J mol^-1.
    t0 : float or ndarray
    Reference temperature for p0 value, K.

    Returns
    -------
    p0_a, p0_b : float
    a (intercept, Pa) and b (slope, 1000/K) linear regression parameters for
    temperature dependence of log10(vapor pressure).

    '''
    p0_a = 1/np.log(10) * ((del_enth/(R*t0)) + np.log(p0))
    p0_b = -del_enth/(1000*np.log(10)*R)
    return p0_a, p0_b

def calcp0(a, b, temp):
    '''given regression line parameters, calculate vapor pressure at given
    temperature.'''
    log_p0 = a + b*(1000./temp)
    p0 = pow(10, log_p0)
    return p0

def plot_evap(x, molec_data, r_data, series_labels, xlabel):
    fig, ax = plt.subplots()
    if series_labels is not None:
        for i in np.arange(molec_data.shape[1]):
            ax.plot(x, molec_data[:, i], label=series_labels[i])
        ax.legend(loc='lower right', ncol=3)
    else:
        ax.plot(x, molec_data[:, i])

    ax.set_ylabel("quantity compound / molec")
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, np.max(molec_data[0, :])*1.1)

    ax2 = ax.twinx()
    ax2.plot(x[:-1], r_data[:-1], 'k--', label='radius (right axis)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax2.set_ylabel("particle radius / m")
    ax2.set_ylim(0, r_data[0]*1.1)
    ax2.legend()

    return fig, (ax, ax2)

def efold_time(t, y):
    '''Calculate "e-fold time" given timeseries and corresponding y values.

    Return -1 if e-fold time not reached in given series.

    Assume that "e-fold time" is sensical values for data series (i.e.,
    monotonic decay) and that value at t=inf is 0.'''
    efold_a = np.where(y <= 1./e*y[0])[0]
    if efold_a.size > 0:
        efold_time = t[efold_a[0]]
    else:
        efold_time = -1
    return efold_time


def analyze_evap(cmpds, comp, r, t, num, temp, makefig=False, complabels=None,
                 has_water=False, xh2o=None):
    '''All-in-one function to run kinetics model and plot and return output.
    
    Parameters
    ----------
    cmpds : list
    List of compounds in particle, each represented as a dict of parameters.
    Required keys in each dict include: `name`, `Dg` (gas-phase diffusivity,
    m^2 s^-1), `Ma` or `M` (molecular mass or molar mass, kg molec^-1 or kg
    mol^-1), `rho` (density, kg m^-3), `cinf` (gas-phase concentration at
    infinite distance from particle surface), `p298` and `delh` or `p0_a` and
    `p0_b` (vapor pressure at 298.15 K, Pa, and delta H of vaporization,
    J mol^-1; or linear fit parameters to calculate intercept, Pa, and slope,
    1000/T, for temperature dependence of vapor pressure).

    comp : list or numpy.ndarray
    List or 1D array of initial composition of particle, given as relative molar
    abundance of each compound in `cmpds`. Does not need to be normalized to
    1 or total number of molecules.

    r : float
    Starting particle radius, m.

    t : float
    Total integration time, s.

    num : int
    Total number of integrated time points, including t=0.

    temp : float
    Temperature of evaporation, K.

    makefig : boolean (optional, default False)
    Whether to make default output plot, using `plot_evap()`. If True, output
    dict is updated with `plot_evap()` output (`fig`, (`ax`, `ax2`)) stored
    under the `evap_fig` key.

    complabels : list (optional, default None)
    List of strings to use in legend of a produced figure, corresponding to
    each compound in `cmpds`. If None, use `name` in each compound dict.

    has_water : boolean (optional, default False)
    Toggle whether the presence of water is considered in calculating
    evaporation (separate from list of `components`).

    xh2o : float (optional, default None)
    Fixed mole fraction of water to include in particle. Only considered if
    `has_water` is True.

    Returns
    -------
    output_dict : dict
    Dictionary of all output. Keys include `evap_a` (array of molecular
    quantities at each timestep from `evaporate()`, `t_a` (array of timesteps
    corresponding to `evap_a` output), `r_a` (array of particle radius at each
    timestep, calculated from `calcr()`), `efold_dict` (dict of characteristic
    e-folding evaporation time for each component, calculated from
    `efold_time()`), and `evap_fig` (only present if `makefig` is True).
    
    '''
    output_dict = dict()

    # calc initial total num molecules
    avg_rho = np.average([x['rho'] for x in cmpds],
                         weights=comp) # kg m^-3
    total_mass = 4./3.*pi * r**3 * avg_rho # kg
    # use either Ma or M
    get_Ma = lambda cmpd: cmpd['Ma'] if 'Ma' in cmpd else cmpd['M']/N_A
    avg_molec_mass = np.average([get_Ma(x) for x in cmpds],
                                weights=comp) # kg molec^-1
    total_molec = total_mass/avg_molec_mass # molec

    if type(comp)==list:
        comp = np.array(comp)
    ncomp = comp/comp.sum() # make sure composition is normalized to 1.
    molec_init = ncomp*total_molec

    # set up and integrate ODE
    t_a, dt_evap = np.linspace(start=0, stop=t, num=num, retstep=True)
    evap_a = evaporate(cmpds, ninit=molec_init, T=temp, num=num, dt=dt_evap,
                       has_water=has_water, xh2o=xh2o)
    output_dict.update({'t_a': t_a, 'evap_a': evap_a})

    # back out radius timeseries

    r_a = calcr(cmpds, evap_a, has_water, xh2o)
    output_dict.update({'r_a': r_a})

    # plot
    if complabels is None:
        complabels = [x['name'] for x in cmpds]

    if makefig:
        xlabel = "time / h"
        fig, (ax, ax2) = plot_evap(x=output_dict['t_a']/3600,
                                   molec_data=evap_a,
                                   r_data=r_a,
                                   series_labels=complabels,
                                   xlabel=xlabel)
        output_dict.update({'evap_fig': (fig, (ax, ax2))})
    else:
        output_dict.update({'evap_fig': None})

    # e-folding times, converted from seconds to hours
    efold_dict = {l: efold_time(output_dict['t_a']/3600, evap_a[:, i])
                  for i, l in enumerate(complabels)}
    output_dict.update({'efold_dict': efold_dict})

    return output_dict
