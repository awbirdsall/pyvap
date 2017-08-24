from __future__ import division
import pytest
import pyvap as pv
import numpy as np
from math import pi

# global constants
N_A = 6.022140857e+23
R = 8.3144598

@pytest.fixture(scope='module')
def unitcmpd():
    # some nice round numbers for convenience
    unitcmpd = {"name": "unitcmpd",
                "Dg": 1e-6,
                "M": 0.100,
                "rho": 1000.0,
                "cinf": 0,
                "p298": 1.00e-2,
                "delh": 100.0e+3}
    return unitcmpd

@pytest.fixture(scope='module')
def peg3():
    peg3 = {"name": "PEG-3",
            "Dg": 5.95e-6,
            "M": 0.1502,
            "rho": 1108.0,
            "cinf": 0,
            "p298": 6.68e-2,
            "delh": 78.3e+3}
    return peg3

@pytest.fixture(scope='module')
def peg6():
    peg6 = {"name": "PEG-6",
            "Dg": 4.26e-6,
            "M": 0.2823,
            "rho": 1180.0,
            "cinf": 0,
            "p298": 3.05e-5,
            "delh": 102.1e+3}
    return peg6

@pytest.fixture(scope='module')
def evaporate_binary_output(peg3, peg6):
    evaporate_binary_output = pv.evaporate([peg3, peg6],
                                           ninit=[N_A/100, N_A/100],
                                           T=298,
                                           num=100,
                                           dt=10,
                                           has_water=False, xh2o=None)
    return evaporate_binary_output

@pytest.fixture(scope='module')
def analyze_evap_binary_output(peg3, peg6):
    analyze_evap_binary_output = pv.analyze_evap([peg3, peg6],
                                                 [0.5, 0.5],
                                                 5e-6,
                                                 3600,
                                                 100,
                                                 298,
                                                 makefig=False,
                                                 complabels=None,
                                                 has_water=False, xh2o=None)
    return analyze_evap_binary_output

def test_calcv_calcs_volume(unitcmpd):
    # consider 1 g of material, which is 1/100 of a mole for mw 100 g/mol.
    n = N_A/100
    v = pv.calcv([unitcmpd], np.array([[n]]), has_water=False, xh2o=None)
    # volume should be about 1 mL, which is 1 cm^3 or 1e-6 m^3
    assert np.allclose(v, 1e-6)

def test_calcr_calcs_radius(unitcmpd):
    # consider 4/3*pi*(1e-2 m)^3 of material, or 4/3*pi grams, which should
    # result in a radius of 1e-2 m.
    grams = 4/3*pi
    moles = N_A/100 * grams
    r = pv.calcr([unitcmpd], np.array([[moles]]), has_water=False, xh2o=None)
    assert np.allclose(r, 1e-2)

def test_roundtrip_convert_p0_enth_single_temp(unitcmpd):
    # Check conversion from p298 and del_enth to a and b linear parameters
    # correctly gives p298 back when calculating vapor pressure at 298.15 K.
    # Does not check whether conversion gets temperature dependence right.
    p298 = unitcmpd['p298']
    del_enth = unitcmpd['delh']
    t0 = 298.15
    p0_a, p0_b = pv.convert_p0_enth_a_b(p298, del_enth, t0)
    assert np.allclose(pv.calcp0(p0_a, p0_b, t0), p298)

def test_efoldtime_calcs_efold_time(unitcmpd):
    # test efolding time calculation for dummy exponential decay with
    # characteristic time tao in fact returns tao
    t = np.arange(0,100)
    tao = 50
    assert pv.efold_time(t, np.exp(-t/tao)) == tao

def test_dn_calcs_binary_particle(peg3, peg6):
    # directly calculate dn for peg3 and peg6 particle. Regression test result
    # against result of 24/08/2017.
    t = 0
    y = np.array([N_A/100, N_A/100])
    Ma = np.array([peg3['M']/N_A, peg6['M']/N_A])
    rho = np.array([peg3['rho'], peg6['rho']])
    cinf = np.array([peg3['cinf'], peg6['cinf']])
    po = np.array([peg3['p298'], peg6['p298']])
    Dg = np.array([peg3['Dg'], peg6['Dg']])
    T = 298.15
    dn = pv.dn(t, y, Ma, rho, cinf, po, Dg, T, has_water=False, xh2o=None)
    assert np.allclose(dn, np.array([ -5.84598631e+12,  -1.91105772e+09]))

def test_evaporate_output_correct_numpts(evaporate_binary_output):
    assert evaporate_binary_output.shape[0] == 100

def test_analyze_evap_correct_keys(analyze_evap_binary_output):
    # test stability of interface of analyze_evap output: produces exactly
    # following keys (when makefig is False)
    model_keys = ['evap_fig', 'efold_dict', 'r_a', 't_a', 'evap_a']
    for k in analyze_evap_binary_output:
        assert k in model_keys
    for mk in model_keys:
        assert mk in analyze_evap_binary_output.keys()

def test_analyze_evap_binary_output_final_value(analyze_evap_binary_output):
    # regression test against final result of analyze_evap_binary_output
    # evaporation array from 24/08/2017
    assert np.allclose(analyze_evap_binary_output['evap_a'][-1],
                       np.array([4.13439828e2, 8.28251764e11]))
