pyvap: kinetic model of particle evaporation
============================================

pyvap models evaporation of components from a multicomponent spherical
particle. The model follows the treatment of Maxwellian flux given in
*Atmospheric Chemistry & Physics* by Seinfeld and Pandis. The model is
written in Python 3, with Scipy's ode solver.

Usage
-----

Example code to model the evaporation of an equimolar mixture of
polyethylene trimer and hexamer (PEG-3 and PEG-6) in a particle with
starting radius 5 micron, and create a figure of the output:

::

    import pyvap
    import matplotlib.pyplot as plt
    import numpy as np

    # define model parameters
    peg3 = {"name": "PEG-3",
            "Dg": 5.95e-6,
            "M": 0.1502,
            "rho": 1108.0,
            "cinf": 0,
            "p298": 6.68e-2,
            "delh": 78.3e+3}
    peg6 = {"name": "PEG-6",
            "Dg": 4.26e-6,
            "M": 0.2823,
            "rho": 1180.0,
            "cinf": 0,
            "p298": 3.05e-5,
            "delh": 102.1e+3}
    cmpds = [peg3, peg6]
    # equimolar mixture of PEG-3 and PEG-6
    comp = [0.5, 0.5]
    r_init = 5e-6 # starting radius, m
    time = 2*60*60 # integration time, s
    numpts = 2000 # number of points to integrate
    temp = 298 # K

    # run model
    model = pyvap.analyze_evap(cmpds, comp, r_init, time, numpts, temp,
                               makefig=True)

    # display generated figure during, e.g., interactive ipython session
    plt.show()

    # save generated figure
    evap_fig, (evap_ax0, evap_ax1) = model["evap_fig"]
    evap_fig.savefig("evaporation.png")

    # save csv of evaporation model output (no. of molecules of each component,
    # at each timestep)
    np.savetxt("evaporation.csv", model["evap_a"], delimiter=",")

Installation
------------

Install using ``pip``.

Install most recent Github commit (stability not guaranteed):

::

    pip install git+https://github.com/awbirdsall/pyvap

Dependencies
------------

Developed for Python 3.5.

Requires ``numpy``, ``scipy``, and ``matplotlib>=1.5`` (automatically
handled if using ``pip`` to install). I recommend using
```conda`` <http://conda.pydata.org/docs/index.html>`__ to install the
Scipy stack on a Windows machine if ``pip`` is having issues.

Development
-----------

Posting issues or pull requests to the github page is welcome!
