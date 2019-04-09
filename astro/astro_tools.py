# -*- coding: utf-8 -*-
"""
Astro functions (
"""


import os 
import astropy.io.fits as fits
import operator

from scipy import (special, log10, array, sqrt, sin, 
                   exp, log, average, 
                   arange, meshgrid, std)
from numpy.random import random_sample     
import numpy as np  

def CC(z, H0=67.3, WM=0.315, WV=0.685, v = 0):
    """ Cosmo calculator given the cosmological parameters z, H0, WM, WV returns
    the luminosity distance, angular seperation, and possibly many more. Based on
    http://www.astro.ucla.edu/~wright/CC.python
    """

    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr

    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    n = 5000
    i = arange(n)
    
    if not hasattr(z, '__iter__'):
        z = np.array([float(z)])

    zage_Gyra = np.array([])
    for zs in z:
        az = 1.0 / (1 + zs)
        a = az * (i + 0.5) / n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = sum(1./adot)
        zage = az*age/n
        zage_Gyr = (Tyr/H0)*zage
        zage_Gyra = np.append(zage_Gyra, zage_Gyr)

    if v == 'age':
        return zage_Gyra
        
    DTT, DCMR = 0.0, 0.0
    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    a = az+(1-az)*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    DTT = sum(1./adot)
    DCMR = sum(1./(a*adot))

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR
    # tangential comoving distance
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0: ratio =  0.5*(exp(x)-exp(-x))/x
        else: ratio = sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL
    # comoving volume computation
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0: ratio = (0.125*(np.exp(2.*x)-np.exp(-2.*x))-x/2.)/(x*x*x/3.)
        else: ratio = (x/2. - np.sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/5. + (2./105.)*y*y
    VCM = ratio*DCMR*DCMR*DCMR/3.
    V_Gpc = 4.*np.pi*((0.001*c/H0)**3)*VCM
    DL_cm = DL_Mpc * 3.08568E24

    if v == 1:
        print ('\tH_0 = %1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = %1.2f' % WV + ', z = ' + '%1.3f' % z)
        print ('\tIt is now %1.3f' % age_Gyr + ' Gyr since the Big Bang.')
        print ('\tAge at redshift z was %1.3f' % zage_Gyr + ' Gyr.')
        print ('\tLight travel time was %1.3f' % DTT_Gyr + ' Gyr.')
        print ('\tComoving radial distance is \t%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % DCMR_Gyr + ' Gly.')
        print ('\tComoving volume within redshift z ' + '%1.1f' % V_Gpc + ' Gpc^3.')
        print ('\tAngular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or %1.1f' % DA_Gyr + ' Gly.')
        print ('\tAngular scale of %.2f' % kpc_DA + ' kpc/".')
        print ('\tLuminosity distance D_L is %1.1f' % DL_Mpc + ' Mpc or ' + '%1.4e' % DL_cm + ' cm.')
        print ('\tDistance modulus, m-M, is %1.2f mag' % (5*log10(DL_Mpc*1e6)-5))
        print ('\tK-correction for equal effective wavelength %1.2f' %(-2.5*log10(1+z)))
        return(DL_Mpc, kpc_DA)
    elif v == 2:
        return(DL_Mpc, kpc_DA)
    else:
        return(DL_Mpc)
