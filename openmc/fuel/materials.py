"""
openmc.fuel.materials
==============

This module contains common functions to determine materials properties
"""
# 
import numpy as np 
import math 
#
akfuel = np.array([1.05, 0., 0., 0., 2150., 73.15])
akclad = np.array([7.51, 2.09e-2, -1.45e-5, 7.67e-9])
#
#functions for water properties at 15.5 MPa
#cubic polynomial for thermal conductivity, max err=0.0204%
# 15.5 Mpa,  280 < T < 340, k in w/m-C, T in C
def fcond(t):        
    k = 8.9182016e-01 + t*(-2.1996892e-03 + t*(9.9347652e-06 + t*(-2.0862471e-08)))
    return k

#cubic polynomial for density as a function of temperature, max err=0.0692%
#15.5 Mpa,  280 < T < 340, rho in Kg/M^3, T in C
def fdens(t):
    den = 5.9901166e+03 + t*(-5.1618182e+01 + t*(1.7541848e-01 + t*(-2.0613054e-04)))
    return den

# cubic polynomial for density as a function of enthalpy, max err=0.0112%
#  15.5 Mpa,  280 < T(h) < 340, rho in Kg/M^3, input h in J/Kg
def fdensh(h):
    hk=0.001*h
    den = 1.4532039e+03 + hk*(-1.1243975 + hk*(7.4502004e-04 + hk*(-2.3216531e-07)))
    return den

# cubic polynomial for enthalpy as a function of temperature, max err=0.0306%
#  15.5 Mpa,  280 < T < 340, output h in J/Kg, T in C
def fenthal(t):
    y = -5.9301427e+03 + t*(6.5488800e+01 + t*(-2.1237562e-01 + t*(2.4941725e-04)))
    h = y*1000.0
    return h


# quartic polynomial for heat cappacity, max error=0.3053%
# 15.5 Mpa,  280 < T < 340, output h in J/Kg-C, T in C
def fhcap(t):
    y = 3.0455749e+03 + t*(-4.0684599e+01 + t*(2.0411250e-01 + t*(-4.5526705e-04 + t*(3.8115453e-07))))
    cp = y*1000.0
    return cp

# wall-to-coolant heat transfer coeffcient in w/m^2-C
def fhtcoef(t,deq,rhou): 
    k = fcond(t)
    mu = fvisco(t)
    cp = fhcap(t)
    pr = cp*mu/k
    re = deq*rhou/mu
    htc = 0.023*k/deq*(pr**0.4)*(re**0.8)
    return htc

def fhtcoef_if97(t,deq,rhou, P): 
    # use IF97PRO
    from iapws import IAPWS97
    # k  = TCPT(P*1.0e+06, t)
    # mu = DVPT(P*1.0e+06, t)
    # cp = HCPT(P*1.0e+06, t)
    k  = IAPWS97(P=P*1.0e+06, T=t)
    mu = IAPWS97(P=P*1.0e+06, T=t)
    cp = IAPWS97(P=P*1.0e+06, T=t)
    #
    pr = cp*mu/k
    re = deq*rhou/mu
    htc = 0.023*k/deq*(pr**0.4)*(re**0.8)
    return htc

def ftemp(h):
    x = h*1.0e-6
    t = -4.993101541248199e+03+x*( 2.570549364678637e+04+x*(-5.174137664233302e+04+ \
        x*( 5.425153902632183e+04+x*(-3.119846426354911e+04+x*( 9.371717734576585e+03+ \
        x*(-1.153820485592518e+03))))))
    return t

# cubic polynomial for viscosity, max err=0.0641%
#  15.5 Mpa,  280 < T < 340, mu in Pa-sec, T in C
def fvisco(t):
    v = 9.0836878e-04 + t*(-7.4542195e-06 + t*(2.3658072e-08 + t*(-2.6398601e-11)))
    return v

def fkf(t):
    k = akfuel[0] + t*(akfuel[1] + t*(akfuel[2] + t*akfuel[3])) + akfuel[4]/(t-akfuel[5])
    return k

# thermal conductivity of Zr in w/m-C, t in K
def fkc(t):
    k = akclad[0] + t*(akclad[1] + t*(akclad[2] + t*akclad[3]))
    return k

def ftfavg(x, nr, dr2):
    #
    xavg = 0.
    i = 1
    a = (x[i + 1-1] - x[i-1])/dr2
    xavg = (7.*x[i-1] + x[i + 1-1])*0.25
    for l in range(2, nr+1):
        area = 2./3.*l*x[l + 1-1]+44./3.*i*x[l-1]+2./3.*(i - 1)*x[l - 1-1]
        xavg += area
        i = l
    #
    nrm1 = nr - 1
    area = (16./3.*nrm1 + 4. + 5./24.)*x[nr + 1-1]+(10./3.*nrm1 + 2.25)*x[nr-1] - (2./3.*nrm1 + 11./24.)*x[nr - 1-1]
    t = (xavg + area)*dr2/(8.*(nr*nr*dr2))
    return t


def fthcon(ftemp, burnup, fden=0.95, gadoln=0.0)
    """ 
    fthcon calculates the fuel thermal conductivity and its derivative with respect to temperature as a function of
    temperature, density, composition and burnup for UO2 fuel
    
    Input
    
    burnup - current local burnup (MWd/MTU)
    ftemp  - current fuel ring temperature (K)
    comp   - PuO2 content of fuel (percent puo2 in total fuel weight)
    fden   - Current fuel density (ratio of actual density to theoretical density)
    fotmtl - oxygen to metal ratio of fuel (atoms oxygen/atoms metal)
    gadoln - weight fraction of gadolinia in the fuel

    Output

    con    - output fuel thermal conductivity (W/(m*K))
    
    Reference:

    (1) Proposed by staff at NFI, Japan, at the May 1997 ANS Topical Meeting on LWR Fuel Performance
        in Portland, OR: (Ohira, K., and N.Itagaki, 1997. "Thermal Conductivity Measurements of High
        Burnup UO2 Pellet and a Benchmark Calculation of Fuel Center Temperature", proceedings pp. 541-549.
    """
    #bug, h, rphonon, elect, base, ax, cx, fm, x, tc, tco, ucon, po, buguo2, frac
    #
    #Burnup in GWd/MTU
    bug = burnup / 1000.0
    #
    # NFI formula (Ohira & Itagaki, ANS LWR Fuel perf. Topical mtg. 1997). 
    # MODIFIED in January 2002 to raise low-burnup thermal conductivity at low temperature
    # and to lower thermal conductivity at very high temp.
    h = 1.0 / (1.0 + 396.0 * math.exp(-6380.0 / ftemp))
    rphonon = 1.0 / (0.0452 + 0.000246 * ftemp + 1.0 * 0.00187 * bug + 1.1599 * gadoln + \
            (1.0 - 0.9 * math.exp(-0.04 * bug)) * 0.038 * (bug ** 0.28) * h)
    elect = (3.50e9 / ftemp ** 2) * math.exp(-16361.0 / ftemp)
    base = rphonon + elect
    # fm is the Lucuta porosity correction factor(applied to 100% TD fuel)
    fm  = fden / (1.0 + 0.5 * (1.0 - fden))
    # NFI base equation is for 95% TD fuel, so multiply by 1.079 to raise to 100% TD fuel conductivity, then multiply by fm
    con = base * fm * 1.079
    # Find uncertainty
    #sigftc = 1.0 # or 0.0  
    #con = con * (1.0 + abs(sigftc * 0.088))
    #
    #fmelt = 3113.15 - 5.0 * burnup / 10000.0 - 4.8 * gadoln * 100 # for UO2 Fuel
    #
    #if (ftemp < ftmelt):
    #    ucon = 0.2 * (1.0 + abs(2.0 - fotmtl) * 10.0)
    #else:
    #    ucon = con / 2.0
    return con

def fthcon_ctf(ftemp):
    """
    fuel thermal conductivity extracted from CTF code
    
    Input
    
       ftemp: fuel temperature in K
    """
    # Formulation
    tt = ftemp - 273.15
    k_ctf = max(1.1038,2.335e+03/(tt + 464.)) + 7.027e-03 * math.exp(1.867e-03 * tt)
    # Unit conversion
    t_btu_kJ = 1.05505585262;
    t_ft_in=12.0;
    t_in_ft=1./t_ft_in;
    t_in_m = 0.0254;
    t_ft_m = t_ft_in * t_in_m;
    t_K_F = 1.8;
    t_F_K = 1. / t_K_F;
    tc_US_SI = t_btu_kJ * 1E3/(3.6E3 * t_ft_m * t_F_K);
    #
    #Thermal Conductivity (SI --> US)
    #[W/(m K)] <--> [btu/(h ft F)]
    #
    k_ctf_SI = k_ctf * tc_US_SI
    return k_ctf_SI

    
