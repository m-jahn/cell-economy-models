#
# PHOTOAUTOTROPHIC CELL MODEL
# version: 2.0
# subversion: 'light and CO2 limitation dynamic'
# date: 2018-08-27
# author: Michael Jahn
# affiliation: Science for Life Laboratory (KTH), Stockholm, Sweden
# based on: R. Burnap, 2015, Molenaar et al., 2009
# ported from GAMS to python GEKKO
# characteristics: protein economy model of a photoautotrophic cell
#
#
# LIBRARIES ------------------------------------------------------------

from gekko import GEKKO
import json
import re
import math
import numpy as np
import pandas as pd


# PARAMETERS -----------------------------------------------------------
#
# Implicitly reuse all global parameters from steady state model
# Explicitly use optimal alpha mass fractions and starting concentrations


# INITIALIZE DYNAMIC MODEL ---------------------------------------------
# alternative: server='http://xps.apmonitor.com'
def dynamic(rs, a ,c):
    
    m = GEKKO(remote = True, server = 'http://xps.apmonitor.com')
    m.options.IMODE = 5
    m.options.REDUCE = 1
    m.options.MAX_ITER = 100
    m.time = time
    
    
    # VARIABLES --------------------------------------------------------
    #
    # list of catalytic rates v for all enzymes
    v = pd.Series(
        [m.Var(value = 1, lb = 0, ub = 100, name = 'v_' + i) for i in enz],
        index = enz)
    
    # list of alpha = fraction of ribosomes engaged in synthesis of protein
    a = pd.Series(
        [m.Param(value = a[i].value, name = 'a_' + i) for i in pro],
        index = pro)
        
    # list of concentration of all components (enzymes and metabolites)
    c = pd.Series(
        [m.Var(value = c[i][1], lb = 0, ub = 500, name = 'c_' + i) for i in pro + met],
        index = pro + met)
    
    # optional protein utilization, µ dependent: u = c - reserve rs
    u = pd.Series(
        [m.Var(value = 1, lb = 0, ub = 1, name = 'u_' + i) for i in enz],
        index = enz)
        
    # hv is the time-dependent light intensity
    hv = m.Param(value = light, name = 'hv')
    
    # growth rate, here simply equals rate of the ribosome
    mu = m.Var(value = 1, name = 'mu')
    
    # biomass accumulated over time with initial value
    bm = m.Var(value = 1, name = 'bm')
    
    
    # EQUATIONS --------------------------------------------------------
    #
    # time-dependent differential equation for change in protein_conc
    # with the 'error' being a logistic term approximating -1 and 1
    # for big differences between a and c, and 0 for small differences
    m.Equations([c[i].dt() == v['RIB']*(a[i]-c[i])/(a[i]+c[i]) for i in pro])
    
    # metabolite mass balance
    m.Equations([sum(stoich.loc[i]*v) - mu*c[i] == 0 for i in met])
    
    # growth rate mu equals rate of the ribosome v_RIB when a_pro = c_pro = 1,
    m.Equation(mu == v['RIB'])
    
    # utilized enzyme fraction = total enzyme - reserve
    m.Equations([u[i] == c[i]-rs[i]*(1-mu/mumax) for i in enz])
    
    # biomass accumulation over time
    m.Equation(bm.dt() == mu*bm)
    
    
    # Michaelis-Menthen type enzyme kinetics
    m.Equation(v['LHC'] == kcat['LHC']*u['LHC']*hv**hc['LHC']/(Km['LHC']**hc['LHC'] + hv**hc['LHC'] + (hv**(2*hc['LHC']))/Ki))
    m.Equation(v['PSET'] == kcat['PSET']*u['PSET']*c['hvi']**hc['PSET']/(c['hvi']**hc['PSET'] + Km['PSET']**hc['PSET']))
    m.Equation(v['CBM'] == kcat['CBM']*u['CBM']*c['nadph']*sub**hc['CBM']*c['atp']/(c['nadph']*sub**hc['CBM']*c['atp'] + KmNADPH*c['atp'] + KmATP*c['nadph'] + KmATP*sub**hc['CBM'] + Km['CBM']**hc['CBM']*c['nadph']))
    m.Equation(v['LPB'] == kcat['LPB']*u['LPB']*c['pre']**hc['LPB']/(Km['LPB']**hc['LPB'] + c['pre']**hc['LPB']))
    m.Equation(v['RIB'] == kcat['RIB']*u['RIB']*c['pre']**hc['RIB']/(Km['RIB']**hc['RIB'] + c['pre']**hc['RIB']))
    
    
    # SOLVING ----------------------------------------------------------
    #
    # solving maximizing specific growth rate; 
    # objective is always minimized, so that we have 
    # to state -1*obj to maximize it
    m.Obj(-mu)
    m.solve()
    
    # collect results and
    return(result("dynamic", m, v, a, c, u))
