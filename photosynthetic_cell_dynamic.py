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
import matplotlib.pyplot as plt


# PARAMETERS -------------------------------------------------------
#
# Implicitly reuse all global parameters from steady state model


# VARIABLES --------------------------------------------------------


# initialize model
m = GEKKO()
m.time = time


# list of catalytic rates v for all enzymes
v = pd.Series(
    [m.Var(value=1, lb=0, ub=100, name='v_'+i) for i in enz],
    index=enz)

# list of alpha=fraction of ribosomes engaged in synthesis of protein    
a = pd.Series(
    [m.Param(value=a_optim[i].value, name='a_'+i) for i in pro],
    index=pro)
    
# list of concentration of all components (enzymes and metabolites)
c = pd.Series(
    [m.Var(value=0.2, lb=0, ub=100, name='c_'+i) for i in pro+met],
    index=pro+met)


# hv is the time-dependent light intensity
hv = m.Param(value=light, name='hv')

# growth rate, in this model equals rate of the ribosome
mu = m.Var(value=1, name='mu')

# optional parameter: ribosome translational capacity
sigma = m.Param(value=0.4, name='sigma')

# biomass accumulated over time
#bm = m.Var(value=1, name='bm')


# EQUATIONS --------------------------------------------------------
#
# time-dependent differential equation for change in protein_conc
# with the 'error' being a logistic term approximating -1 and 1
# for big differences between a and c, and 0 for small differences
m.Equations([c[i].dt() == c[i]*sigma*(a[i]-c[i])/(a[i]+c[i]) for i in pro])

# metabolite mass balance
m.Equations([sum(stoich.loc[i]*v) - mu*c[i] == 0 for i in met])

# growth rate mu equals rate of the ribosome v_RIB when a_pro = c_pro = 1,
m.Equation(mu == v['RIB'])

# biomass accumulation over time with initial biomass concentration 0.1
#m.Equation(bm.dt() == mu*bm)

# Michaelis-Menthen type enzyme kinetics
m.Equation(v['LHC'] == kcat['LHC']*c['LHC']*hv**hc['LHC']/(Km['LHC']**hc['LHC'] + hv**hc['LHC'] + (hv**(2*hc['LHC']))/Ki))
m.Equation(v['PSET'] == kcat['PSET']*c['PSET']*c['hvi']**hc['PSET']/(c['hvi']**hc['PSET'] + Km['PSET']**hc['PSET']))
m.Equation(v['CBM'] == kcat['CBM']*c['CBM']*c['nadph']*sub**hc['CBM']*c['atp']/(c['nadph']*sub**hc['CBM']*c['atp'] + KmNADPH*c['atp'] + KmATP*c['nadph'] + KmATP*sub**hc['CBM'] + Km['CBM']**hc['CBM']*c['nadph']))
m.Equation(v['LPB'] == kcat['LPB']*c['LPB']*c['pre']**hc['LPB']/(Km['LPB']**hc['LPB'] + c['pre']**hc['LPB']))
m.Equation(v['RIB'] == kcat['RIB']*c['RIB']*c['pre']**hc['RIB']/(Km['RIB']**hc['RIB'] + c['pre']**hc['RIB']))


# SOLVING ----------------------------------------------------------
#
# Solve using mode 4, which simulates a process using
# the differential equations and parameters
# no degrees of freedom allowed (no optimized variables)
m.options.IMODE=4
m.solve()


# COLLECTING RESULTS -----------------------------------------------
#
# collect results in pandas data frame and save
with open(m.path+'//results.json') as f:
    result = pd.DataFrame(json.load(f))

result.to_csv('/home/michael/Documents/SciLifeLab/Resources/Models/GEKKO/cyano/result_dynamic_DN.csv')
