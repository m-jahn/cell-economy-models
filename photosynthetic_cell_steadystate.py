#
# PHOTOAUTOTROPHIC CELL MODEL
# version: 2.0
# subversion: 'light and CO2 limitation steady state'
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


# PARAMETERS -------------------------------------------------------

# organize variables in sets to simplify indexing
enz = ['LHC', 'PSET', 'CBM', 'LPB', 'RIB']      # enzymes
pro = enz+['MAI']                               # proteins
met = ['hvi', 'atp', 'nadph', 'pre', 'lip']     # metabolites
mem = ['cpm', 'thy']                            # membrane lipids
thyP = ['LHC', 'PSET']                          # thylakoid membrane located proteins
#cmpP = []                                      # cytoplasmic membrane located proteins


# enzyme kinetic parameters as pandas series 
# (kcat, Km, Hill coefficient)
kcat = pd.Series([172, 35, 6, 6, 11], index=enz)
Km = pd.Series([58, 108, 229, 13, 128], index=enz)
hc = pd.Series([2.0043, 1.3989, 2.2477, 0.6744, 0.8659], index=enz)


KmATP = 1   # affinity constant of CBM for ATP
KmNADPH = 1 # affinity constant of CBM for NADPH


# specific surface area of membrane located components
#spA = pd.Series([1, 1, 1, 1], index=mem)


# reaction stoichiometry matrix of met x enz
# as pandas data frame
stoich = pd.DataFrame([
    [1, -1,  0,  0,  0], 
    [0,  1, -1,  0,  0],
    [0,  1, -1,  0,  0],
    [0,  0,  1, -1, -1],
    [0,  0,  0,  1,  0]],
    index=met,
    columns=enz)


# define starting values
sub = 100                                                           # initial substrate concentration, CO2/HCO3-
Ki = 5000                                                           # light inhibition constant for photosystems
light = 100.0/1.5**np.array(range(0,10))                            # light in % max intensity, log decrease
#light = np.array([5.0]*25+[50.0]*26)                                # light in % max intensity, step
#light = np.array([5.0]*24+[50.0]*24+[5.0]*24+[50.0]*24+[5.0]*25)    # light in % max intensity, pulse
time = np.linspace(0, (len(light)-1)*2, len(light))                 # time as a function of light step number. Fewer steps is faster computation


# VARIABLES --------------------------------------------------------


# initialize model
m = GEKKO(remote=True, server='http://xps.apmonitor.com') # alternative: server='http://xps.apmonitor.com'
m.options.IMODE = 5
m.options.REDUCE = 1
m.options.MAX_ITER = 500
m.time = time


# variables calculated by solver to meet the constraints of 
# equations and parameters
# syntax: v = m.Var([starting value], [lb], [ub], [integer], [name]):

# list of catalytic rates v for all enzymes
v = pd.Series(
    [m.Var(value=1, lb=0, ub=100, name='v_'+i) for i in enz],
    index=enz)

# list of alpha = fraction of ribosomes engaged in synthesis of protein
a = pd.Series(
    [m.Var(value=1, lb=0, ub=1, name='a_'+i) for i in pro],
    index=pro)
    
# list of concentration of all components (enzymes and metabolites)
c = pd.Series(
    [m.Var(value=1, lb=0, ub=100, name='c_'+i) for i in pro+met+mem],
    index=pro+met+mem)


# hv is the time-dependent light intensity
hv = m.Param(value=light, name='hv')

# growth rate as variable that becomes the objective function
mu = m.Var(value=1, name='mu')

# biomass accumulated over time with initial value
bm = m.Var(value=1, name='bm')

# volume-to-surface ratio, increases with sphericity of a cell
#beta = m.Var(value=1, lb=0.1, ub=10)


# EQUATIONS --------------------------------------------------------
#
# equations constrain the solution space using parameters;
# they outline the topology of the model

# alpha is fraction of ribosomes engaged in synthesis of protein x
m.Equation(sum(a) == 1)

# protein mass balance: left side, rate of ribosome dedicated to 
# synthesis of each protein, right side, growth rate times protein conc
m.Equations([a[i]*v['RIB'] - mu*c[i] == 0 for i in pro])

# metabolite mass balance: left side, production of metabolites by
# the respective enzyme, right side, growth rate times metabolite conc
m.Equations([sum(stoich.loc[i]*v) - mu*c[i] == 0 for i in met])

# biomass accumulation over time
m.Equation(bm.dt() == mu*bm)

# Michaelis-Menthen type enzyme kinetics
m.Equation(v['LHC'] == kcat['LHC']*c['LHC']*hv**hc['LHC']/(Km['LHC']**hc['LHC'] + hv**hc['LHC'] + (hv**(2*hc['LHC']))/Ki))
m.Equation(v['PSET'] == kcat['PSET']*c['PSET']*c['hvi']**hc['PSET']/(c['hvi']**hc['PSET'] + Km['PSET']**hc['PSET']))
m.Equation(v['CBM'] == kcat['CBM']*c['CBM']*c['nadph']*sub**hc['CBM']*c['atp']/(c['nadph']*sub**hc['CBM']*c['atp'] + KmNADPH*c['atp'] + KmATP*c['nadph'] + KmATP*sub**hc['CBM'] + Km['CBM']**hc['CBM']*c['nadph']))
m.Equation(v['LPB'] == kcat['LPB']*c['LPB']*c['pre']**hc['LPB']/(Km['LPB']**hc['LPB'] + c['pre']**hc['LPB']))
m.Equation(v['RIB'] == kcat['RIB']*c['RIB']*c['pre']**hc['RIB']/(Km['RIB']**hc['RIB'] + c['pre']**hc['RIB']))


# OPTIONAL CONSTRAINTS
#
# total intracellular protein concentration is constrained
m.Equation(sum(c[pro]) == 1)

# membrane composition (c_protein should not exceed c_lipid)
# and minimal amount of membrane is 0.1
m.Equation(             0.1 <= c['cpm'])
m.Equation(sum(c[thyP])+0.1 <= c['thy'])

# lipid balance: lipids are sum of cytoplasmic and thylakoid membrane
m.Equation(sum(c[mem]) == c['lip'])

# fix the mass fraction of maintenance proteins (or others)
m.Equation(a['MAI'] == 0.3)
#m.Equation(a['LHC'] >= 0.256*1.2)

# cell volume is determined by beta and the cytoplasmic 
# membrane surface. The volume is a constant, bot not the surface
#m.Equation(beta*(spA['cpm']*c['cpm']+spA['cpmP']*c['cpmP'])) == 1)


# SOLVING ----------------------------------------------------------
#
# solving maximizing specific growth rate; 
# objective is always minimized, so that we have 
# to state -1*obj to maximize it
m.Obj(-mu)
m.solve()


# COLLECTING RESULTS -----------------------------------------------
#
# assign new variable with optimal proteome mass fraction alpha
# and steday state starting concentrations 
# that can be reused by dynamic model
a_optim = a
c_start = c


# collect results in pandas data frame and save
result = pd.DataFrame(m.load_results())
result.to_csv('/home/micha/Documents/Code/models/GEKKO/cell-economy-models/result_steady_state.csv')
#result.to_csv('/home/michael/Documents/SciLifeLab/Resources/Models/GEKKO/cyano/result_steady_state.csv')

