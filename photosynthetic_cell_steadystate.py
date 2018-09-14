#
# PHOTOAUTOTROPHIC CELL MODEL
# version: 1.0
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
import re
import math
import numpy as np
import pandas as pd


# define function that runs one steady state optimization of
# the model with the given parameters
def optim(substrate, light, photoinhib=4351):
    
    # initialize model
    m = GEKKO(remote=True, server='http://xps.apmonitor.com') # alternative: server='http://xps.apmonitor.com'
    m.options.IMODE=3
    m.options.REDUCE=1
    
    
    # PARAMETERS -------------------------------------------------------
    
    # organize variables in sets to simplify indexing
    enz = ['LHC', 'PSET', 'CBM', 'LPB', 'RIB']      # enzymes
    pro = enz+['MAI']                               # proteins
    metab = ['hvi', 'atp', 'nadph', 'pre', 'lip']   # metabolites
    mem = ['cpm', 'thy']                            # membrane lipids
    thyP = ['LHC', 'PSET']                          # thylakoid membrane located proteins
    #cmpP = []                                      # cytoplasmic membrane located proteins
    
    
    
    # Gekko parameters are usually fixed single values or series
    # enzyme kinetic paraemters as pandas series 
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
        index=metab,
        columns=enz)
    
    
    # define constants. Usually single fixed values (integer, float)
    sub = m.Const(substrate, 'initial substrate concentration, CO2/HCO3-')
    hv  = m.Const(light, 'extracellular light irradiation')
    Ki = m.Const(photoinhib, 'light inhibition constant for photosystems')
    
    
    # VARIABLES --------------------------------------------------------
    
    # calculated by solver to meet the constraints of equations and params
    # format: v = m.Var(1, [lb], [ub], [integer], [name]):
    
    # list of catalytic rates v for all enzymes
    v = pd.Series(
        [m.Var(value=x, lb=0, ub=100) for x in [1]*len(enz)],
        index=enz)
    
    # list of alpha=fraction of ribosomes engaged in synthesis of protein
    a = pd.Series(
        [m.Var(value=x, lb=0, ub=1) for x in [1]*len(pro)],
        index=pro)
    
    # list of concentration of all components (enzymes and metabolites)
    c = pd.Series(
        [m.Var(value=x, lb=0, ub=100) for x in [1]*len(pro+metab+mem)],
        index=pro+metab+mem)
    
    
    # growth rate as variable that becomes the objective function
    mu = m.Var(value=1)
    # volume-to-surface ratio, increases with sphericity of a cell
    #beta = m.Var(value=1, lb=0.1, ub=10)
    
    
    # EQUATIONS --------------------------------------------------------
    # equations constrain the solution space using parameters;
    # they outline the topology of the model
    
    # alpha is fraction of ribosomes engaged in synthesis of protein x
    m.Equation(sum(a) == 1)
    
    # protein mass balance: left side, rate of ribosome dedicated to 
    # synthesis of each protein, right side, growth rate times protein conc
    eq1 = [m.Equation(a[i]*v['RIB'] - mu*c[i] == 0) for i in pro]
    
    # metabolite mass balance: left side, production of metabolites by
    # the respective enzyme, right side, growth rate times metabolite conc
    eq2 = [m.Equation(sum(stoich.loc[i]*v) - mu*c[i] == 0) for i in metab]
    
    
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
    m.Equation(a['MAI'] == 0.3)
    
    # cell volume is determined by beta and the cytoplasmic 
    # membrane surface. The volume is a constant, bot not the surface
    #m.Equation(beta*(spA['cpm']*c['cpm']+spA['cpmP']*c['cpmP'])) == 1)
    
    
    # OBJECTIVE --------------------------------------------------------
    
    # specific growth rate; objective is always minimized, so that we have 
    # to state -1*obj to maximize it
    m.Obj(-1*mu)
    
    
    # SOLVING ----------------------------------------------------------
    # Solve simulation
    m.solve()
    
    # return results of optimization as pandas data frame
    # with 'long' column format suitable for R
    return pd.DataFrame([
        mu.value*(1+len(a)+len(c)+len(v)), 
        [sub.value]*(1+len(a)+len(c)+len(v)), 
        [hv.value]*(1+len(a)+len(c)+len(v)),
        ['mu']+['a']*len(a)+['c']*len(c)+['v']*len(v), 
        ['model']+a.index.tolist()+c.index.tolist()+v.index.tolist(),
        mu.value+[re.sub('\[|\]', '', str(i.value)) for i in pd.concat([a, c, v])]
    ]).T


# ITERATIVE SOLVING OF MODEL -------------------------------------------
#
# iterate over changing parameters and collect results in
# pandas data frame
result = pd.DataFrame([])
for i in [100/(2**x) for x in np.linspace(0,7,15)]: #
    result = result.append(
        optim(substrate=100, light=i, photoinhib=5000)
    )

# add some column names and export csv table
result.columns = ['mu', 'substrate', 'light', 'variable', 'component', 'concentration']
result['concentration'] = result['concentration'].astype('float')
result.to_csv('/home/michael/Documents/SciLifeLab/Resources/Models/GEKKO/cyano/result.csv')

