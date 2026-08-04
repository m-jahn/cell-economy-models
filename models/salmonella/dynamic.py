#
# CELLULAR ECONOMY MODEL FOR SALMONELLA
# version: 1.0
# subversion: "flagella expression trade-offs, dynamic"
# initial setup: 2024-09-04
# author: Michael Jahn
# affiliation: Max PLanck Unit for the Science of Pathogens
# based on: Jahn et al., 2018, Sharma et al., 2019, Molenaar et al., 2009
# implemented in Python using the GEKKO optimization package
# characteristics: resource allocation model with motility
#
#
# LIBRARIES ------------------------------------------------------------

from gekko import GEKKO
import pandas as pd
import numpy as np
from models import common


# INITIALIZE DYNAMIC MODEL ----------------------------------------
def simulate(time, time_init, dist_init, c_init, c_ub, a_fla, kcat=None, Km=None, hc=None, c=None, remote=False):

    m = GEKKO(remote = remote)
    m.options.IMODE = 6  # dynamic optimization
    m.options.REDUCE = 1
    m.options.MAX_ITER = 2000
    m.options.RTOL = 1e-5
    m.options.OTOL = 1e-5
    m.options.SCALING = 1
    m.options.SOLVER = 3
    m.options.TIME_SHIFT = 1
    m.time = time

    # organize variables in sets to simplify indexing
    enz = ["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla"]     # enzymes
    exc = ["ex_e", "ex_a"]                                      # exchange reactions "ex_c", "ex_l",
    pro = enz + ["Oth"]                                         # proteins
    met = ["cin", "cpre", "aa", "lip", "e"]                     # metabolites
    mem = ["cpm"]                                               # membrane compartment(s)
    memP = ["Tra", "Etc", "Fla"]                                # membrane located proteins
    cytP = ["Cbn", "Aab", "Rib", "Lpb", "Oth"]                  # cytoplasm located proteins

    # PARAMETERS --------------------------------------------------------
    #
    # NOTE: details on all parameter estimates in table 'parameters.csv'
    #
    # enzyme kinetic parameters as pandas series
    # kcat [molec/s], Km [mM], Hill coefficient [dimensionless]
    if kcat is None:
        kcat = pd.Series([200, 500, 100, 10, 22, 20, 250 * 100], index = enz)
    if Km is None:
        Km = pd.Series([20, 0.05, 0.03, 1, 1, 0.5, 1.0], index = enz)
    if hc is None:
        hc = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index = enz)

    # protein size [1000 aa]
    pro_size = pd.Series([1, 10, 10, 5, 7.5, 2, 7.2e3, 0.35], index = pro)

    # reaction stoichiometry matrix of met x enz
    stoich = pd.DataFrame([
        # Tra  Cbn  Etc  Aab  Rib  Lpb  Fla  ex_e ex_a   #
        [ 1,  -1,   0,   0,   0,   0,   0,   0,   0 ],   # cin
        [ 0,   2,  -1,  -2,   0,  -8,   0,   0,   0 ],   # cpre
        [ 0,   0,   0,   1,  -1,   0,   0,   0,  -1 ],   # aa
        [ 0,   0,   0,   0,   0,   1,   0,   0,   0 ],   # lip
        [ 0,   2,  30,  -2,  -3,  -8,  -1,  -1,   0 ]],  # e
        index = met,
        columns = enz + exc)

    # VARIABLES --------------------------------------------------------
    #
    # length of the cell [µm]
    length = m.Var(value=3, lb=2, ub=5, name = "length")

    # radius of the cell [µm]
    radius = m.Var(value=0.5, lb=0.25, ub=0.75, name = "radius")

    # volume of cylindrical cell (cylinder) [µm^3]
    volume = m.Var(value=8, lb=0, ub=10, name = "volume")

    # surface area of cylindrical cell  [µm^2]
    surface = m.Var(value=25, lb=0, ub=100, name = "surface")

    # density of the cell in 1000 aa / µm3 (PMID: 31690234)
    density = m.Param(value=8e6, name = "density")

    # max swimming speed in optimal conditions [µm * s^-1]
    v_swimmax = m.Param(value = 35, name = "v_swimmax")

    # actual swimming speed [µm * s^-1]
    v_swim = m.Var(value = 10, lb = 0, ub = 35, name = "v_swim")

    # distance to the source of the substrate
    distance = m.Var(value = dist_init, lb = 0, ub = dist_init, name = "distance")

    # specific surface area of membrane located components [µm^2]
    spA = pd.Series([1e-5, 5e-5, 1e-4, 5e-4], index = mem + memP)

    # area of membrane proteins as fraction of total surface
    surface_pro = m.Var(value=0.5, lb=0, ub=1, name = "surface_pro")

    # area of membrane lipids as fraction of total surface
    surface_lip = m.Var(value=0.5, lb=0, ub=1, name = "surface_lip")

    # list of catalytic rates v for all enzymes
    v = pd.Series(
        [m.Var(value = 1e5, lb = 0, ub = 2e7, name = "v_" + i) for i in enz + exc],
        index = enz + exc)

    # list of alpha = fraction of ribosomes engaged in synthesis of protein
    a = pd.Series(
        [m.Var(value = 0.1, lb = 0, ub = 1, name = "a_" + i) for i in pro],
        index = pro)

    # list of concentration of all components (enzymes and metabolites)
    c = pd.Series(
        [m.Var(value = c[i] if c is not None else 1e3, lb = 0, ub = c_ub[i], name = "c_" + i) for i in pro + met + mem],
        index = pro + met + mem)

    # cex is (time and location dependent) substrate concentration [mM]
    cex = m.Var(value = 1, lb = 0, ub = 10, name = "cex")
    
    # growth rate as variable that is to be maximized [h^-1]
    mu = m.Var(value = 0.1, lb = 0.1, ub = 1.5, name = "mu")

    # fraction of utilized protein space (total aa content / allowed aa content)
    utilization = m.Var(value=0.1, lb=0, ub=1, name = "utilization")

    # EQUATIONS --------------------------------------------------------
    #
    # equations constrain the solution space using parameters;
    # they outline the topology of the model

    # alpha is fraction of ribosomes engaged in synthesis of protein x
    m.Equation(sum(a) == 1)

    # protein mass balance: left side, rate of ribosome dedicated to
    # synthesis of each protein, right side, growth rate times protein conc
    m.Equations([a[i] * v["Rib"] - mu * c[i] * pro_size[i] == 0 for i in pro])

    # metabolite mass balance: left side, production of metabolites by
    # the respective enzyme, right side, growth rate times metabolite conc
    m.Equations([sum(stoich.loc[i] * v) - mu * c[i] == 0 for i in met])

    # Michaelis-Menthen type enzyme kinetics (V in molec s^-1 enz^-1)
    m.Equation(v["Tra"] == kcat["Tra"]*c["Tra"]*cex**hc["Tra"]/(Km["Tra"]**hc["Tra"] + cex**hc["Tra"]))
    m.Equation(v["Cbn"] == kcat["Cbn"]*c["Cbn"]*c["cin"]**hc["Cbn"]/(Km["Cbn"]**hc["Cbn"] + c["cin"]**hc["Cbn"]))
    m.Equation(v["Etc"] == kcat["Etc"]*c["Etc"]*c["cpre"]**hc["Etc"]/(Km["Etc"]**hc["Etc"] + c["cpre"]**hc["Etc"]))
    m.Equation(v["Aab"] == kcat["Aab"]*c["Aab"]*c["cpre"]**hc["Aab"]/(Km["Aab"]**hc["Aab"] + c["cpre"]**hc["Aab"]))
    m.Equation(v["Rib"] == kcat["Rib"]*c["Rib"]*c["aa"]**hc["Rib"]/(Km["Rib"]**hc["Rib"] + c["aa"]**hc["Rib"]))
    m.Equation(v["Lpb"] == kcat["Lpb"]*c["Lpb"]*c["cpre"]**hc["Lpb"]/(Km["Lpb"]**hc["Lpb"] + c["cpre"]**hc["Lpb"]))
    m.Equation(v["Fla"] == kcat["Fla"]*c["Fla"])

    # CELLULAR CONSTRAINTS
    #
    # cell volume [µm3] of the rod is determined by length and radius
    # average: l = 2-5 µm, r = 0.25-0.75 µm (PMID: 25664339)
    m.Equation(np.pi * radius**2 * length == volume)

    # cell surface [µm2] of the rod is determined by length and radius
    m.Equation(2 * np.pi * radius * (radius + length) == surface)

    # fraction of utilized protein space (total aa content / allowed aa content)
    m.Equation(utilization == sum(c[cytP] * pro_size[cytP]) / (volume * density))

    # total intracellular protein mass is constrained by volume and density
    m.Equation(sum(c[cytP] * pro_size[cytP]) <= volume * density)

    # area of membrane proteins as fraction of total surface
    m.Equation(surface_pro == sum(c[memP] * spA[memP]) / surface)

    # area of membrane lipids as fraction of total surface
    m.Equation(surface_lip == sum(c[mem] * spA[mem]) / surface)

    # membrane composition is constrained by total membrane surface area
    m.Equation(sum(c[mem + memP] * spA) == surface)

    # membrane proteins shall not exceed a certain fraction of total membrane components
    m.Equation(surface_pro <= surface_lip)

    # lipid balance: lipids are sum of cytoplasmic and other membranes
    m.Equation(sum(c[mem]) == c["lip"])

    # fix the mass fraction of maintenance proteins (or others)
    m.Equation(a["Oth"] == 0.35)

    # force production of flagella
    m.Equation(a["Fla"] == a_fla)

    # swimming speed [µm s^-1], based on max speed with optimal number of flagella
    m.Equation(v_swim == v_swimmax * c["Fla"] / (c["Fla"] + 10))

    # distance to source [µm], calculated from swimming speed and drift efficiency
    m.Equation(distance.dt() == -v_swim * 0.01 * 3600)

    # external substrate concentration depends on distance to source
    m.Equation(cex == c_init * m.erfc( distance / (2 * m.sqrt(600 * time_init))))

    # RATE OF CHANGE CONSTRAINTS
    #
    # limit how fast cell dimensions can change in [µm/h], to prevent unrealistic step changes
    m.Equation(length.dt() <= mu * length * 0.5)
    m.Equation(length.dt() >= 0)
    m.Equation(radius.dt() <= mu * radius * 0.5)
    m.Equation(radius.dt() >= 0)

    # limit growth rate change to only increase when moving up gradient
    m.Equations([mu.dt() >= 0])

    # optional constraint on protein concentration changes based on biological limits, can make simulation more stable
    # lower bound: proteins can decrease at most by growth-related dilution
    # m.Equations([c[i].dt() >= -mu * c[i] for i in pro])


    # SOLVING ----------------------------------------------------------
    #
    # objective: maximize specific growth rate
    m.Maximize(mu)
    m.solve()

    # convert rates from [1000 s^-1] to [h^-1]
    result = common.result("steady_state", m, v, a, c, c[pro])
    rate_cols = ["v_" + i.lower() for i in enz + exc]
    result.table[["mu"] + rate_cols] = result.table[["mu"] + rate_cols].apply(lambda x: x * 3.6)

    # return results
    return(result)
