#
# CELLULAR ECONOMY MODEL FOR SALMONELLA
# version: 1.0
# subversion: "flagella expression trade-offs, steady state"
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


# INITIALIZE STEADY STATE MODEL ----------------------------------------
#def simulate(time, c_ex, c_ub, reserve, mumax, remote = False):

# STARTING PARAMETERS --> to be removed in the end
remote = True
cex = [0, 0.75, 1.0, 5.0, 10.0, 50.0, 100.0]        # [µM]
time = np.arange(0, len(cex), 1)                    #  [h]
c_ub_pro = pd.Series([5e5, 5e5, 5e5, 5e5, 5e5, 5e5, 1e0, 1e6], index=["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla", "Oth"])
c_ub_met = pd.Series([5e5, 1e6, 1e5, 1e4, 1e6], index=["cin", "cpre", "aa", "lip", "e"])
c_ub_mem = pd.Series([1e5], index=["cpm"])
c_ub = pd.concat([c_ub_pro, c_ub_met, c_ub_mem])
mumax = 1.0

m = GEKKO(remote = remote)
m.options.IMODE = 5
m.options.REDUCE = 1
m.options.MAX_ITER = 1000
m.time = time


# organize variables in sets to simplify indexing
enz = ["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla"]     # enzymes
pro = enz + ["Oth"]                                           # proteins
met = ["cin", "cpre", "aa", "lip", "e"]                     # metabolites
mem = ["cpm"]                                               # membrane compartment(s)
memP = ["Tra", "Etc", "Fla"]                                # membrane located proteins
cytP = ["Cbn", "Aab", "Rib", "Lpb", "Oth"]                  # cytoplasm located proteins


# PARAMETERS --------------------------------------------------------
#
# enzyme kinetic parameters as pandas series
# kcat [molec/s], Km [µM], Hill coefficient [dimensionless]
kcat = pd.Series([50, 10, 50, 10, 20, 50, 50], index = enz)
Km = pd.Series([25, 20, 20, 100, 15, 35, 10], index = enz)
hc = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index = enz)

# protein size [aa]; more details about size estimation in suppl. tables
pro_size = pd.Series([1500, 2000, 10000, 20000, 7500, 2000, 4e6, 1000], index = pro)

# protein reserve (inactive proteins) [molec]
reserve = pd.Series([0, 0, 0, 0, 0, 0, 0], index = enz)

# reaction stoichiometry matrix of met x enz
# more details on stoichiometry estimates in suppl. table
# 1 rotation of flagellum costs ~1200 protons equivalent to 400 ATP (PMID: 35881430)
stoich = pd.DataFrame([
    # Tra  Cbn  Etc  Aab  Rib  Lpb  Fla    #
    [ 1,  -1,   0,   0,   0,   0,   0 ],   # cin
    [ 0,   2,  -1,  -2,   0,  -1,   0 ],   # cpre
    [ 0,   0,   0,   1,  -1,   0,   0 ],   # aa
    [ 0,   0,   0,   0,   0,   1,   0 ],   # lip
    [ 0,   2,  30,  -1,  -3,  -1,   0 ]],  # e
    index = met,
    columns = enz)

# VARIABLES --------------------------------------------------------
#
# length of the cell [µm]
length = m.Var(value=3, lb=2, ub=5, name = "length")

# radius of the cell [µm]
radius = m.Var(value=0.5, lb=0.25, ub=0.75, name = "radius")

# volume of cylindrical cell (cylinder) [µm^3]
volume = m.Var(value=1, lb=0, ub=10, name = "volume")

# surface area of cylindrical cell  [µm^2]
surface = m.Var(value=10, lb=0, ub=100, name = "surface")

# density of the cell in aa / µm3 (PMID: 31690234)
density = m.Param(value=8e9, name = "density")

# specific surface area of membrane located components [µm^2]
# phospholipid = 1 nm^2 (1e-6), transporter = 50 nm^2 (5e-5), ETC complexes xx nm^2, flagellum = 500 nm^2 (5e-4)
spA = pd.Series([1e-3, 5e-5, 1e-4, 5e-4], index = mem + memP)

# list of catalytic rates v for all enzymes
v = pd.Series(
    [m.Var(value = 1, lb = 0, ub = 1e6, name = "v_" + i) for i in enz],
    index = enz)

# list of alpha = fraction of ribosomes engaged in synthesis of protein
a = pd.Series(
    [m.Var(value = 1, lb = 0, ub = 1, name = "a_" + i) for i in pro],
    index = pro)

# list of concentration of all components (enzymes and metabolites)
c = pd.Series(
    [m.Var(value = 1, lb = 0, ub = c_ub[i], name = "c_" + i) for i in pro + met + mem],
    index = pro + met + mem)

# cx is (time dependent) substrate concentration [µM]
cx = m.Param(value = cex, name = "cex")

# growth rate as variable that is to be maximized [h^-1]
mu = m.Var(value = 1, name = "mu")

# biomass accumulated over time with initial value [fold change]
bm = m.Var(value = 1, name = "bm")


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

# biomass accumulation over time
m.Equation(bm.dt() == mu*bm)

# Michaelis-Menthen type enzyme kinetics (V in molec s^-1 enz^-1)
m.Equation(v["Tra"] == kcat["Tra"]*c["Tra"]*cx**hc["Tra"]/(Km["Tra"]**hc["Tra"] + cx**hc["Tra"]))
m.Equation(v["Cbn"] == kcat["Cbn"]*c["Cbn"]*c["cin"]**hc["Cbn"]/(Km["Cbn"]**hc["Cbn"] + c["cin"]**hc["Cbn"]))
m.Equation(v["Etc"] == kcat["Etc"]*c["Etc"]*c["cpre"]**hc["Etc"]/(Km["Etc"]**hc["Etc"] + c["cpre"]**hc["Etc"]))
m.Equation(v["Aab"] == kcat["Aab"]*c["Aab"]*c["cpre"]**hc["Aab"]/(Km["Aab"]**hc["Aab"] + c["cpre"]**hc["Aab"]))
m.Equation(v["Rib"] == kcat["Rib"]*c["Rib"]*c["aa"]**hc["Rib"]/(Km["Rib"]**hc["Rib"] + c["aa"]**hc["Rib"]))
m.Equation(v["Lpb"] == kcat["Lpb"]*c["Lpb"]*c["cpre"]**hc["Lpb"]/(Km["Lpb"]**hc["Lpb"] + c["cpre"]**hc["Lpb"]))
m.Equation(v["Fla"] == kcat["Fla"]*c["Fla"]*c["cpre"]**hc["Fla"]/(Km["Fla"]**hc["Fla"] + c["cpre"]**hc["Fla"]))


# CELLULAR CONSTRAINTS
#
# cell volume [µm3] of the rod is determined by length and radius
# average: l = 2-5 µm, r = 0.25-0.75 µm (PMID: 25664339)
m.Equation(np.pi * radius**2 * length == volume)

# cell surface [µm2] of the rod is determined by length and radius
m.Equation(2 * np.pi * radius * (radius + length) == surface)

# total intracellular protein mass is constrained by volume and density
m.Equation(sum(c[cytP] * pro_size[cytP]) <= volume * density)

# membrane composition is constrained by total membrane surface area
m.Equation(sum(c[mem + memP] * spA) == surface)

# membrane proteins shall shall not exceed a certain fraction of total membrane components
m.Equation(sum(c[memP] * spA[memP]) <= sum(c[mem] * spA[mem]))

# lipid balance: lipids are sum of cytoplasmic and other membranes
m.Equation(sum(c[mem]) == c["lip"])

# fix the mass fraction of maintenance proteins (or others)
m.Equation(a["Oth"] == 0.25)


# SOLVING ----------------------------------------------------------
#
# solving maximizing specific growth rate;
# objective is always minimized, so that we have
# to state -1*obj to maximize it
m.Obj(-mu)
m.solve()

# collect results and return
#return(common.result("steady_state", m, v, a, c, u))
result = common.result("steady_state", m, v, a, c, a)
print(result.table[["time", "cex", "mu", "v_tra", "v_cbn", "v_etc", "v_aab", "v_rib", "v_lpb", "v_fla"]])
print(result.table[["time", "cex", "mu", "a_tra", "a_cbn", "a_etc", "a_aab", "a_rib", "a_lpb", "a_fla", "a_oth"]])
print(result.table[["time", "cex", "mu", "c_tra", "c_cbn", "c_etc", "c_aab", "c_rib", "c_lpb", "c_fla", "c_oth"]])
print(result.table[["time", "cex", "mu", "c_cin", "c_cpre", "c_aa", "c_lip", "c_e"]])
print(result.table[["length", "radius", "surface", "volume", "density"]])
result.table.to_csv("result.csv")

# special diagnostics:
# total aa content / allowed aa content
sum(result.table.loc[2, ["c_tra", "c_cbn", "c_etc", "c_aab", "c_rib", "c_lpb", "c_fla", "c_oth"]] * pro_size.to_list())/(density[1] * volume[1]) * 100

# area of membrane proteins as fraction of total surface
sum(result.table.loc[2, ["c_tra", "c_etc", "c_fla"]] * spA.to_list()[1:]) / result.table.loc[2, ["surface"]] * 100

# area of membrane lipids as fraction of total surface
sum(result.table.loc[2, ["c_lip"]] * spA.to_list()[0])/ result.table.loc[2, ["surface"]] * 100

# Retrieve shadow prices (Lagrange mutipliers)
# only for local solve with IPOPT!
# lam = np.loadtxt(m.path + '/apm_lam.txt')
# print(lam.shape)
