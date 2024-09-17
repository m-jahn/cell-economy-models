# ----------------------------------------------------------
# Run cell economy model for Salmonella typhimurium
# (pathogenic motile bacterium, heterotrophic)
# ----------------------------------------------------------


# 1. import libraries and model(s)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import importlib
from glob import glob
from models.salmonella import steadystate
from models import common


# 2. define initial parameters
# ----------------------------
importlib.reload(steadystate)
remote = True
c_ex = np.round(2 ** np.arange(-6, 1, 0.5), 3)
time = np.arange(0, len(c_ex), 1)

# upper boundaries improve solve speed and success rate
enz = ["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla"]
c_ub_pro = pd.Series([1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 0, 1e7], index=enz + ["Oth"])
c_ub_met = pd.Series([1e5, 1e5, 1e5, 1e6, 1e6], index=["cin", "cpre", "aa", "lip", "e"])
c_ub_mem = pd.Series([1e6], index=["cpm"])
c_ub = pd.concat([c_ub_pro, c_ub_met, c_ub_mem])


# 3. run parameter sampling
# -------------------------
# perform parameter sampling to find stable sets and
# improve solver performance (typical problem is over-constrainment)
n_iterations = 0
n_solves = 0
n_flag = 0
outdir = "results/salmonella/sampling/"


while n_solves < 1 or n_iterations <= 20:
    n_iterations += 1
    try:
        iter = "{0:03d}".format(n_iterations) # "{0:02.1f}".format(n_iterations)
        kcat = pd.Series(common.randomize([200, 500, 100, 10, 22, 20, 400]), index=enz)
        Km = pd.Series(common.randomize([20, 0.05, 0.03, 1, 1, 0.5, 1]), index=enz)
        hc = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index=enz)
        result_ss = steadystate.simulate(time, c_ex, c_ub, n_flag, kcat, Km, hc,remote)
        result_ss.table.to_csv(outdir + "steady_state_iter_" + iter + ".csv")
        kinetic_params = pd.DataFrame({"kcat" : kcat, "Km" : Km, "hc" : hc})
        kinetic_params.to_csv(outdir + "kinetic_params_iter_" + iter + ".csv")
        n_solves += 1
    except:
        print("\n-------\nmodel not solvable, trying next parameter set")


# 4. run model simulations
# ------------------------
#
# import desired parameter set
df_top_params = pd.read_csv("results/salmonella/sampling/top/kinetic_params_iter_019.csv", index_col=0)
kcat = df_top_params.kcat
Km = df_top_params.Km
hc = df_top_params.hc

# loop through different values of a variable
outdir = "results/salmonella/flagella/"
for n_flag in [0, 1, 2, 3, 4]:
    try:
        result_ss = steadystate.simulate(time, c_ex, c_ub, n_flag, kcat, Km, hc, remote)
        result_ss.table.to_csv(outdir + "steady_state_flag_" + str(n_flag) + ".csv")
    except:
        print("\n-------\nmodel not solvable, trying next parameter set")


# 5. import result tables
# -----------------------
df_steadystate = []

for file in glob(outdir + "steady*.csv"):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("flag_[0-9]+", file)[0]
    df = df.query("time != 0")
    df_steadystate.append(df)

# combine into one df
df_combined = pd.concat(df_steadystate, ignore_index=True)


# 6. visualize results
# --------------------
# set seaborn style
sns.set_theme(style="ticks", font_scale=0.75)

# create canvas
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

# 5.1 growth rate, biomass, enzyme concentrations
# -----------------------------------------------
common.subplots(df_combined, xvar="time", yvar="mu", pos=1, ylim=[0, 1.0], title="growth rate")
common.subplots(df_combined, xvar="time", yvar="bm", pos=2, ylim=[0, 100.0], title="biomass acquisition")
common.subplots(df_combined, xvar="time", yvar="a_tra", pos=3, ylim=[0, 0.5], title="carbon transport")
common.subplots(df_combined, xvar="time", yvar="a_cbn", pos=4,ylim=[0, 0.5], title="carbon metabolism")
common.subplots(df_combined, xvar="time", yvar="a_etc", pos=5, ylim=[0, 0.5], title="electron transport chain")
common.subplots(df_combined, xvar="time", yvar="a_aab", pos=6, ylim=[0, 0.5], title="amino acid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="a_rib", pos=7, ylim=[0, 0.5], title="ribosomes")
common.subplots(df_combined, xvar="time", yvar="a_lpb", pos=8, ylim=[0, 0.5], title="lipid biosynthesis")

plt.savefig(outdir + "enzymes.png", dpi=182)


# 5.2  physicochemical properties
# -------------------------------
df_combined["density"] = df_combined["density"] / 1e6

plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

common.subplots(df_combined, xvar="time", yvar="length", pos=1, ylim=[1, 6], title="length [µm]")
common.subplots(df_combined, xvar="time", yvar="radius", pos=2, ylim=[0.2, 0.8], title="radius [µm]")
common.subplots(df_combined, xvar="time", yvar="surface", pos=3, ylim=[0, 35], title="surface [µm^2]")
common.subplots(df_combined, xvar="time", yvar="volume", pos=4, ylim=[0, 10.0], title="volume [µm^3]")
common.subplots(df_combined, xvar="time", yvar="density", pos=5, ylim=[0, 10], title="density [10^9 aa µm^3]")
common.subplots(df_combined, xvar="time", yvar="utilization", pos=6, ylim=[0, 1.1], title="utilization")
common.subplots(df_combined, xvar="time", yvar="surface_pro", pos=7, ylim=[0, 1.0], title="relative area of mem proteins")
common.subplots(df_combined, xvar="time", yvar="surface_lip", pos=8, ylim=[0, 1.0], title="relative area of mem lipids")

plt.savefig(outdir + "properties.png", dpi=182)

# 5.3 enzymatic rates
# -------------------
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

common.subplots(df_combined, xvar="time", yvar="v_tra", pos=1, ylim=[0, 2e6], title="V carbon transport")
common.subplots(df_combined, xvar="time", yvar="v_cbn", pos=2, ylim=[0, 2e6], title="V carbon metabolism")
common.subplots(df_combined, xvar="time", yvar="v_etc", pos=3, ylim=[0, 2e6], title="V electron transport chain")
common.subplots(df_combined, xvar="time", yvar="v_aab", pos=4, ylim=[0, 2e6], title="V amino acid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="v_rib", pos=5, ylim=[0, 2e6], title="V ribosomes")
common.subplots(df_combined, xvar="time", yvar="v_lpb", pos=6, ylim=[0, 2e6], title="V lipid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="v_fla", pos=7, ylim=[0, 2e6], title="V flagellum")

plt.savefig(outdir + "rates.png", dpi=182)
