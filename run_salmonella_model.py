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
from models.salmonella import dynamic
from models import common


# 2. define initial parameters
# ----------------------------
importlib.reload(steadystate)
remote = True
c_ex = np.round(2 ** np.arange(-6, 1, 0.5), 3)
time = np.arange(0, len(c_ex), 1)
max_retries = 3

# define sets
enz = ["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla"]
pro = enz + ["Oth"]
met = ["cin", "cpre", "aa", "lip", "e"]
mem = ["cpm"]

# upper boundaries
c_ub_pro = pd.Series([1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e3, 1e7], index=pro)
c_ub_met = pd.Series([1e5, 1e5, 1e5, 1e6, 1e6], index=met)
c_ub_mem = pd.Series([1e6], index=mem)
c_ub = pd.concat([c_ub_pro, c_ub_met, c_ub_mem])


# 3. run parameter sampling
# -------------------------
# perform parameter sampling to find stable sets and
# improve solver performance (typical problem is over-constrainment)
n_iterations = 0
n_solves = 0
outdir = "results/salmonella/sampling/"

while n_solves < 1 and n_iterations <= 10:
    n_iterations += 1
    iter = "{0:03d}".format(n_iterations)  # "{0:02.2f}".format(0.0)
    kcat = pd.Series(common.randomize([200, 500, 100, 10, 22, 20, 4e4]), index=enz)
    Km = pd.Series(common.randomize([20, 0.05, 0.03, 1, 1, 0.5, 1.0]), index=enz)
    hc = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index=enz)
    try:
        for a_fla in np.arange(0, 0.05, 0.01):
            result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
            result_ss.table.to_csv(outdir + "steady_state_iter_" + iter + "_flag_" + str(a_fla) + ".csv")
        kinetic_params = pd.DataFrame({"kcat": kcat, "Km": Km, "hc": hc})
        kinetic_params.to_csv(outdir + "kinetic_params_iter_" + iter + ".csv")
        n_solves += 1
    except:
        print("\n-------\nmodel not solvable, trying next parameter set")


# 4. run steady state model simulations
# -------------------------------------
#
# import desired parameter set
df_top_params = pd.read_csv("results/salmonella/sampling/top/kinetic_params_iter_01.csv", index_col=0)
kcat = df_top_params.kcat
Km = df_top_params.Km
hc = df_top_params.hc


# 4.1 simulate substrate limitation with different amount of flagella
outdir = "results/salmonella/c_limitation/"
retries = 0
for a_fla in np.arange(0, 0.06, 0.01):
    while retries <= max_retries:
        try:
            result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
            result_ss.table.to_csv(outdir + "steady_state_flag_" + "{0:02.2f}".format(a_fla) + ".csv")
            break
        except:
            print(f"\n---\nmodel not solvable, varying parameter (retry {retries})")
            kcat["Fla"] = kcat["Fla"] + np.random.normal(1) / 100
            retries += 1
    retries = 0


# 4.2 simulate substrate limitation with and without ATP cost for flagella (adjust stoich matrix)
outdir = "results/salmonella/rotation/"
c_ex = np.round(2 ** np.arange(-3, 0, 0.5), 3)
time = np.arange(0, len(c_ex), 1)
kcat["Fla"] = 400 * 100
retries = 0
for a_fla in np.arange(0, 0.06, 0.01):
    while retries <= max_retries:
        try:
            result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
            result_ss.table.to_csv(outdir + "steady_state_flag_" + "{0:02.2f}".format(a_fla) + ".csv")
            break
        except:
            print("\n-------\nmodel not solvable, varying parameter")
            kcat["Fla"] = kcat["Fla"] + np.random.normal(1) / 100
            retries += 1
    retries = 0


# 4.3 import result tables
df_steadystate = []

for file in sorted(glob(outdir + "steady_state*.csv")):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("(flag_[0-9]+\\.[0-9]+(\\_no_ATP)?)", file)[0][0]
    df = df.query("time != 0")
    df_steadystate.append(df)

# combine into one df
df_steadystate = pd.concat(df_steadystate, ignore_index=True)


# 4.4 visualize results
# set seaborn style
sns.set_theme(style="ticks", font_scale=0.75)


# 4.4.1 growth rate, relative enzyme concentrations
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_steadystate, xvar="time", yvar="mu", pos=1, ylim=[0, 1.0], title="growth rate")
common.subplots(df_steadystate, xvar="time", yvar="a_tra", pos=2, ylim=[0, 0.5], title="carbon transport")
common.subplots(df_steadystate, xvar="time", yvar="a_cbn", pos=3,ylim=[0, 0.01], title="carbon metabolism")
common.subplots(df_steadystate, xvar="time", yvar="a_etc", pos=4, ylim=[0, 0.01], title="electron transport chain")
common.subplots(df_steadystate, xvar="time", yvar="a_aab", pos=5, ylim=[0, 0.5], title="amino acid biosynthesis")
common.subplots(df_steadystate, xvar="time", yvar="a_rib", pos=6, ylim=[0, 0.1], title="ribosomes")
common.subplots(df_steadystate, xvar="time", yvar="a_lpb", pos=7, ylim=[0, 0.01], title="lipid biosynthesis")
common.subplots(df_steadystate, xvar="time", yvar="a_fla", pos=8, ylim=[0, 0.5], title="flagella biosynthesis")

plt.savefig(outdir + "enzymes.png", dpi=182)
plt.savefig(outdir + "enzymes.svg")


# 4.4.2  physicochemical properties
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_steadystate, xvar="time", yvar="length", pos=1, ylim=[1, 6], title="length [µm]")
common.subplots(df_steadystate, xvar="time", yvar="radius", pos=2, ylim=[0.2, 0.8], title="radius [µm]")
common.subplots(df_steadystate, xvar="time", yvar="surface", pos=3, ylim=[0, 35], title="surface [µm^2]")
common.subplots(df_steadystate, xvar="time", yvar="volume", pos=4, ylim=[0, 10.0], title="volume [µm^3]")
common.subplots(df_steadystate, xvar="time", yvar="surface_pro", pos=5, ylim=[0, 1.0], title="relative area of mem proteins")
common.subplots(df_steadystate, xvar="time", yvar="surface_lip", pos=6, ylim=[0, 1.0], title="relative area of mem lipids")
common.subplots(df_steadystate, xvar="time", yvar="utilization", pos=7, ylim=[0, 1.1], title="utilization")
common.subplots(df_steadystate, xvar="time", yvar="distance", pos=8, ylim=[0, 10000], title="distance to C source [µm]")

plt.savefig(outdir + "properties.png", dpi=182)
plt.savefig(outdir + "properties.svg")


# 4.4.3 enzymatic rates
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_steadystate, xvar="time", yvar="v_tra", pos=1, ylim=[0, 1e7], title="V carbon transport")
common.subplots(df_steadystate, xvar="time", yvar="v_cbn", pos=2, ylim=[0, 1e7], title="V carbon metabolism")
common.subplots(df_steadystate, xvar="time", yvar="v_etc", pos=3, ylim=[0, 1e7], title="V electron transport chain")
common.subplots(df_steadystate, xvar="time", yvar="v_aab", pos=4, ylim=[0, 1e7], title="V amino acid biosynthesis")
common.subplots(df_steadystate, xvar="time", yvar="v_rib", pos=5, ylim=[0, 1e7], title="V ribosomes")
common.subplots(df_steadystate, xvar="time", yvar="v_lpb", pos=6, ylim=[0, 1e7], title="V lipid biosynthesis")
common.subplots(df_steadystate, xvar="time", yvar="v_fla", pos=7, ylim=[0, 2e7], title="V flagellum")
common.subplots(df_steadystate, xvar="time", yvar="v_swim", pos=8, ylim=[0, 50], title="V swim [µm / s]")

plt.savefig(outdir + "rates.png", dpi=182)
plt.savefig(outdir + "rates.svg")


# 4.4.4 flagella with/without energy cost
df_cost = df_steadystate[df_steadystate["time"] == 3.0]
df_cost["a_fla"] = df_cost["a_fla"] * 100
df_cost["energy cost"] = df_cost["iteration"].str.contains("no_ATP")
df_cost["energy cost"] = df_cost["energy cost"].apply(
    lambda x: "without ATP cost" if x else "with ATP cost"
)
plt.figure(figsize=[8, 3.5])
plt.subplot(1, 2, 1)
plt.subplots_adjust(bottom=0.15)
plt.title("growth rate with increasing flagella", loc="left", fontsize=10)
ax = sns.barplot(
    x=df_cost["a_fla"].astype(str),
    y=df_cost["mu"],
    hue=df_cost["energy cost"],
    palette=sns.color_palette("flare", n_colors=2),
)
ax.set(xlabel="% protein to flagellum", ylabel="growth rate [h^-1]")
plt.legend(title="", fontsize="6", loc="lower right")
plt.grid(axis="both")

plt.subplot(1, 2, 2)
plt.title("flagellar activity with increasing flagella", loc="left", fontsize=10)
ax = sns.barplot(
    x=df_cost["a_fla"].astype(str),
    y=df_cost["v_fla"],
    hue=df_cost["energy cost"],
    palette=sns.color_palette("flare", n_colors=2),
)
ax.set(xlabel="% protein to flagellum", ylabel="V flagellum")
plt.legend(title="", fontsize="6", loc="lower right")
plt.grid(axis="both")
plt.savefig(outdir + "energy_vs_protein_cost.png", dpi=182)
plt.savefig(outdir + "energy_vs_protein_cost.svg")


# 5. run dynamic model simulations
# -------------------------------------
#
# 5.1 simulate swimming at variable speed, depending on number of flagella
importlib.reload(dynamic)
outdir = "results/salmonella/swimming/"
c_init = 5.0 # [mM]
dist_init = 8000 # [µm]
time_init = 3 * 3600 # [sec] only relevant for substrate gradient
time =  np.concatenate([[0, 0.1], np.arange(0.5, 10, 0.5)]) # [h]
retries = 0
# common.diffusion_model(x = dist_init, t = time_init, D = 600, C0 = c_init)
for a_fla in [0.001, 0.021, 0.051]: #np.arange(0.00, 0.06, 0.01):
    while retries <= max_retries:
        try:
            result_ss = dynamic.simulate(time, time_init, dist_init, c_init, c_ub, a_fla, kcat, Km, hc, remote)
            result_ss.table.to_csv(outdir + "dynamic_flag_" + "{0:02.2f}".format(a_fla) + ".csv")
            break
        except:
            print("\n-------\nmodel not solvable, varying parameter")
            c_init = c_init + np.random.normal(1) / 100
            retries += 1
    retries = 0


# 5.2 import result tables
df_dynamic = []

for file in sorted(glob(outdir + "dynamic*.csv")):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("(flag_[0-9]+\\.[0-9]+(\\_no_ATP)?)", file)[0][0]
    df = df.query("time != 0")
    df_dynamic.append(df)

# combine into one df
df_combined = pd.concat(df_dynamic, ignore_index=True)
df_combined = df_combined[df_combined["iteration"].str.contains("0\\.0[0-5]")]


# 5.3 visualize results
# set seaborn style
sns.set_theme(style="ticks", font_scale=0.75)


# 5.3.1 growth rate, biomass, enzyme concentrations
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_combined, xvar="time", yvar="mu", pos=1, ylim=[0, 1.0], title="growth rate")
common.subplots(df_combined, xvar="time", yvar="a_tra", pos=2, ylim=[0, 0.5], title="carbon transport")
common.subplots(df_combined, xvar="time", yvar="a_cbn", pos=3,ylim=[0, 0.01], title="carbon metabolism")
common.subplots(df_combined, xvar="time", yvar="a_etc", pos=4, ylim=[0, 0.01], title="electron transport chain")
common.subplots(df_combined, xvar="time", yvar="a_aab", pos=5, ylim=[0, 0.5], title="amino acid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="a_rib", pos=6, ylim=[0, 0.1], title="ribosomes")
common.subplots(df_combined, xvar="time", yvar="a_lpb", pos=7, ylim=[0, 0.01], title="lipid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="a_fla", pos=8, ylim=[0, 0.5], title="flagella biosynthesis")

plt.savefig(outdir + "enzymes.png", dpi=182)
plt.savefig(outdir + "enzymes.svg")


# 5.3.2  physicochemical properties
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_combined, xvar="time", yvar="length", pos=1, ylim=[1, 6], title="length [µm]")
common.subplots(df_combined, xvar="time", yvar="radius", pos=2, ylim=[0.2, 0.8], title="radius [µm]")
common.subplots(df_combined, xvar="time", yvar="surface", pos=3, ylim=[0, 35], title="surface [µm^2]")
common.subplots(df_combined, xvar="time", yvar="volume", pos=4, ylim=[0, 10.0], title="volume [µm^3]")
common.subplots(df_combined, xvar="time", yvar="surface_pro", pos=5, ylim=[0, 1.0], title="relative area of mem proteins")
common.subplots(df_combined, xvar="time", yvar="surface_lip", pos=6, ylim=[0, 1.0], title="relative area of mem lipids")
common.subplots(df_combined, xvar="time", yvar="utilization", pos=7, ylim=[0, 1.1], title="utilization")
common.subplots(df_combined, xvar="time", yvar="distance", pos=8, ylim=[0, 10000], title="distance to C source [µm]")

plt.savefig(outdir + "properties.png", dpi=182)
plt.savefig(outdir + "properties.svg")


# 5.3.3 enzymatic rates
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)

common.subplots(df_combined, xvar="time", yvar="v_tra", pos=1, ylim=[0, 1e7], title="V carbon transport")
common.subplots(df_combined, xvar="time", yvar="v_cbn", pos=2, ylim=[0, 1e7], title="V carbon metabolism")
common.subplots(df_combined, xvar="time", yvar="v_etc", pos=3, ylim=[0, 1e7], title="V electron transport chain")
common.subplots(df_combined, xvar="time", yvar="v_aab", pos=4, ylim=[0, 1e7], title="V amino acid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="v_rib", pos=5, ylim=[0, 1e7], title="V ribosomes")
common.subplots(df_combined, xvar="time", yvar="v_lpb", pos=6, ylim=[0, 1e7], title="V lipid biosynthesis")
common.subplots(df_combined, xvar="time", yvar="v_fla", pos=7, ylim=[0, 2e7], title="V flagellum")
common.subplots(df_combined, xvar="time", yvar="v_swim", pos=8, ylim=[0, 50], title="V swim [µm / s]")

plt.savefig(outdir + "rates.png", dpi=182)
plt.savefig(outdir + "rates.svg")
