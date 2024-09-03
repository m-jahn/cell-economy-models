# ----------------------------------------------------------
# Examples to import and run various cellular economy models
# ----------------------------------------------------------
#
# Model 1: Synechocystis sp. PCC6803 (cyanobacterium, photosynthetic)
# ----------------------------------------------------------


# 1. import libraries and model(s)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from glob import glob
from models import synechocystis_steadystate
from models import synechocystis_dynamic


# 2. define initial parameters
sub = 100  # initial substrate concentration, CO2/HCO3-
Ki = 5000  # light inhibition constant for photosystems
mumax = 0.11  # maximum growth rate, used to calculate protein utilization
# optional list of concentration upper bounds
ub_pro = pd.Series(
    [1, 1, 1, 1, 1, 1], index=["LHC", "PSET", "CBM", "LPB", "RIB", "MAI"]
)
ub_met = pd.Series([90, 25, 25, 5, 1], index=["hvi", "atp", "nadph", "pre", "lip"])
ub_mem = pd.Series([1, 1], index=["cpm", "thy"])
ub = pd.concat([ub_pro, ub_met, ub_mem])
reserve = [0.0, 0.0, 0.0, 0.0, 0.0]  # fraction of enzyme reserve in total


# 3. define light conditions

# (A) light in % max intensity, log decrease
# light = 100.0/1.5**np.array(range(0,12))

# (B) light as step change
# light = np.array([5.0]*25+[50.0]*26)

# (C) light coming in pulses
# light = np.array([3.0]*12+[50.0]*6+[3.0]*12+[50.0]*6+[3.0]*13)

# (D) light as smooth day night cycle
# use sine function to simulate one full day at length 2*pi = 6.283,
# so 2 days is 4*pi, and step width = 4*pi/96,
# since sine(x) is between -1 and 1, we rescale by (sine(x)+1)*50 (0 to 100)
light = np.round((np.sin(np.arange(0, 4 * 3.1415, 4 * 3.1415 / 96)) + 1) * 50) + 1

# time as a function of light step number, in hours
time = np.arange(0, len(light) / 2, 0.5)

# output directory
outdir = "results/synechocystis/daynight/"


# 3. run model simulations
#    loop through different values of ribosome reserve
for i in [0.0, 0.05, 0.1, 0.15]:
    reserve[4] = i
    res = "{0:02.2f}".format(i)
    result_ss = synechocystis_steadystate.simulate(
        time, light, sub, ub, reserve, mumax, Ki, remote=True
    )
    result_ss.table.to_csv(outdir + "steady_state_RIB_" + res + ".csv")
    result_dy = synechocystis_dynamic.simulate(
        time, light, sub, reserve, mumax, Ki, a=result_ss.a, c=result_ss.c
    )
    result_dy.table.to_csv(outdir + "dynamic_RIB_" + res + ".csv")


# 4. import result tables
df_steadystate = []
df_dynamic = []

for file in glob(outdir + "steady*0.00.csv"):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["reserve_rib"] = re.findall("RIB_[0-9]+\\.[0-9]+", file)[0]
    df_steadystate.append(df)

for file in glob(outdir + "dynamic*.csv"):
    df = pd.read_csv(file)
    df["type"] = "dynamic"
    df["reserve_rib"] = re.findall("RIB_[0-9]+\\.[0-9]+", file)[0]
    df_dynamic.append(df)

df_combined = pd.concat(df_steadystate + df_dynamic, ignore_index=True)
df_light = df_combined.query("type=='steady_state'")
df_light.loc[0, "hv"] = 0.0
df_light.loc[df_light.shape[0] - 1, "hv"] = 0.0
df_light["hv"] = df_light["hv"] / 100


# 5. visualize results
# generalized plotting function
def subplots(rows=3, cols=2, pos=1, ylim=[0, 1], title="", xvar="time", yvar="mu"):
    plt.subplot(rows, cols, pos)
    plt.axis([0, 50, ylim[0], ylim[1]])
    plt.title(title, loc="left", fontsize=10)
    plt.fill_between(
        x=df_light["time"],
        y1=df_light["hv"] * ylim[1],
        color="grey",
        alpha=0.2,
        linewidth=0,
    )
    sns.lineplot(
        x=df_combined[xvar],
        y=df_combined[yvar],
        hue=df_combined["type"] + " " + df_combined["reserve_rib"],
    )
    plt.legend(title="", fontsize="6", loc="lower right")


# set seaborn style
sns.set_theme(style="whitegrid", font_scale=0.75)

# create canvas
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

# growth rate, biomass, ribosomes, central carbon metabolism,
# light harvesting complex, 
subplots(pos=1, ylim=[0, 0.15], title="growth rate", xvar="time", yvar="mu")
subplots(pos=2, ylim=[0, 100], title="biomass acquisition", xvar="time", yvar="bm")
subplots(pos=3, ylim=[0, 0.5], title="ribosome mass fraction", xvar="time", yvar="c_rib")
subplots(pos=4, ylim=[0, 0.5], title="carbon metabolism mass fraction", xvar="time", yvar="c_cbm") 
subplots(pos=5, ylim=[0, 0.5], title="light harvesting mass fraction", xvar="time", yvar="c_lhc")
subplots(pos=6, ylim=[0, 0.5], title="photosystems mass fraction", xvar="time", yvar="c_pset")

plt.savefig(outdir + "pyplot_result.png", dpi = 182)