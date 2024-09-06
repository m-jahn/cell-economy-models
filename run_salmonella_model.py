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
from glob import glob
from models.salmonella import steadystate


# 2. define initial parameters
# ----------------------------
remote = True                                        # remote or local solver
c_ex = [0, 0.75, 1.0, 5.0, 10.0, 50.0, 100.0]        # substrate concentration [µM]
time = np.arange(0, len(c_ex), 1)                    # time in [s]

# upper boundaries improve solve speed and success rate
c_ub_pro = pd.Series([5e5, 5e5, 5e5, 5e5, 5e5, 5e5, 1e0, 1e6], index=["Tra", "Cbn", "Etc", "Aab", "Rib", "Lpb", "Fla", "Oth"])
c_ub_met = pd.Series([5e5, 5e6, 1e5, 1e5, 1e6], index=["cin", "cpre", "aa", "lip", "e"])
c_ub_mem = pd.Series([1e5], index=["cpm"])
c_ub = pd.concat([c_ub_pro, c_ub_met, c_ub_mem])

# output directory
outdir = "results/salmonella/clim/"


# 3. run model simulations
# ------------------------
#    loop through different values of a variable
for i in [0.0]:
    iteration = "{0:02.2f}".format(i)
    result_ss = steadystate.simulate(time, c_ex, c_ub, remote)
    result_ss.table.to_csv(outdir + "steady_state_iter" + iteration + ".csv")


# 4. import result tables
# -----------------------
df_steadystate = []

for file in glob(outdir + "steady*.csv"):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("iter[0-9]+\\.[0-9]+", file)[0]
    df = df.query("time != 0")
    df_steadystate.append(df)

# combine into one df
df_combined = pd.concat(df_steadystate, ignore_index=True)

# convert rates from s^-1 to h^-1
rate_cols = [i for i in list(df_combined.columns) if i.startswith("v_")]
df_combined[["mu"] + rate_cols] = df_combined[["mu"] + rate_cols].apply(lambda x: x * 3.6e3)


# 5. visualize results
# --------------------
# generalized plotting function
def subplots(df, xvar="time", yvar="mu", rows=4, cols=2, pos=1, ylim=[0, 1], title=""):
    plt.subplot(rows, cols, pos)
    plt.axis([0, df.shape[0], ylim[0], ylim[1]])
    plt.title(title, loc="left", fontsize=10)
    plt.fill_between(
        x=df["time"],
        y1=df["cex"] / max(df["cex"])* ylim[1],
        color="grey",
        alpha=0.2,
        linewidth=0,
    )
    sns.lineplot(
        x=df[xvar],
        y=df[yvar],
        hue=df["type"],
    )
    plt.legend(title="", fontsize="6", loc="lower right")


# set seaborn style
sns.set_theme(style="whitegrid", font_scale=0.75)

# create canvas
plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

# growth rate, biomass, enzyme concentrations
subplots(df_combined, xvar="time", yvar="mu", pos=1, ylim=[0, 1.2], title="growth rate")
subplots(df_combined, xvar="time", yvar="bm", pos=2, ylim=[0, 1.2], title="biomass acquisition")
subplots(df_combined, xvar="time", yvar="a_rib", pos=3, ylim=[0, 0.5], title="ribosomes")
subplots(df_combined, xvar="time", yvar="a_cbn", pos=4, ylim=[0, 0.5], title="carbon metabolism")
subplots(df_combined, xvar="time", yvar="a_aab", pos=5, ylim=[0, 0.5], title="amino acid biosynthesis")
subplots(df_combined, xvar="time", yvar="a_etc", pos=6, ylim=[0, 0.5], title="electron transport chain")
subplots(df_combined, xvar="time", yvar="a_tra", pos=7, ylim=[0, 0.5], title="carbon transport")
subplots(df_combined, xvar="time", yvar="a_fla", pos=8, ylim=[0, 0.5], title="flagellum")

plt.savefig(outdir + "enzymes.png", dpi = 182)


# 6. plot physicochemical properties
# ----------------------------------
df_combined["density"] = df_combined["density"]/1e9

plt.figure(figsize=[8, 8])
plt.subplots_adjust(wspace=0.5, hspace=0.5)

# growth rate, biomass, enzyme concentrations
subplots(df_combined, xvar="time", yvar="length", pos=1, ylim=[1, 5], title="length [µm]")
subplots(df_combined, xvar="time", yvar="radius", pos=2, ylim=[0.2, 0.7], title="radius [µm]")
subplots(df_combined, xvar="time", yvar="surface", pos=3, ylim=[3, 6], title="surface [µm^2]")
subplots(df_combined, xvar="time", yvar="volume", pos=4, ylim=[0.3, 1.0], title="volume [µm^3]")
subplots(df_combined, xvar="time", yvar="density", pos=5, ylim=[0, 10], title="density [10^9 aa µm^3]")
#subplots(df_combined, xvar="time", yvar="", pos=6, ylim=[0, 0.5], title="")

plt.savefig(outdir + "properties.png", dpi = 182)

# special diagnostics:
# total aa content / allowed aa content
#sum(result.table.loc[2, ["c_tra", "c_cbn", "c_etc", "c_aab", "c_rib", "c_lpb", "c_fla", "c_oth"]] * pro_size.to_list())/(density[1] * volume[1]) * 100

# area of membrane proteins as fraction of total surface
#sum(result.table.loc[2, ["c_tra", "c_etc", "c_fla"]] * spA.to_list()[1:]) / result.table.loc[2, ["surface"]] * 100

# area of membrane lipids as fraction of total surface
#sum(result.table.loc[2, ["c_lip"]] * spA.to_list()[0])/ result.table.loc[2, ["surface"]] * 100

