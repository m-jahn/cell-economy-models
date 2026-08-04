import pandas as pd  # type: ignore
import numpy as np  # type: ignore
import matplotlib.pyplot as plt
import seaborn as sns
import math


# define class result where model results are collected
class result:
    def __init__(self, name, model, v, a, c, u):
        self.name = name
        self.model = model
        self.table = pd.DataFrame(model.load_results())
        self.v = v
        self.a = a
        self.c = c
        self.u = u


# function to randomize parameters within defined boundaries
def randomize(n, lb=0.5, ub=2.0):
    res = []
    for i in n:
        if i > 10:
            res += [float(np.random.randint(i * lb, i * ub, 1)[0])]
        else:
            res += [float(np.round(np.random.uniform(i * lb, i * ub), 2))]
    return res


# generalized plotting function
def subplots(
    df, xvar="time", yvar="mu", rows=4, cols=2, pos=1, ylim=None, title="", cex=False
):
    if xvar not in df.columns or yvar not in df.columns:
        return None
    df_sub = df[df["iteration"] == list(df["iteration"])[0]]
    ylim = [0, max(df[yvar]) * 1.2] if not ylim else ylim
    palette = sns.color_palette("YlOrBr", n_colors=len(set(df["iteration"])))
    plt.subplot(rows, cols, pos)
    plt.axis([0, max(df[xvar]), ylim[0], ylim[1]])
    plt.title(title, loc="left", fontsize=10)
    if cex:
        plt.fill_between(
            x=df_sub["time"],
            y1=ylim[0] + df_sub["cex"] / max(df_sub["cex"]) * np.diff(ylim),
            color="grey",
            alpha=0.2,
            linewidth=0,
        )
    sns.lineplot(x=df[xvar], y=df[yvar], hue=df["iteration"], palette=palette)
    plt.legend(title="", fontsize="6", loc="lower right")
    plt.grid(axis="both")


# plot growth rate, relative enzyme concentrations
def plot_enzymes(df, outdir):
    sns.set_theme(style="ticks", font_scale=0.75)
    # figure setup
    plt.figure(figsize=[8, 8])
    plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)
    # subplots
    subplots(df, xvar="time", yvar="mu", pos=1, title="growth rate")
    subplots(df, xvar="time", yvar="a_tra", pos=2, title="carbon transport")
    subplots(df, xvar="time", yvar="a_cbn", pos=3, title="carbon metabolism")
    subplots(df, xvar="time", yvar="a_etc", pos=4, title="electron transport chain")
    subplots(df, xvar="time", yvar="a_aab", pos=5, title="amino acid biosynthesis")
    subplots(df, xvar="time", yvar="a_rib", pos=6, title="ribosomes")
    subplots(df, xvar="time", yvar="a_lpb", pos=7, title="lipid biosynthesis")
    subplots(df, xvar="time", yvar="a_fla", pos=8, title="flagella biosynthesis")
    # export
    plt.savefig(outdir + "enzymes.png", dpi=182)
    plt.savefig(outdir + "enzymes.svg")


# plot  physicochemical properties
def plot_properties(df, outdir):
    sns.set_theme(style="ticks", font_scale=0.75)
    # figure setup
    plt.figure(figsize=[8, 8])
    plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)
    # subplots
    subplots(df, xvar="time", yvar="length", pos=1, title="length [µm]")
    subplots(df, xvar="time", yvar="radius", pos=2, title="radius [µm]")
    subplots(df, xvar="time", yvar="surface", pos=3, title="surface [µm^2]")
    subplots(df, xvar="time", yvar="volume", pos=4, title="volume [µm^3]")
    subplots(df, xvar="time", yvar="surface_pro", pos=5, title="area [mem_pro]")
    subplots(df, xvar="time", yvar="surface_lip", pos=6, title="area [mem lip]")
    subplots(df, xvar="time", yvar="utilization", pos=7, title="utilization")
    subplots(df, xvar="time", yvar="distance", pos=8, title="distance [µm]")
    # export
    plt.savefig(outdir + "properties.png", dpi=182)
    plt.savefig(outdir + "properties.svg")


def plot_rates(df, outdir):
    sns.set_theme(style="ticks", font_scale=0.75)
    # figure setup
    plt.figure(figsize=[8, 8])
    plt.subplots_adjust(wspace=0.5, hspace=0.7, top=0.925, bottom=0.075)
    # subplots
    subplots(df, xvar="time", yvar="v_tra", pos=1, title="V carbon transport")
    subplots(df, xvar="time", yvar="v_cbn", pos=2, title="V carbon metabolism")
    subplots(df, xvar="time", yvar="v_etc", pos=3, title="V ETC")
    subplots(df, xvar="time", yvar="v_aab", pos=4, title="V amino acid biosyn.")
    subplots(df, xvar="time", yvar="v_rib", pos=5, title="V ribosomes")
    subplots(df, xvar="time", yvar="v_lpb", pos=6, title="V lipid biosyn.")
    subplots(df, xvar="time", yvar="v_fla", pos=7, title="V flagellum")
    subplots(df, xvar="time", yvar="v_swim", pos=8, title="V swim [µm / s]")
    # export
    plt.savefig(outdir + "rates.png", dpi=182)
    plt.savefig(outdir + "rates.svg")


# Diffusion model for a simple source boundary with constant concentration C0
# and semi-infinite boundary (no reflection)
# reference: Crank, J. The Mathematics of Diffusion; Oxford University Press: New York, NY, USA, 1979.
# See Section 2.4.2, page 21
# Model: C(x,t) = C0 * erfc( x / (2 * sqrt(D * t)) ), with
#   C = concentration [mM] at time t [s] and location x [µm]
#   C0 = initial concentration [mM]
#   erfc = error function complement, erfc z = 1 - erf z
#   D = diffusion coefficient of glucose in water = 600 µm^2 / s (Bionumbers ID:104089)
def diffusion_model(x, t=3600, D=600, C0=5.0):
    C = C0 * math.erfc(x / (2 * math.sqrt(D * t)))
    return C
