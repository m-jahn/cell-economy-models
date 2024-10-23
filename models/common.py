import pandas as pd
import numpy as np
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
def subplots(df, xvar="time", yvar="mu", rows=4, cols=2, pos=1, ylim=[0, 1], title="", cex=False):
    df_sub = df[df["iteration"] == list(df["iteration"])[0]]
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
    sns.lineplot(
        x=df[xvar],
        y=df[yvar],
        hue=df["iteration"],
        palette=palette
    )
    plt.legend(title="", fontsize="6", loc="lower right")
    plt.grid(axis="both")


# Diffusion model for a simple source boundary with constant concentration C0
# and semi-infinite boundary (no reflection)
# reference: Crank, J. The Mathematics of Diffusion; Oxford University Press: New York, NY, USA, 1979. 
# See Section 2.4.2, page 21
# Model: C(x,t) = C0 * erfc( x / (2 * sqrt(D * t)) ), with
#   C = concentration [mM] at time t [s] and location x [µm]
#   C0 = initial concentration [mM]
#   erfc = error function complement, erfc z = 1 - erf z
#   D = diffusion coefficient of glucose in water = 600 µm^2 / s (Bionumbers ID:104089)
def diffusion_model(x, t=3600, D=600, C0=5):
    C = C0 * math.erfc( x / (2 * math.sqrt(D * t)))
    return(C)
