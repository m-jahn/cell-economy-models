import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


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
def subplots(df, xvar="time", yvar="mu", rows=4, cols=2, pos=1, ylim=[0, 1], title=""):
    df_sub = df[df["iteration"] == list(df["iteration"])[0]]
    plt.subplot(rows, cols, pos)
    plt.axis([0, max(df[xvar]), ylim[0], ylim[1]])
    plt.title(title, loc="left", fontsize=10)
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
    )
    plt.legend(title="", fontsize="6", loc="lower right")
    plt.grid(axis="both")
