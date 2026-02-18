# ----------------------------------------------------------
# Run cell economy model for Salmonella typhimurium
# (pathogenic motile bacterium, heterotrophic)
# ----------------------------------------------------------


# 1. import libraries and model(s)
import pandas as pd # type: ignore
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore
import os
import re
import importlib
from math import log2
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


# 3. Model fitting with proteomics data
# -------------------------------------
# The model was parametrized with kinetic data from literature where available.
# (see table 'models/salmonella/parameters.csv').
# In addition to that, we use MS derived protein abundances that were summed up for
# each functional sector of the cellular economy model.
# The kinetic parameters are randomly sampled and adjusted such that the error
# between model prediction and actual experimentally observed size of proteome sectors
# is minimized.
#
# Strategy: fit data to the Null mutant condition (no flagella expression)
# First import proteomics data:
df_mf = pd.read_csv("data/tables/sector_mass_fractions.tsv", delimiter="\t")
df_mf = (df_mf.groupby(["condition", "sector_short"])
    .agg("mean")
    .reset_index()
    .query("condition == 'EM16223'")
    .filter(["sector_short", "mass_fraction", "mean_growth_rate"])
)

# perform parameter sampling to find stable sets and
# improve solver performance (typical problem is over-constrainment)
n_iterations = 0
outdir = "results/salmonella/sampling/"
kcat = pd.Series([185, 102, 22, 122, 7.5, 5.5, 4e4], index=enz)
Km = pd.Series([2.5, 0.02, 0.08, 0.3, 1.50, 0.35, 1.0], index=enz)
hc = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index=enz)
a_fla = 0.002
best_error = np.inf

while n_iterations <= 50:
    n_iterations += 1
    iter = "{0:03d}".format(n_iterations)
    try:
        result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
        mass_predicted = result_ss.table.query("cex == 1.0").filter(["a_aab", "a_cbn", "a_etc", "a_fla", "a_lpb", "a_oth","a_rib", "a_tra"])
        mass_measured = df_mf.set_index("sector_short").loc[["Aab", "Cbn", "Etc", "Fla", "Lpb", "Oth", "Rib", "Tra"], "mass_fraction"]
        mu_pred = result_ss.table.query("cex == 1.0")["mu"].values[0]
        df_sampling = pd.DataFrame({"sector": mass_measured.index, "predicted": mass_predicted.values.flatten(), "measured": mass_measured.values.flatten()})
        df_sampling = pd.concat([df_sampling, pd.DataFrame({"sector": ["mu"], "predicted": mu_pred, "measured": [df_mf["mean_growth_rate"].values[0]]})], ignore_index=True)
        df_sampling["abs_error"] = abs(df_sampling["predicted"] - df_sampling["measured"])
        df_sampling["rel_error"] = list(map(lambda x: log2(x[0] / x[1]), zip(df_sampling["predicted"], df_sampling["measured"])))
        error = np.sum(abs(df_sampling["rel_error"]))
        if  (error < best_error and mu_pred < 2.0):
            best_error = error.copy()
            msg_error = f"""
            ---
            iteration {iter} with flagella {a_fla:.3f} and final
            growth rate {mu_pred:.2f}
            has error {error:.2f}
            """
            print(msg_error)
            # save current best parameter set and adjust kinetic parameters
            result_ss.table.to_csv(outdir + "steady_state_iter_" + iter + "_flag_" + str(a_fla) + ".csv")
            df_kinetic_params = pd.DataFrame({"kcat": kcat, "Km": Km, "hc": hc})
            df_kinetic_params.to_csv(outdir + "kinetic_params.csv")
            # write msg to log_file
            with open(outdir + "parameter_sampling.log", "a") as log_file:
                log_file.write(msg_error + "\n")
                log_file.write(str(df_sampling) + "\n")
                log_file.write(str(df_kinetic_params) + "\n")
            kcat = pd.Series(common.randomize(kcat, 0.85, 1.15), index=enz)
            Km = pd.Series(common.randomize(Km, 0.85, 1.15), index=enz)
        else:
            print(f"iteration {iter} with growth rate {mu_pred:.3f} has error {error:.2f}, not improving")
            raise ValueError("not improving")
    except:
        print("\n---\nmodel not solvable, trying next parameter set")
        kcat = pd.Series(common.randomize([185, 102, 22, 122, 7.5, 5.5, 4e4], 0.75, 1.25), index=enz)
        Km = pd.Series(common.randomize([2.5, 0.02, 0.08, 0.3, 1.50, 0.35, 1.0], 0.75, 1.25), index=enz)


# 4. run steady state model simulations
# -------------------------------------
#
# import desired parameter set
df_top_params = pd.read_csv("results/salmonella/sampling/top/kinetic_params_2026.csv", index_col=0)
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

# import results
df_climitation = []
for file in sorted(glob(outdir + "steady_state*.csv")):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("(flag_[0-9]+\\.[0-9]+)", file)[0]
    df = df.query("time != 0")
    df_climitation.append(df)

# combine into one df and plot
df_climitation = pd.concat(df_climitation, ignore_index=True)
common.plot_enzymes(df_climitation, outdir)
common.plot_properties(df_climitation, outdir)
common.plot_rates(df_climitation, outdir)


# 4.2 simulate substrate limitation with and without rotational ATP cost for flagella
# (set kcat of flagella to 0, but force protein cost)
outdir = "results/salmonella/rotation/"
for k, v in {"ATP": 4e4, "no_ATP": 0}.items():
    kcat["Fla"] = v
    retries = 0
    for a_fla in np.arange(0, 0.06, 0.01):
        while retries <= max_retries:
            try:
                result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
                result_ss.table.to_csv(outdir + "steady_state_flag_" + "{0:02.2f}".format(a_fla) + f"_{k}.csv")
                break
            except:
                print("\n---\nmodel not solvable, varying parameter")
                kcat["Fla"] = kcat["Fla"] + abs(np.random.normal(1) / 100)
                retries += 1
        retries = 0
# set kcat back to orig value
kcat["Fla"] = 4e4

# import results
df_rotation = []
for file in sorted(glob(outdir + "steady_state*.csv")):
    df = pd.read_csv(file)
    df["type"] = "steady_state"
    df["iteration"] = re.findall("(flag_[0-9]+\\.[0-9]+(\\_no)?_ATP)", file)[0][0]
    df = df.query("time != 0")
    df_rotation.append(df)

# combine into one df and plot
df_rotation = pd.concat(df_rotation, ignore_index=True)
common.plot_enzymes(df_rotation, outdir)
common.plot_properties(df_rotation, outdir)
common.plot_rates(df_rotation, outdir)


# 4.3 plot barchart of flagella with/without energy cost
df_cost = df_rotation[df_rotation["cex"] == 1.0]
df_cost["a_fla"] = round(df_cost["a_fla"] * 100)
df_cost["energy cost"] = df_cost["iteration"].str.contains("no_ATP")
df_cost["energy cost"] = df_cost["energy cost"].apply(
    lambda x: "without ATP cost (no rotation)" if x else "with ATP cost (rotation)"
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
ax.set(xlabel="% protein to flagella", ylabel="growth rate [h^-1]")
plt.legend(title="", fontsize="6", loc="lower right")
plt.grid(axis="both")

plt.subplot(1, 2, 2)
plt.title("relative growth penalty by rotation", loc="left", fontsize=10)
ax = sns.barplot(
    x=df_cost["a_fla"].astype(str),
    y=(max(df_cost["mu"]) - df_cost["mu"]) / max(df_cost["mu"]) * 100,
    hue=df_cost["energy cost"],
    palette=sns.color_palette("flare", n_colors=2),
)
ax.set(xlabel="% protein to flagella", ylabel="% reduction in growth")
plt.legend(title="", fontsize="6", loc="lower right")
plt.grid(axis="both")
plt.savefig(outdir + "energy_vs_protein_cost.png", dpi=182)
plt.savefig(outdir + "energy_vs_protein_cost.svg")


# 5. run dynamic model simulations
# -------------------------------------
#
# 5.1 simulate swimming at variable speed, depending on number of flagella
importlib.reload(dynamic)
for dist_init in [8500]:
    outdir = f"results/salmonella/swimming/{dist_init}/"
    os.makedirs(outdir, exist_ok=True)
    c_init = 5.0 # [mM] max substrate conc at gradient boundary
    time_init = 3 * 3600 # [sec] time for establishing substrate gradient
    time =  np.concatenate([[0, 0.1], np.arange(0.5, 10, 0.5)]) # [h]
    retries = 0
    for a_fla in np.concatenate([[0.005], np.arange(0.01, 0.06, 0.01)]):
        # pre-run steady state model to find good starting values for dynamic model (fix length and radius!)
        c_ex = round(common.diffusion_model(dist_init, time_init, 600, c_init), 3) # [mM]
        result_ss = steadystate.simulate(time, c_ex, c_ub, a_fla, kcat, Km, hc, remote)
        c_start = result_ss.c.apply(lambda x: round(x[1], 3))
        while retries <= max_retries:
            try:
                result_dy = dynamic.simulate(time, time_init, dist_init, c_init, c_ub, a_fla, kcat, Km, hc, c_start, remote)
                result_dy.table.to_csv(outdir + "dynamic_flag_" + "{0:02.3f}".format(a_fla) + ".csv")
                break
            except:
                print("\n---\nmodel not solvable, varying parameter")
                c_init = c_init + round(abs(np.random.normal(1)) / 100, 3)
                retries += 1
        retries = 0
    # 
    # import result tables
    df_dynamic = []
    for file in sorted(glob(outdir + "dynamic*.csv")):
        df = pd.read_csv(file)
        df["type"] = "steady_state"
        df["iteration"] = re.findall("(flag_[0-9]+\\.[0-9]+)", file)[0]
        df = df.query("time != 0")
        df_dynamic.append(df)
    # 
    # combine into one df and plot
    df_combined = pd.concat(df_dynamic, ignore_index=True)
    common.plot_enzymes(df_combined, outdir)
    common.plot_properties(df_combined, outdir)
    common.plot_rates(df_combined, outdir)
    #
    # export table with substrate gradient
    pd.DataFrame({
        "time_h": time_init / 3600,
        "substrate_initial_mM": round(c_init, 1),
        "dist_um": list(range(0, dist_init + 1, 100)),
        "substrate_mM": [common.diffusion_model(d, time_init, 600, round(c_init, 1)) for d in range(0, dist_init + 1, 100)]
    }).to_csv(outdir + "substrate_gradient.csv", index=False)
