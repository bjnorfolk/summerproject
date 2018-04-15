import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.table import Table

data = Table.read("Bastar_candidates_with_abundances.csv")

ba = np.vstack([data["BA_FE_4554"], data["BA_FE_4934"]])
ba_err = np.vstack([data["BA_FE_4554_err"], data["BA_FE_4934_err"]])

sr = np.vstack([data["SR_FE_4077"], data["SR_FE_4215"]])
sr_err = np.vstack([data["SR_FE_4077_err"], data["SR_FE_4215_err"]])

sr_avg = np.nanmean(sr, axis=0)
ba_avg = np.nanmean(ba, axis=0)

ba_err[(ba_err == 0) + ~np.isfinite(ba_err)] = 0.20
sr_err[(sr_err == 0) + ~np.isfinite(sr_err)] = 0.20

ba_err = np.nanmean(ba_err, axis=0)
sr_err = np.nanmean(sr_err, axis=0)

hs_ls = ba_avg - sr_avg
hs_ls_err = np.sqrt(0.2**2 +  ba_err**2 + sr_err**2)

fig, ax = plt.subplots(figsize=(4, 4))

ax.scatter(data["cannon_m_h_1"], hs_ls,
           facecolor="#000000", edgecolor="#000000", linewidths=1, s=5, zorder=5,
           label="_")
ax.errorbar(data["cannon_m_h_1"], hs_ls,
            xerr=data["cannon_m_h_err_1"], yerr=hs_ls_err,
            fmt=None, ecolor="#BBBBBB", elinewidth=1, label="_")


cristallo = Table.read("cristallo_data_new.csv")
fehs = np.array(cristallo.dtype.names[1:]).astype(float)

masses = cristallo["Mass/[Fe/H]"]

colors = {
    "2.0": "k",
    "3.0": "g",
}

for i, mass in enumerate(masses):
    model_hs_ls = [cristallo[feh][i] for feh in cristallo.dtype.names[1:]]
    ax.plot(fehs, model_hs_ls, label=r"${{{:.1f}}}\,M_\odot$".format(mass),
            lw=5, alpha=0.75)

amanda = Table.read("amanda_yields.csv")
fehs = np.array(amanda.dtype.names[1:]).astype(float)

masses = amanda["mass/[Fe/H]"]

colors = {
    "2.0": "k",
    "3.0": "m",
}

for i, mass in enumerate(masses):
    model_hs_ls = [amanda[feh][i] for feh in amanda.dtype.names[1:]]
    ax.plot(fehs, model_hs_ls, label=r"${{{:.1f}}}\,M_\odot$".format(mass),
            lw=5, alpha=0.75)


plt.legend()
ax.set_xlim(-2.3, 0.4)
ax.set_ylim(-1, 1.5)

ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(MaxNLocator(6))

ax.legend(ncol=1, loc="upper left", frameon=True, fontsize=10)
ax.set_xlabel(r"Metallicity [Fe/H]")
ax.set_ylabel(r"Surface [Ba/Sr]")

fig.tight_layout()

fig.savefig("yields_test.pdf", dpi=300, overwrite=True)
