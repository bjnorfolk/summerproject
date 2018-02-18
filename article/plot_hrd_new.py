import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.table import Table

data = Table.read("Bastar_candidates.csv")

fig, ax = plt.subplots(figsize=(4, 4))

ax.scatter(data["cannon_teff_1"], data["cannon_logg_1"],
           facecolor="#000000", edgecolor="#000000", linewidths=1, s=5, zorder=100)
ax.errorbar(data["cannon_teff_1"], data["cannon_logg_1"],
            xerr=data["cannon_teff_err_1"], yerr=data["cannon_logg_err_1"],
            fmt=None, ecolor="#666666", elinewidth=1)

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(4, 0)
ax.set_yticks([0, 1, 2, 3, 4])

lhs = 5400


ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(MaxNLocator(6))


ax.set_xlabel(r"Effective temperature $T_{\rm eff}$ $({\rm K})$")
ax.set_ylabel(r"Surface gravity $\log{g}$")

fig.tight_layout()
fig.savefig("hrd_new.pdf", dpi=300, overwrite=True)
