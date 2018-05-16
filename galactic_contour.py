import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import (colors, cm, ticker)

from mpl_utils import mpl_style

plt.style.use(mpl_style)

stars = Table.read("dynamics_complete_test.csv")
lamosts = Table.read("lamost_gaia_result.fits")

xs = stars['GAL_LONG']
ys = stars['GAL_LAT']
xl = lamosts['GAL_LONG']
yl = lamosts['GAL_LAT']

grid_xs = np.linspace(0, 360, 160)
grid_xl = np.linspace(0, 360, 160)
grid_ys = np.linspace(-70, 90, 120)
grid_yl = np.linspace(-70, 90, 120)

grid_s, _, _ = np.histogram2d(xs, ys, bins=[grid_xs, grid_ys])
grid_l, _, _ = np.histogram2d(xl, yl, bins=[grid_xl, grid_yl])

div = np.divide(grid_s, grid_l, out=np.zeros_like(grid_s), where=grid_l!=0)
grid_sl = div
grid_sl = np.transpose(grid_sl)

grid_slm = np.ma.masked_where(grid_sl < 0.001, grid_sl)

fig, axes = plt.subplots(3, 1)
axes[0].plot(xs, ys,s=1, alpha=0.5, facecolor='#FF9309', edgecolor='none')
axes[0].set_xlabel(r'$l$ $[^\circ]$')
axes[0].set_ylabel(r'$b$ $[^\circ]$')

axes[1].plot(xl, yl,s=1, alpha=0.5, facecolor='#FF9309', edgecolor='none')
axes[1].set_xlabel(r'$l$ $[^\circ]$')
axes[1].set_ylabel(r'$b$ $[^\circ]$')

cmap = cm.viridis
cmap.set_bad(color='white')
axes[2].pcolormesh(grid_xs, grid_ys, grid_slm, cmap=cmap, edgecolors='None',
    norm=colors.LogNorm(vmin=0+0.001, vmax=0.1))
axes[2].colorbar().set_label(r'\textrm{Candidate fraction}', rotation=270)
axes[2].set_xlabel(r'$l$ $[^\circ]$')
axes[2].set_ylabel(r'$b$ $[^\circ]$')

for ax in axes:
    ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))

fig.tight_layout()

plt.show()

plt.savefig('DynamicFigures/'+'testcontoursub'+'.png')
