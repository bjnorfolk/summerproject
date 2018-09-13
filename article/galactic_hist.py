import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import (colors, cm, ticker)

from mpl_utils import mpl_style

plt.style.use(mpl_style)

stars = Table.read("dynamics_complete_test.csv")
lamosts = Table.read("lamost_gaia_result.fits")


s = stars['GAL_LAT']
l = lamosts['GAL_LAT']

b = np.absolute(s)
b_l = np.absolute(l)
s_h, bin_edges = np.histogram(b, 89)
l_h, bin_edges = np.histogram(b_l, 89)
div = np.divide(s_h, l_h)
menStd = np.sqrt(s_h)
menStd_f = ((s_h + np.sqrt(s_h))/(l_h + np.sqrt(l_h))) - (s_h/l_h)

fig, axes = plt.subplots(2, 1, figsize=(8, 12))
axes[0].bar(bin_edges[:-1], s_h, width = 1, facecolor='#FF9309', alpha=0.5,
	                               yerr=menStd)
axes[0].set_xlim(min(bin_edges), max(bin_edges))
axes[0].set_xlabel(r'$|b|$ $[^\circ]$')
axes[0].set_ylabel(r'$N$')

axes[1].bar(bin_edges[:-1], div, width = 1, facecolor='#FF9309', alpha=0.5,
	                               yerr=menStd_f)
axes[1].set_xlim(min(bin_edges), max(bin_edges))
axes[1].set_xlabel(r'$|b|$ $[^\circ]$')
axes[1].set_ylabel(r'$N/N$ $_{LAMOST}$')

for ax in axes:
    ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))

fig.tight_layout()
plt.show()

fig.savefig('DynamicFigures/'+'histogram'+'.png')
fig.savefig("DynamicFigures/histogram.pdf", dpi=150)
plt.close()