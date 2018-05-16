import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib.colors as colors


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

fig = plt.figure()  # create a figure object
ax1 = fig.add_subplot(3, 1, 1)  # create an axes object in the figure
ax1.fig.plot(xs, ys,s=1, alpha=0.5, facecolor='#FF9309', edgecolor='none')
ax1.fig.set_xlabel('Galactic Longitude')
ax1.fig.set_ylabel('Galactic Lattitude')

ax2 = fig.add_subplot(3, 1, 2)
ax2.fig.plot(xl, yl,s=1, alpha=0.5, facecolor='#FF9309', edgecolor='none')
ax2.fig.set_xlabel('Galactic Longitude')
ax2.fig.set_ylabel('Galactic Lattitude')

ax3 = fig.add_subplot(3, 1, 3)
cmap = ax3.cm.viridis
cmap.set_bad(color='white')
ax3.fig.pcolormesh(grid_xs, grid_ys, grid_slm, cmap=cmap, edgecolors = 'None', norm=colors.LogNorm(vmin=0+0.001, vmax=0.1))
ax3.fig.colorbar().set_label('Candidate fraction', rotation=270)
ax3.fig.set_xlabel('Galactic Longitude')
ax3.fig.set_ylabel('Galactic Lattitude')


plt.savefig('DynamicFigures/'+'testcontoursub'+'.png')
