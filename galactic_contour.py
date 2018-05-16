import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

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

cmap = plt.cm.viridis
cmap.set_bad(color='white')
plt.pcolormesh(grid_xs, grid_ys, grid_slm, cmap=cmap, edgecolors = 'None', vmin =0 , vmax = 0.1)
plt.colorbar()

plt.title('Contour Percentage')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Lattitude')
plt.savefig('DynamicFigures/'+'testcontour'+'.png')



