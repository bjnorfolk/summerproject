import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from mpl_utils import mpl_style
from matplotlib import (colors, cm, ticker)
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
import gala.dynamics as gd
import gala.integrate as gi
import gala.potential as gp
from gala.units import galactic

stars = Table.read("velocity_dynamics.csv")
lamost = Table.read("lamostUVW_dynamics.csv")
Ba_velx = []
Ba_vely = []
Ba_velz = []
Ba_y = []
Sr_velx = []
Sr_vely = []
Sr_velz = []
Sr_y = []
Barium_velx = []
Barium_vely = []
Barium_velz = []
Barium_y = []

for index, star in enumerate(stars):

  icrs = coord.ICRS(ra=star["ra_1"] * u.deg,
                  dec=star["dec_1"] * u.deg,
                  distance=1.0 / star["parallax_1"] * u.kpc,
                  pm_ra_cosdec=star["pmra_1"]*u.mas/u.yr,
                  pm_dec=star["pmdec_1"]*u.mas/u.yr,
                  radial_velocity=star['radial_velocity_1']*u.km/u.s)

  """
  icrs_err = coord.ICRS(ra=0*u.deg, dec=0*u.deg, distance=6*u.kpc,
                      pm_ra_cosdec=0.009*u.mas/u.yr,
                      pm_dec=0.009*u.mas/u.yr,
                      radial_velocity=0.1*u.km/u.s)

  """

  v_sun = coord.CartesianDifferential([0, 0, 0]*u.km/u.s)
  gc_frame = coord.Galactocentric(galcen_distance=8.3*u.kpc,
                                z_sun=0*u.pc,
                                galcen_v_sun=v_sun)


  gc = icrs.transform_to(gc_frame)

  w0 = gd.PhaseSpacePosition(gc.data)

  if (star["Ba_only_candidate_1"]=='true' and star["Ba_Sr_candidate_1"]=='false'):

   Ba_velx.append(w0.v_x.to(u.km/u.s).value)
   Ba_vely.append(w0.v_y.to(u.km/u.s).value)
   Ba_velz.append(w0.v_z.to(u.km/u.s).value)
   Ba_y.append(np.sqrt(Ba_velx[-1]**2+Ba_velz[-1]**2))

   
  if (star["Sr_only_candidate_1"]=='true' and star["Ba_Sr_candidate_1"]=='false'):

   Sr_velx.append(w0.v_x.to(u.km/u.s).value)
   Sr_vely.append(w0.v_y.to(u.km/u.s).value)
   Sr_velz.append(w0.v_z.to(u.km/u.s).value)
   Sr_y.append(np.sqrt(Sr_velx[-1]**2+Sr_velz[-1]**2))

   
  if (star["Ba_Sr_candidate_1"]=='true'):

   Barium_velx.append(w0.v_x.to(u.km/u.s).value)
   Barium_vely.append(w0.v_y.to(u.km/u.s).value)
   Barium_velz.append(w0.v_z.to(u.km/u.s).value)
   Barium_y.append(np.sqrt(Barium_velx[-1]**2+Barium_velz[-1]**2))


xs = stars['V_vel(km/s)_1']
ys = np.sqrt(stars['U_vel(km/s)_1']**2+stars['W_vel(km/s)_1']**2)
xl = lamost['V_vel']
yl = np.sqrt(lamost['U_vel']**2+lamost['W_vel']**2)


grid_xs = np.linspace(-310, 310, 620)
grid_xl = np.linspace(-310, 310, 620)
grid_ys = np.linspace(0, 400, 400)
grid_yl = np.linspace(0, 400, 400)

grid_s, _, _ = np.histogram2d(xs, ys, bins=[grid_xs, grid_ys])
grid_l, _, _ = np.histogram2d(xl, yl, bins=[grid_xl, grid_yl])
div = np.divide(grid_s, grid_l, out=np.zeros_like(grid_s), where=grid_l!=0)
grid_sl = div
grid_sl = np.transpose(grid_sl)

grid_slm = np.ma.masked_where(grid_sl < 0.001, grid_sl)

thinx = np.linspace(-85,85,10000) 
thiny = np.sqrt(85**2-thinx**2)
thickx = np.linspace(-180,180,100000) 
thicky = np.sqrt(180**2-thickx**2)

fig, ax = plt.subplots()
cmap = cm.Greys
cmap.set_bad(color='white')
pmesh = ax.pcolormesh(grid_xs, grid_ys, grid_slm, 
    cmap=cmap, edgecolors='None', norm=colors.LogNorm(vmin=0+0.001, vmax=0.1))
cbar = plt.colorbar(pmesh)
cbar.set_label(r'\textrm{Candidate fraction}', rotation=270)

ax.scatter(Ba_vely, Ba_y,label='Ba', s=1, alpha=0.5, facecolor='#030CED', edgecolor='none',
                rasterized=True)
ax.scatter(Sr_vely, Sr_y, label='Sr', s=1, alpha=0.5, facecolor='#ED28BF', edgecolor='none',
                rasterized=True)
ax.scatter(Barium_vely, Barium_y, label='BaSr', s=1, alpha=0.5, facecolor='#B50B0E', edgecolor='none',
                rasterized=True)
ax.scatter(thinx, thiny, s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none',
                rasterized=True)
ax.scatter(thickx, thicky, s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none',
                rasterized=True)
ax.set_xlabel('V (km/s)')
ax.set_ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')

plt.xlim(-310, 310)
plt.ylim(0, 400)
plt.legend(loc='upper left')
plt.text(30, 10, 'Thin', fontsize=12)
plt.text(90,75, 'Thick', fontsize=12)
plt.text(150, 150, 'Halo', fontsize=12)
plt.tight_layout()

plt.savefig('DynamicFigures/'+'toomre'+'.png')
