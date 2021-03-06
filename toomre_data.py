# Some imports we'll need later:

# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Gala
from gala.mpl_style import mpl_style
plt.style.use(mpl_style)
import gala.dynamics as gd
import gala.integrate as gi
import gala.potential as gp
from gala.units import galactic

stars = Table.read("dynamics_complete.csv")

controls = Table.read("lamost_gaia_result.fits")
L = len(stars)
N=1

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
C_velx = []
C_vely = []
C_velz = []
C_y = []


U = np.zeros((L, N))
U = U.flatten(L)
V = np.zeros((L, N))
V = V.flatten(L)
W = np.zeros((L, N))
W = W.flatten(L)

L = len(controls)

U_lamost = np.zeros((L, N))
U_lamost = U_lamost.flatten(L)
V_lamost = np.zeros((L, N))
V_lamost = V_lamost.flatten(L)
W_lamost = np.zeros((L, N))
W_lamost = W_lamost.flatten(L)

for index, star in enumerate(stars):

  icrs = coord.ICRS(ra=star["ra"] * u.deg,
                  dec=star["dec"] * u.deg,
                  distance=1.0 / star["parallax"] * u.kpc,
                  pm_ra_cosdec=star["pmra"]*u.mas/u.yr,
                  pm_dec=star["pmdec"]*u.mas/u.yr,
                  radial_velocity=star['radial_velocity']*u.km/u.s)

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

  if (star["Ba_only_candidate"]=='true' and star["Ba_Sr_candidate"]=='false'):

   Ba_velx.append(w0.v_x.to(u.km/u.s).value)
   Ba_vely.append(w0.v_y.to(u.km/u.s).value)
   Ba_velz.append(w0.v_z.to(u.km/u.s).value)
   Ba_y.append(np.sqrt(Ba_velx[-1]**2+Ba_velz[-1]**2))

   
  if (star["Sr_only_candidate"]=='true' and star["Ba_Sr_candidate"]=='false'):

   Sr_velx.append(w0.v_x.to(u.km/u.s).value)
   Sr_vely.append(w0.v_y.to(u.km/u.s).value)
   Sr_velz.append(w0.v_z.to(u.km/u.s).value)
   Sr_y.append(np.sqrt(Sr_velx[-1]**2+Sr_velz[-1]**2))

   
  if (star["Ba_Sr_candidate"]=='true'):

   Barium_velx.append(w0.v_x.to(u.km/u.s).value)
   Barium_vely.append(w0.v_y.to(u.km/u.s).value)
   Barium_velz.append(w0.v_z.to(u.km/u.s).value)
   Barium_y.append(np.sqrt(Barium_velx[-1]**2+Barium_velz[-1]**2))

  U[index]= w0.v_x.to(u.km/u.s).value
  V[index]= w0.v_y.to(u.km/u.s).value
  W[index]= w0.v_z.to(u.km/u.s).value


stars["U_vel(km/s)"] = U
stars["V_vel(km/s)"] = V
stars["W_vel(km/s)"] = W
stars.write("velocity_dynamics.fits", overwrite=True)
stars.write("velocity_dynamics.csv", overwrite=True) 


for index, control in enumerate(controls):

  icrs = coord.ICRS(ra=control["ra"] * u.deg,
                  dec=control["dec"] * u.deg,
                  distance=1.0 / control["parallax"] * u.kpc,
                  pm_ra_cosdec=control["pmra"]*u.mas/u.yr,
                  pm_dec=control["pmdec"]*u.mas/u.yr,
                  radial_velocity=control['radial_velocity']*u.km/u.s)

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

  C_velx.append(w0.v_x.to(u.km/u.s).value)
  C_vely.append(w0.v_y.to(u.km/u.s).value)
  C_velz.append(w0.v_z.to(u.km/u.s).value)
  C_y.append(np.sqrt(C_velx[-1]**2+C_velz[-1]**2))

  U_lamost[index]= w0.v_x.to(u.km/u.s).value
  V_lamost[index]= w0.v_y.to(u.km/u.s).value
  W_lamost[index]= w0.v_z.to(u.km/u.s).value

controls["U_vel"] = U_lamost
controls["V_vel"] = V_lamost
controls["W_vel"] = W_lamost 
controls.write("lamostUVW_dynamics.fits", overwrite=True)
controls.write("lamostUVW_dynamics.csv", overwrite=True)   


thinx = np.linspace(-85,85,10000) 
thiny = np.sqrt(85**2-thinx**2)
thickx = np.linspace(-180,180,100000) 
thicky = np.sqrt(180**2-thickx**2)

fig, ax = plt.subplots()
ax.scatter(C_vely, C_y,s=1, alpha=0.5, facecolor='#FF9309', edgecolor='none')
ax.scatter(thinx, thiny,s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none')
ax.scatter(thickx, thicky,s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none')
ax.scatter(Ba_vely, Ba_y,label='Ba',s=5, alpha=0.5, facecolor='#B50B0E', edgecolor='none',zorder=2)
ax.scatter(Sr_vely, Sr_y, label='Sr',s=5, alpha=0.5, facecolor='#030CED', edgecolor='none',zorder=2)
ax.scatter(Barium_vely, Barium_y, label='BaSr',s=5, alpha=0.5, facecolor='#ED28BF', edgecolor='none',zorder=2)

plt.legend(loc='upper left')
plt.text(30, 10, 'Thin', fontsize=12)
plt.text(90,75, 'Thick', fontsize=12)
plt.text(150, 150, 'Halo', fontsize=12)
plt.xlim(-310, 310)
plt.ylim(0, 400)
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
plt.savefig('DynamicFigures/'+'toomre'+'.png')


