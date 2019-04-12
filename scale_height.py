# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

# Gala
from gala.mpl_style import mpl_style
plt.style.use(mpl_style)
import gala.dynamics as gd
import gala.integrate as gi
import gala.potential as gp
from gala.units import galactic

stars = Table.read("velocity_dynamics.csv")
L = len(stars)
N=1

x = np.zeros((L, N))
x = x.flatten(L)
y = np.zeros((L, N))
y = y.flatten(L)
z = np.zeros((L, N))
z = z.flatten(L)

for index, star in enumerate(stars):
  
  icrs = coord.ICRS(ra=star["ra_1_1"] * u.deg,
                  dec=star["dec_1_1"] * u.deg,
                  distance=1.0 / star["parallax_1"] * u.kpc,
                  pm_ra_cosdec=star["pmra_1"]*u.mas/u.yr,
                  pm_dec=star["pmdec_1"]*u.mas/u.yr,
                  radial_velocity=star['radial_velocity_1']*u.km/u.s)

  v_sun = coord.CartesianDifferential([0, 0, 0]*u.km/u.s)
  gc_frame = coord.Galactocentric(galcen_distance=8.3*u.kpc,
                                z_sun=0*u.pc,
                                galcen_v_sun=v_sun)


  gc = icrs.transform_to(gc_frame)

  w0 = gd.PhaseSpacePosition(gc.data)

  x[index]= w0.x.value
  y[index]= w0.y.value
  z[index]= w0.z.value 

stars["x(kpc)"] = x
stars["y(kpc)"] = y
stars["z(kpc)"] = z
stars.write("dynamics_finalv2.fits", overwrite=True)
stars.write("dynamics_finalv2.csv", overwrite=True) 


