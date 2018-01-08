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

potential = gp.MilkyWayPotential()

icrs = coord.ICRS(ra=337.9 * u.deg, dec=2.16 * u.deg, distance=178* u.kpc, pm_ra_cosdec=51.69*u.mas/u.yr, pm_dec= -328.64*u.mas/u.yr, radial_velocity= -210.4*299792.458*u.km/u.s)

v_sun = coord.CartesianDifferential([11.1, 250, 7.25]*u.km/u.s)
gc_frame = coord.Galactocentric(galcen_distance=8.3*u.kpc, z_sun=0*u.pc, galcen_v_sun=v_sun)


gc = icrs.transform_to(gc_frame)

w0 = gd.PhaseSpacePosition(gc.data)
orbit = potential.integrate_orbit(w0, dt=-0.5*u.Myr, n_steps=2000)


#print("orbit_energy: {}".format(np.mean(orbit.energy().to(u.km * u.km / u.s / u.s)).value)
#print("angular_mom_z: {}".format(np.mean(orbit.angular_momentum()[2].to(u.kpc * u.km / u.s)).value)
