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

#d = Table.read("Updated_88_data_gaia.fits")
#d.write("Updated_88_data_gaia.csv")

#data=pd.read_csv('Updated_88_data_gaia.csv')

#csv=np.array(data)

stars = Table.read("Updated_88_data_gaia.fits")


for index, star in enumerate(stars):

  potential = gp.MilkyWayPotential()

  icrs = coord.ICRS(ra=star["ra"] * u.deg,
                  dec=star["dec"] * u.deg,
                  distance=1.0 / star["parallax"] * u.kpc,
                  pm_ra_cosdec=star["pmra"]*u.mas/u.yr,
                  pm_dec=star["pmdec"]*u.mas/u.yr,
                  radial_velocity=-0*u.km/u.s)

  """
  icrs_err = coord.ICRS(ra=0*u.deg, dec=0*u.deg, distance=6*u.kpc,
                      pm_ra_cosdec=0.009*u.mas/u.yr,
                      pm_dec=0.009*u.mas/u.yr,
                      radial_velocity=0.1*u.km/u.s)

  """
  v_sun = coord.CartesianDifferential([11.1, 250, 7.25]*u.km/u.s)
  gc_frame = coord.Galactocentric(galcen_distance=8.3*u.kpc,
                                z_sun=0*u.pc,
                                galcen_v_sun=v_sun)


  gc = icrs.transform_to(gc_frame)

  w0 = gd.PhaseSpacePosition(gc.data)
  orbit = potential.integrate_orbit(w0, dt=-0.5*u.Myr, n_steps=2000)

  fig = orbit.plot()
  plt.show()
  raise a

"""
plt.plot(orbit.t, orbit.spherical.distance, marker='None')

per, per_times = orbit.pericenter(return_times=True)
apo, apo_times = orbit.apocenter(return_times=True)

for t in per_times:
    plt.axvline(t.value, color='#67a9cf')

for t in apo_times:
    plt.axvline(t.value, color='#ef8a62')

plt.xlabel('$t$ [{0}]'.format(orbit.t.unit.to_string('latex')))
plt.ylabel('$r$ [{0}]'.format(orbit.x.unit.to_string('latex')))


n_samples = 256
dist = np.random.normal(icrs.distance.value, icrs_err.distance.value,
                        n_samples) * icrs.distance.unit

pm_ra_cosdec = np.random.normal(icrs.pm_ra_cosdec.value,
                                icrs_err.pm_ra_cosdec.value,
                                n_samples) * icrs.pm_ra_cosdec.unit

pm_dec = np.random.normal(icrs.pm_dec.value,
                          icrs_err.pm_dec.value,
                          n_samples) * icrs.pm_dec.unit

rv = np.random.normal(icrs.radial_velocity.value, icrs_err.radial_velocity.value,
                      n_samples) * icrs.radial_velocity.unit

ra = np.full(n_samples, icrs.ra.degree) * u.degree
dec = np.full(n_samples, icrs.dec.degree) * u.degree


icrs_samples = coord.ICRS(ra=ra, dec=dec, distance=dist,
                          pm_ra_cosdec=pm_ra_cosdec,
                          pm_dec=pm_dec, radial_velocity=rv)


icrs_samples.shape


gc_samples = icrs_samples.transform_to(gc_frame)


w0_samples = gd.PhaseSpacePosition(gc_samples.data)


orbit_samples = potential.integrate_orbit(w0_samples, dt=-0.5*u.Myr, n_steps=10000)
orbit_samples.shape


pers = []
apos = []
eccs = []
for n in range(n_samples):
    orbit = orbit_samples[:,n]
    pers.append(orbit.pericenter())
    apos.append(orbit.apocenter())
    eccs.append(orbit.eccentricity())

pers = u.Quantity(pers)
apos = u.Quantity(apos)
eccs = u.Quantity(eccs)


fig, axes = plt.subplots(1, 3, figsize=(12, 4))

axes[0].hist(pers, bins='auto')
axes[0].set_xlabel('pericenter [kpc]')

axes[1].hist(apos, bins='auto')
axes[1].set_xlabel('apocenter [kpc]')

axes[2].hist(eccs, bins='auto')
axes[2].set_xlabel('eccentricity');
"""


