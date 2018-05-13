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
ecc = np.zeros((L, N))
ecc = ecc.flatten(L)
apocenter = np.zeros((L, N))
apocenter = apocenter.flatten(L)

orbit_energy = np.zeros((L, N))
orbit_energy = orbit_energy.flatten(L)

angular_mom_z = np.zeros((L, N))
angular_mom_z = angular_mom_z.flatten(L)

L = len(controls)

U_lamost = np.zeros((L, N))
U_lamost = U.flatten(L)
V_lamost = np.zeros((L, N))
V_lamost = V.flatten(L)
W_lamost = np.zeros((L, N))
W_lamost = W.flatten(L)
ecc_lamost = np.zeros((L, N))
ecc_lamost = ecc.flatten(L)
apocenter_lamost = np.zeros((L, N))
apocenter_lamost = apocenter_lamost.flatten(L)

orbit_energy_lamost = np.zeros((L, N))
orbit_energy_lamost = orbit_energy_lamost.flatten(L)

angular_mom_z_lamost = np.zeros((L, N))
angular_mom_z_lamost = angular_mom_z_lamost.flatten(L)

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

  #fig = orbit.plot()
  #plt.show()
  if (star["Ba_only_candidate"]=='true' and star["Ba_Sr_candidate"]=='false'):
   #vel = orbit.vel_components.keys() 
   #v_x = orbit.v_x
   #v_y = orbit.v_y
   #v_z = orbit.v_z
   Ba_velx.append(w0.v_x.to(u.km/u.s).value)
   Ba_vely.append(w0.v_y.to(u.km/u.s).value)
   Ba_velz.append(w0.v_z.to(u.km/u.s).value)
   Ba_y.append(np.sqrt(Ba_velx[-1]**2+Ba_velz[-1]**2))
   
  if (star["Sr_only_candidate"]=='true' and star["Ba_Sr_candidate"]=='false'):
   #vel = orbit.vel_components.keys() 
   #v_x = orbit.v_x
   #v_y = orbit.v_y
   #v_z = orbit.v_z
   Sr_velx.append(w0.v_x.to(u.km/u.s).value)
   Sr_vely.append(w0.v_y.to(u.km/u.s).value)
   Sr_velz.append(w0.v_z.to(u.km/u.s).value)
   Sr_y.append(np.sqrt(Sr_velx[-1]**2+Sr_velz[-1]**2))
   
  if (star["Ba_Sr_candidate"]=='true'):
   #vel = orbit.vel_components.keys() 
   #v_x = orbit.v_x
   #v_y = orbit.v_y
   #v_z = orbit.v_z
   Barium_velx.append(w0.v_x.to(u.km/u.s).value)
   Barium_vely.append(w0.v_y.to(u.km/u.s).value)
   Barium_velz.append(w0.v_z.to(u.km/u.s).value)
   Barium_y.append(np.sqrt(Barium_velx[-1]**2+Barium_velz[-1]**2))
  U[index]= w0.v_x.to(u.km/u.s).value
  V[index]= w0.v_y.to(u.km/u.s).value
  W[index]= w0.v_z.to(u.km/u.s).value
  apocenter[index]= orbit.apocenter().value
  ecc[index]= orbit.eccentricity()
  orbit_energy[index] = np.mean(orbit.energy().to(u.km * u.km / u.s / u.s)).value
  angular_mom_z[index] = np.mean(orbit.angular_momentum()[2].to(u.kpc * u.km / u.s)).value

stars["apocenter"] = apocenter
stars["ecc"] = ecc
stars["orbit_energy"] = orbit_energy
stars["angular_mom_z"] = angular_mom_z
stars["U_vel(km/s)"] = U
stars["V_vel(km/s)"] = V
stars["W_vel(km/s)"] = W
stars.write("velocity_dynamics.fits", overwrite=True)
stars.write("velocity_dynamics.csv", overwrite=True) 
  #odict_keys(['v_x', 'v_y', 'v_z'])
 #plt.plot(v_z, sqrt(v_x**2+v_y**2))
#plt.plot(Ba_vely, Ba_y, marker='D')
'''
for index, control in enumerate(controls):

  potential = gp.MilkyWayPotential()

  icrs = coord.ICRS(ra=control["ra"] * u.deg,
                  dec=control["dec"] * u.deg,
                  distance=1.0 / control["parallax"] * u.kpc,
                  pm_ra_cosdec=control["pmra"]*u.mas/u.yr,
                  pm_dec=control["pmdec"]*u.mas/u.yr,
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
  
  C_velx.append(w0.v_x.to(u.km/u.s).value)
  C_vely.append(w0.v_y.to(u.km/u.s).value)
  C_velz.append(w0.v_z.to(u.km/u.s).value)
  C_y.append(np.sqrt(C_velx[-1]**2+C_velz[-1]**2))
  U_lamost[index]= w0.v_x.to(u.km/u.s).value
  V_lamost[index]= w0.v_y.to(u.km/u.s).value
  W_lamost[index]= w0.v_z.to(u.km/u.s).value
  apocenter_lamost[index]= orbit.apocenter().value
  ecc_lamost[index]= orbit.eccentricity()
  orbit_energy_lamost[index] = np.mean(orbit.energy().to(u.km * u.km / u.s / u.s)).value
  angular_mom_z_lamost[index] = np.mean(orbit.angular_momentum()[2].to(u.kpc * u.km / u.s)).value

controls["apocenter"] = apocenter_lamost
controls["ecc"] = ecc_lamost
controls["orbit_energy"] = orbit_energy_lamost
controls["angular_mom_z"] = angular_mom_z_lamost
controls["U_vel"] = U_lamost
controls["V_vel"] = V_lamost
controls["W_vel"] = W_lamost   
'''

thinx = np.linspace(-85,85,10000) 
thiny = np.sqrt(85**2-thinx**2)
thickx = np.linspace(-180,180,100000) 
thicky = np.sqrt(180**2-thickx**2)

fig, ax = plt.subplots()
#ax.scatter(C_vely, C_y,s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none')
ax.scatter(thinx, thiny,s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none')
ax.scatter(thickx, thicky,s=1, alpha=0.5, facecolor='#1B0000', edgecolor='none')
ax.scatter(Ba_vely, Ba_y,label='Ba',s=5, alpha=0.5, facecolor='#B50B0E', edgecolor='none',zorder=2)
ax.scatter(Sr_vely, Sr_y, label='Sr',s=5, alpha=0.5, facecolor='#030CED', edgecolor='none',zorder=2)
ax.scatter(Barium_vely, Barium_y, label='BaSr',s=5, alpha=0.5, facecolor='#ED28BF', edgecolor='none',zorder=2)

plt.legend(loc='upper left')
plt.text(-25, 0, 'Thin', fontsize=12)
plt.text(-90,75, 'Thick', fontsize=12)
plt.text(-150, 150, 'Halo', fontsize=12)
plt.xlim(310, -310)
plt.ylim(0, 400)
plt.title('Toomre Diagram')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
plt.savefig('testtoomre'+'.png')

#`fig, ax = plt.subplots()`

#`ax.scatter`
#because `plt.scatter` just plots it on the most recently created figure, which may or may not be the figure you wanted it to plot on (e.g., if you are making multiple figures at one time)
"""
plt.plot(orbit.t, orbit.spherical.distance, marker='None')

per, per_times = orbit.pericenter(return_times=True)
apo, apo_times = orbit.apocenter(return_times=True)

for t in per_times:
    plt.axvline(t.value, color='#67a9cf')
()
for t in apo_times:
    plt.axvline(t.value, color='#ef8a62')

plt.xlabel('$t$ [{0}]'.format(orbit.t.unit.to_string('latex')))
plt.ylabel('$r$ [{0}]'.format(orbit.x.unit.to_string('latex')))


n_samples = 968
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


