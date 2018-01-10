import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table

import lamost
import utils
import pandas as pd

plt.ioff()
catalog = lamost.load_catalog()
wavelengths = lamost.common_wavelengths
N, P = (len(catalog), wavelengths.size)

# Open the data arrays
all_observed_flux = np.memmap(
    os.path.join(lamost.LAMOST_PATH, "observed_flux.memmap"),
    mode="r", dtype='float32', shape=(N, P))

all_observed_ivar = np.memmap(
    os.path.join(lamost.LAMOST_PATH, "observed_ivar.memmap"),
    mode="r", dtype='float32', shape=(N, P))

all_model_flux = np.memmap(
    os.path.join(lamost.LAMOST_PATH, "model_flux.memmap"),
    mode="r", dtype="float32", shape=(N, P))

d = Table.read("Tc_catalog.fits")
d.write("Tc_catalog.csv",overwrite=True)

data=pd.read_csv('Tc_catalog.csv')

csv=np.array(data)
#Barium
Ba4554_wavelength = csv[:,29] #
Ba4554_amp = csv[:,13]
Ba4554_amperr = csv[:,14]
#Strontium
St4077_wavelength = csv[:,31] #
St4077_amp = csv[:,15]
St4077_amperr = csv[:,16]

#Technetium
Tc4049_wavelength = csv[:,33] #
Tc4049_amp = csv[:,17]
Tc4049_amperr = csv[:,18]

Tc4238_wavelength = csv[:,35] #
Tc4238_amp = csv[:,19]
Tc4238_amperr = csv[:,20]

Tc4262_wavelength = csv[:,37] #
Tc4262_amp = csv[:,21]
Tc4262_amperr = csv[:,22]

Tc4297_wavelength = csv[:,39] #
Tc4297_amp = csv[:,23]
Tc4297_amperr = csv[:,24]

Tc5924_wavelength = csv[:,41] #
Tc5924_amp = csv[:,25]
Tc5924_amperr = csv[:,26]


#Sodium
Na5889_wavelength = csv[:,43] #
Na5889_amp = csv[:,27]
Na5889_amperr = csv[:,28]





snrg = csv[:,12]
chi = csv[:,9]
starID=csv[:,0]

teff= csv[:,1]
surfg= csv[:,2]
met= csv[:,3]

#Barium Numbers
N_Ba4554 = (np.abs(Ba4554_amp/(Ba4554_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba4554_amp < -0.05) \
       * (np.abs(Ba4554_wavelength - 4554) < 2)\
       * (snrg > 30)

#Stontium Lines
N_St4077 = (np.abs(St4077_amp/(St4077_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (St4077_amp < -0.05) \
       * (np.abs(St4077_wavelength - 4077) < 2) \
       * (St4077_amperr.astype(float) > 0)\
       * (snrg > 30)


#Technetium Lines
N_Tc4049 = (np.abs(Tc4049_amp/(Tc4049_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Tc4049_amp < -0.05) \
       * (np.abs(Tc4049_wavelength - 4049) < 2) \
       * (Tc4049_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Tc4238 = (np.abs(Tc4238_amp/(Tc4238_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Tc4238_amp < -0.05) \
       * (np.abs(Tc4238_wavelength - 4238) < 2) \
       * (Tc4238_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Tc4262 = (np.abs(Tc4262_amp/(Tc4262_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Tc4262_amp < -0.05) \
       * (np.abs(Tc4262_wavelength - 4262) < 2) \
       * (Tc4262_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Tc4297 = (np.abs(Tc4297_amp/(Tc4297_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Tc4297_amp < -0.05) \
       * (np.abs(Tc4297_wavelength - 4297) < 2) \
       * (Tc4297_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Tc5924 = (np.abs(Tc5924_amp/(Tc5924_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Tc5924_amp < -0.05) \
       * (np.abs(Tc5924_wavelength - 5924) < 2) \
       * (Tc5924_amperr.astype(float) > 0)\
       * (snrg > 30)

#Sodium
N_Na5889 = (np.abs(Na5889_amp/(Na5889_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Na5889_amp < -0.05) \
       * (np.abs(Na5889_wavelength - 5889) < 2) \
       * (Na5889_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Bariumstars = N_St4077*N_Ba4554
N_Tc = N_St4077*N_Ba4554*N_Tc4049*N_Tc4238*N_Tc4262*N_Tc4297*N_Tc5924
N_TcAll = N_Tc4049*N_Tc4238*N_Tc4262*N_Tc4297*N_Tc5924
N_Na = N_St4077*N_Ba4554*N_Na5889

print("Barium stars: {}".format(sum(N_Bariumstars)))
print("Sodium 5889: {}".format(sum(N_Na5889)))
print("Technetium 4049: {}".format(sum(N_Tc4049)))
print("Technetium 4238: {}".format(sum(N_Tc4238)))
print("Technetium 4262: {}".format(sum(N_Tc4262)))
print("Technetium 4297: {}".format(sum(N_Tc4297)))
print("Technetium 5924: {}".format(sum(N_Tc5924)))
print("Technetium all: {}".format(sum(N_TcAll)))
print("Barium star Tc Matches: {}".format(sum(N_Tc)))
print("Barium star Na Matches: {}".format(sum(N_Na)))


"""

for index, star in enumerate(csv):
    if not N_Tc[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

    fig.axes[0].set_xlim(4000, 4320)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(4049, c="#666666", zorder=-1)
    fig.axes[1].axvline(4238, c="#666666", zorder=-1)
    fig.axes[1].axvline(4262, c="#666666", zorder=-1)
    fig.axes[1].axvline(4297, c="#666666", zorder=-1)
    string="Barium star Tc"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Tc_pics/' +str(index) +'Tc' +'.png')

    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)
    fig.axes[0].set_xlim(5824, 6024)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(5924, c="#666666", zorder=-1)
    string="Barium star Tc5824"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Tc_pics/' +str(index) + 'Tcb' +'.png')

    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)
    string="Barium star full"+starID[index] + "Index: " + str(index)
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Tc_pics/'+str(index)+'zfull'+'.png')
    plt.close("all")

    del fig
    print(index)

for index, star in enumerate(csv):
    if not N_Tc[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]

    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)   
    fig.axes[0].set_xlim(5789, 5989)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(5889, c="#666666", zorder=-1)
    string="Na5889"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Na_pics/' +str(index) +'Na' +'.png')
 
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)   
    fig.axes[0].set_xlim(4000, 4600)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(4077, c="#666666", zorder=-1)
    fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    string="Ba4554 and St4077"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Na_pics/' +str(index) +'StBa' +'.png')

    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)
    string="Barium star full"+starID[index] + "Index: " + str(index)
    fig.suptitle(string)
    fig.savefig('Goodfigures/Bariumstar_Na_pics/'+str(index)+'zfull'+'.png')
    plt.close("all")

    del fig
    print(index)
"""

print("Barium stars: {}".format(sum(N_Bariumstars)))
print("Sodium 5889: {}".format(sum(N_Na5889)))
print("Technetium 4049: {}".format(sum(N_Tc4049)))
print("Technetium 4238: {}".format(sum(N_Tc4238)))
print("Technetium 4262: {}".format(sum(N_Tc4262)))
print("Technetium 4297: {}".format(sum(N_Tc4297)))
print("Technetium 5924: {}".format(sum(N_Tc5924)))
print("Technetium all: {}".format(sum(N_TcAll)))
print("Barium star Tc Matches: {}".format(sum(N_Tc)))
print("Barium star Na Matches: {}".format(sum(N_Na)))

