import numpy as np
import os
import matplotlib.pyplot as plt

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

data=pd.read_csv('mycatalog.csv')

csv=np.array(data)

Li6707_amp = csv[:,13]
Li6707_amperr = csv[:,14]

Ba4554_amp = csv[:,15]
Ba4554_amperr = csv[:,16]
Ba4934_amp = csv[:,17]
Ba4934_amperr = csv[:,18]
Ba5853_amp = csv[:,19]
Ba5853_amperr = csv[:,20]
Ba6141_amp = csv[:,21]
Ba6141_amperr = csv[:,22]
Ba6496_amp = csv[:,23]
Ba6496_amperr = csv[:,24]

Eu4129_amp = csv[:,25]
Eu4129_amperr = csv[:,26]
Eu4205_amp = csv[:,27]
Eu4205_amperr = csv[:,28]
Eu4435_amp = csv[:,29]
Eu4435_amperr = csv[:,30]
Eu4522_amp = csv[:,31]
Eu4522_amperr = csv[:,32]
Eu6645_amp = csv[:,33]
Eu6645_amperr = csv[:,34]

St4077_amp = csv[:,35]
St4077_amperr = csv[:,36]
St4215_amp = csv[:,37]
St4215_amperr = csv[:,38] 

starID=csv[:,0]

for index, star in enumerate(csv[0:10000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Li6707_amperr[index] ==0:
        continue
    if np.abs(Li6707_amp[index]/Li6707_amperr[index]) > 5 and Li6707_amp[index] < -0.025:

        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(6000, 7000)
        fig.axes[1].axvline(6707, c="#666666", zorder=-1)
        string="Li6707"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Li6707_pics/' +str(index) +'.png')

    if Ba4554_amperr[index] ==0:
        continue
    if np.abs(Ba4554_amp[index]/Ba4554_amperr[index]) > 3 and Ba4554_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(4000, 5000)
        fig.axes[1].axvline(4554, c="#666666", zorder=-1)
        string="Ba4554"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba4554_pics/' +str(index) +'.png')

    if Ba4934_amperr[index] ==0:
        continue
    if np.abs(Ba4934_amp[index]/Ba4934_amperr[index]) > 3 and Ba4934_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(4500, 5500)
        fig.axes[1].axvline(4934, c="#666666", zorder=-1)
        string="Ba4934"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba4934_pics/' +str(index) +'.png')

    if Ba5853_amperr[index] ==0:
        continue
    if np.abs(Ba5853_amp[index]/Ba5853_amperr[index]) > 3 and Ba5853_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(5200, 6200)
        fig.axes[1].axvline(5853, c="#666666", zorder=-1)
        string="Ba5853"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba5853_pics/' +str(index) +'.png')

    if Ba6141_amperr[index] ==0:
        continue
    if np.abs(Ba6141_amp[index]/Ba6141_amperr[index]) > 3 and Ba6141_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(5500, 6500)
        fig.axes[1].axvline(6141, c="#666666", zorder=-1)
        string="Ba6141"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba6141_pics/' +str(index) +'.png')

    if Ba6496_amperr[index] ==0:
        continue
    if np.abs(Ba6496_amp[index]/Ba6496_amperr[index]) > 3 and Ba6496_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(6000, 7000)
        fig.axes[1].axvline(6496, c="#666666", zorder=-1)
        string="Ba6496"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba6496_pics/' +str(index) +'.png')

    if Eu4129_amperr[index] ==0:
        continue
    if np.abs(Eu4129_amp[index]/Eu4129_amperr[index]) > 3 and Eu4129_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(3500, 4500)
        fig.axes[1].axvline(4129, c="#666666", zorder=-1)
        string="Eu4129"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4129_pics/' +str(index) +'.png')

    if Eu4205_amperr[index] ==0:
        continue
    if np.abs(Eu4205_amp[index]/Eu4205_amperr[index]) > 3 and Eu4205_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(3500, 4500)
        fig.axes[1].axvline(4205, c="#666666", zorder=-1)
        string="Eu4205"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4205_pics/' +str(index) +'.png')

    if Eu4435_amperr[index] ==0:
        continue
    if np.abs(Eu4435_amp[index]/Eu4435_amperr[index]) > 3 and Eu4435_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(4000, 5000)
        fig.axes[1].axvline(4435, c="#666666", zorder=-1)
        string="Eu4435"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4435_pics/' +str(index) +'.png')

    if Eu4522_amperr[index] ==0:
        continue
    if np.abs(Eu4522_amp[index]/Eu4522_amperr[index]) > 3 and Eu4522_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(4000, 5000)
        fig.axes[1].axvline(4522, c="#666666", zorder=-1)
        string="Eu4522"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4522_pics/' +str(index) +'.png')

    if Eu6645_amperr[index] ==0:
        continue
    if np.abs(Eu6645_amp[index]/Eu6645_amperr[index]) > 3 and Eu6645_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(6000, 7000)
        fig.axes[1].axvline(6645, c="#666666", zorder=-1)
        string="Eu6645"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu6645_pics/' +str(index) +'.png')

    if St4077_amperr[index] ==0:
        continue
    if np.abs(St4077_amp[index]/St4077_amperr[index]) > 3 and St4077_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(3500, 4500)
        fig.axes[1].axvline(4077, c="#666666", zorder=-1)
        string="St4077"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/St4077_pics/' +str(index) +'.png')

    if St4215_amperr[index] ==0:
        continue
    if np.abs(St4215_amp[index]/St4215_amperr[index]) > 3 and St4215_amp[index] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(3500, 4500)
        fig.axes[1].axvline(4215, c="#666666", zorder=-1)
        string="St4215"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/St4215_pics/' +str(index) +'.png')
        print(index)

