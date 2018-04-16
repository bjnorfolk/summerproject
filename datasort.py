import numpy as np
import os
import matplotlib.pyplot as plt

import lamost
import utils
import pandas as pd


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


Li6707_amp, Li6707_amperr, Ba4554_amp, Ba4554_amperr, Ba4934_amp, Ba4934_amperr, Ba5853_amp, Ba5853_amperr, Ba6141_amp, Ba6141_amperr, Ba6496_amp, Ba6496_amperr, Eu4129_amp, Eu4129_amperr, Eu4205_amp, Eu4205_amperr, Eu4435_amp, Eu4435_amperr, Eu4522_amp, Eu4522_amperr, Eu6645_amp, Eu6645_amperr, St4077_amp, St4077_amperr, St4215_amp, St4215_amperr = pd.read_csv("my-catalog.csv", index_col=13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38)
    

for index, star in enumerate(catalog[0:10000]):
    
    observed_flux = all_observed_flux[star_index]
    observed_ivar = all_observed_ivar[star_index]
    model_flux = all_model_flux[star_index]


    #try >2 sigma and without abs(p_opt[0]) constraint
    if np.abs(p_opt[0]/p_opt_error[0]) > 5 and p_opt[0] < -0.025:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            model_flux)

        fig.axes[0].set_xlim(4000, 5000)
        fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    
