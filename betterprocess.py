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

"""
d = Table.read("mycatalog.fits")
d.write("mycatalog.csv")
"""
data=pd.read_csv('mycatalog.csv')

csv=np.array(data)

Ba4554_wavelength = csv[:,41] #
Ba4554_amp = csv[:,15]
Ba4554_amperr = csv[:,16]
Ba4934_wavelength = csv[:,43] #
Ba4934_amp = csv[:,17]
Ba4934_amperr = csv[:,18]
Ba5853_wavelength = csv[:,45]
Ba5853_amp = csv[:,19]
Ba5853_amperr = csv[:,20]
Ba6141_wavelength = csv[:,47]
Ba6141_amp = csv[:,21]
Ba6141_amperr = csv[:,22]
Ba6496_wavelength = csv[:,49]
Ba6496_amp = csv[:,23]
Ba6496_amperr = csv[:,24]


Eu4129_wavelength = csv[:,51]
Eu4129_amp = csv[:,25]
Eu4129_amperr = csv[:,26]
Eu4205_wavelength = csv[:,53]
Eu4205_amp = csv[:,27]
Eu4205_amperr = csv[:,28]
Eu4435_wavelength = csv[:,55]
Eu4435_amp = csv[:,29]
Eu4435_amperr = csv[:,30]
Eu4522_wavelength = csv[:,57]
Eu4522_amp = csv[:,31]
Eu4522_amperr = csv[:,32]
Eu6645_wavelength = csv[:,59]
Eu6645_amp = csv[:,33]
Eu6645_amperr = csv[:,34]

St4077_wavelength = csv[:,61] #
St4077_amp = csv[:,35]
St4077_amperr = csv[:,36]
St4215_wavelength = csv[:,63] #
St4215_amp = csv[:,37]
St4215_amperr = csv[:,38] 

snrg = csv[:,12]
chi = csv[:,9]
starID=csv[:,0]

#Barium Numbers
N_Ba4554 = (np.abs(Ba4554_amp/(Ba4554_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba4554_amp < -0.05) \
       * (np.abs(Ba4554_wavelength - 4554) < 2)\
       * (snrg > 30)

N_Ba4934 = (np.abs(Ba4934_amp/(Ba4934_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba4934_amp < -0.05) \
       * (np.abs(Ba4934_wavelength - 4934) < 2)\
       * (snrg > 30)

N_Ba5853 = (np.abs(Ba5853_amp/(Ba5853_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba5853_amp < -0.05) \
       * (np.abs(Ba5853_wavelength - 5853) < 2)\
       * (snrg > 30)

N_Ba6141 = (np.abs(Ba6141_amp/(Ba6141_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba6141_amp < -0.05) \
       * (np.abs(Ba6141_wavelength - 6141) < 2)\
       * (snrg > 30)

N_Ba6496 = (np.abs(Ba6496_amp/(Ba6496_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba6496_amp < -0.05) \
       * (np.abs(Ba6496_wavelength - 6496) < 2)\
       * (snrg > 30)

N_Ba_doublematch = N_Ba4554*N_Ba4934
N_Ba_triplematch = N_Ba4554*N_Ba4934*N_Ba5853
N_Ba_quadmatch = N_Ba4554*N_Ba4934*N_Ba5853*N_Ba6141
N_Ba_multimatch = N_Ba4554*N_Ba4934*N_Ba5853*N_Ba6141*N_Ba6496

print("Number matched by Ba 4554 and Ba 4934: {}".format(sum(N_Ba_doublematch)))
print("Number matched by Ba 4554, 4934, 5853, 6141, 6496: {}".format(sum(N_Ba_multimatch)))

#Europium Lines
N_Eu4129 = (np.abs(Eu4129_amp/(Eu4129_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4129_amp < -0.05) \
       * (np.abs(Eu4129_wavelength - 4129) < 2) \
       * (Eu4129_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Eu4205 = (np.abs(Eu4205_amp/(Eu4205_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4205_amp < -0.05)\
       * (np.abs(Eu4205_wavelength - 4205) < 2) \
       * (Eu4205_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Eu4435 = (np.abs(Eu4435_amp/(Eu4435_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4435_amp < -0.05) \
       * (np.abs(Eu4435_wavelength - 4435) < 2) \
       * (Eu4435_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Eu4522 = (np.abs(Eu4522_amp/(Eu4522_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4522_amp < -0.05) \
       * (np.abs(Eu4522_wavelength - 4522) < 2) \
       * (Eu4522_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Eu6645 = (np.abs(Eu6645_amp/(Eu6645_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu6645_amp < -0.05) \
       * (np.abs(Eu6645_wavelength - 6645) < 2) \
       * (Eu6645_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Eu_doublematch = N_Eu4129*N_Eu4205
N_Eu_triplematch = N_Eu4129*N_Eu4205*N_Eu4435
N_Eu_quadmatch = N_Eu4129*N_Eu4205*N_Eu4435*N_Eu4522
N_Eu_multimatch = N_Eu4129*N_Eu4205*N_Eu4435*N_Eu4522*N_Eu6645

print("Number matched by Eu 4129 and Eu 4205: {}".format(sum(N_Eu_doublematch)))
print("Number matched by Eu 4129, 4205, 4435, 4552, 6645: {}".format(sum(N_Eu_multimatch)))

#Stontium Lines
N_St4077 = (np.abs(St4077_amp/(St4077_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (St4077_amp < -0.05) \
       * (np.abs(St4077_wavelength - 4077) < 2) \
       * (St4077_amperr.astype(float) > 0)\
       * (snrg > 30)

N_St4215 = (np.abs(St4215_amp/(St4215_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (St4215_amp < -0.05) \
       * (np.abs(St4215_wavelength - 4215) < 2) \
       * (St4215_amperr.astype(float) > 0)\
       * (snrg > 30)

N_St_doublematch = N_St4077*N_St4215


print("Number matched by St 4077 and St 4215: {}".format(sum(N_St_doublematch)))

"""
def plot_my_fucking_spectrum(index, xlims, ylims, show_wavelengths=None):

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

    fig.axes[0].set_xlim(xlims)
    fig.axes[1].set_ylim(ylims)
    if show_wavelengths is not None:
        show_wavelengths = np.atleast_1d(show_wavelengths)
        for each in show_wavelengths:
           fig.axes[1].axvline(each, c="#666666", zorder=-1)
    return fig

"""

for index, star in enumerate(csv):
    if not N_Eu_doublematch[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

    fig.axes[0].set_xlim(4029, 4229)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(4129, c="#666666", zorder=-1)
    string="Eu4129"+starID[index] + "Index: " + str(index)
    fig.suptitle(string)
    fig.savefig('Goodfigures/Eu_doublematch_pics/' +str(index) +'.png')
    
    fig.axes[0].set_xlim(4105, 4305)
    fig.axes[1].axvline(4205, c="#666666", zorder=-1)
    string="Eu4205 {} Index: {}".format(starID[index], str(index))
    fig.suptitle(string)
    fig.savefig("Goodfigures/Eu_doublematch_pics/{}.png".format(index))

    plt.close("all")

    del fig
    print(index)



print("Number matched by Ba 4554 and Ba 4934: {}".format(sum(N_Ba_doublematch)))
print("Number matched by Ba 4554, 4934, 5853: {}".format(sum(N_Ba_triplematch)))
print("Number matched by Ba 4554, 4934, 5853, 6141: {}".format(sum(N_Ba_quadmatch)))
print("Number matched by Ba 4554, 4934, 5853, 6141, 6496: {}".format(sum(N_Ba_multimatch)))


print("Number matched by Eu 4129 and Eu 4205: {}".format(sum(N_Eu_doublematch)))
print("Number matched by Eu 4129, 4205, 4435: {}".format(sum(N_Eu_triplematch)))
print("Number matched by Eu 4129, 4205, 4435, 4552: {}".format(sum(N_Eu_quadmatch)))
print("Number matched by Eu 4129, 4205, 4435, 4552, 6645: {}".format(sum(N_Eu_multimatch)))


print("Number matched by St 4077 and St 4215: {}".format(sum(N_St_doublematch)))
