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

Li6707_amp = csv[:,13]
Li6707_amperr = csv[:,14]

Ba4554_wavelength = csv[:,41]
Ba4554_amp = csv[:,15]
Ba4554_amperr = csv[:,16]
Ba4934_wavelength = csv[:,43]
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

St4077_wavelength = csv[:,61]
St4077_amp = csv[:,35]
St4077_amperr = csv[:,36]
St4215_wavelength = csv[:,63]
St4215_amp = csv[:,37]
St4215_amperr = csv[:,38] 

chi = csv[:,9]
starID=csv[:,0]

#Barium Numbers
N_Ba4554 = (np.abs(Ba4554_amp/(Ba4554_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba4554_amp < -0.20) \
       * (np.abs(Ba4554_wavelength - 4554) < 2)

print("Number matched by Ba 4554 filter: {}".format(sum(N_Ba4554)))

N_Ba4934 = (np.abs(Ba4934_amp/(Ba4934_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba4934_amp < -0.20) \
       * (np.abs(Ba4934_wavelength - 4934) < 2)

print("Number matched by Ba 4934 filter: {}".format(sum(N_Ba4934)))

N_Ba5853 = (np.abs(Ba5853_amp/(Ba5853_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba5853_amp < -0.20) \
       * (np.abs(Ba5853_wavelength - 5853) < 2)

print("Number matched by Ba 5853 filter: {}".format(sum(N_Ba5853)))

N_Ba6141 = (np.abs(Ba6141_amp/(Ba6141_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba6141_amp < -0.20) \
       * (np.abs(Ba6141_wavelength - 6141) < 2)

print("Number matched by Ba 6141 filter: {}".format(sum(N_Ba6141)))

N_Ba6496 = (np.abs(Ba6496_amp/(Ba6496_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Ba6496_amp < -0.20) \
       * (np.abs(Ba6496_wavelength - 6496) < 2)

print("Number matched by Ba 6496 filter: {}".format(sum(N_Ba6496)))


#Europium Lines
N_Eu4129 = (np.abs(Eu4129_amp/(Eu4129_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4129_amp < -0.20) \
       * (np.abs(Eu4129_wavelength - 4129) < 2) \
       * (Eu4129_amperr.astype(float) > 0)

print("Number matched by Eu 4129 filter: {}".format(sum(N_Eu4129)))

N_Eu4205 = (np.abs(Eu4205_amp/(Eu4205_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4205_amp < -0.20) \
       * (np.abs(Eu4205_wavelength - 4205) < 2) \
       * (Eu4205_amperr.astype(float) > 0)

print("Number matched by Eu 4205 filter: {}".format(sum(N_Eu4205)))

N_Eu4435 = (np.abs(Eu4435_amp/(Eu4435_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4435_amp < -0.20) \
       * (np.abs(Eu4435_wavelength - 4435) < 2) \
       * (Eu4435_amperr.astype(float) > 0)

print("Number matched by Eu 4435 filter: {}".format(sum(N_Eu4435)))

N_Eu4522 = (np.abs(Eu4522_amp/(Eu4522_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu4522_amp < -0.20) \
       * (np.abs(Eu4522_wavelength - 4522) < 2) \
       * (Eu4522_amperr.astype(float) > 0)

print("Number matched by Eu 4522 filter: {}".format(sum(N_Eu4522)))


N_Eu6645 = (np.abs(Eu6645_amp/(Eu6645_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Eu6645_amp < -0.05) \
       * (np.abs(Eu6645_wavelength - 6645) < 2) \
       * (Eu6645_amperr.astype(float) > 0)

print("Number matched by Eu 6645 filter: {}".format(sum(N_Eu6645)))

#Stontium Lines
N_St4077 = (np.abs(St4077_amp/(St4077_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (St4077_amp < -0.2) \
       * (np.abs(St4077_wavelength - 4077) < 2) \
       * (St4077_amperr.astype(float) > 0)

print("Number matched by St 4077 filter: {}".format(sum(N_St4077)))

N_St4215 = (np.abs(St4215_amp/(St4215_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (St4215_amp < -0.2) \
       * (np.abs(St4215_wavelength - 4215) < 2) \
       * (St4215_amperr.astype(float) > 0)

print("Number matched by St 4215 filter: {}".format(sum(N_St4215)))

#Figure Printing

#Lithium
"""

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Li6707_amperr[index] ==0:
        continue
    if np.abs(Li6707_amp[index]/Li6707_amperr[index]) > 3 and chi[index] < 3 and Li6707_amp[index] < -0.025:

        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(6607, 6807)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(6707, c="#666666", zorder=-1)
        string="Li6707"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Li6707_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)
                
        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 6701, 6713)
    
        p0 = np.array([-0.1, 6707, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Li6707_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Li6707_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)
"""

#Barium
"""
 

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Ba4554_amperr[index] ==0:
        continue
    if np.abs(Ba4554_amp[index]/Ba4554_amperr[index]) > 3 and chi[index] < 3 and Ba4554_amp[index] < -0.20 and np.abs(Ba4554_wavelength[index] - 4554) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4454, 4654)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4554, c="#666666", zorder=-1)
        string="Ba4554"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba4554_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4548, 4560)
    
        p0 = np.array([-0.1, 4554, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Ba4554_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Ba4554_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)


for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Ba4934_amperr[index] ==0:
        continue
    if np.abs(Ba4934_amp[index]/Ba4934_amperr[index]) > 3 and chi[index] < 3 and Ba4934_amp[index] < -0.025 and np.abs(Ba4934_wavelength[index] - 4934) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4834, 5034)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4934, c="#666666", zorder=-1)
        string="Ba4934"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba4934_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4928, 4940)
    
        p0 = np.array([-0.1, 4934, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Ba4934_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Ba4934_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)


for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Ba5853_amperr[index] ==0:
        continue
    if np.abs(Ba5853_amp[index]/Ba5853_amperr[index]) > 3 and chi[index] < 3 and Ba5853_amp[index] < -0.025 and np.abs(Ba5853_wavelength[index] - 5853) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(5753, 5953)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(5853, c="#666666", zorder=-1)
        string="Ba5853"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba5853_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 5847, 5859)
    
        p0 = np.array([-0.1, 5853, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Ba5853_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Ba5853_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)



for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Ba6141_amperr[index] ==0:
        continue
    if np.abs(Ba6141_amp[index]/Ba6141_amperr[index]) > 3 and chi[index] < 3 and Ba6141_amp[index] < -0.025 and np.abs(Ba6141_wavelength[index] - 6141) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(6041, 6241)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(6141, c="#666666", zorder=-1)
        string="Ba6141"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba6141_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 6135, 6147)
    
        p0 = np.array([-0.1, 6141, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Ba6141_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Ba6141_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Ba6496_amperr[index] ==0:
        continue
    if np.abs(Ba6496_amp[index]/Ba6496_amperr[index]) > 3 and chi[index] < 3 and Ba6496_amp[index] < -0.025 and np.abs(Ba6496_wavelength[index] - 6496) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(6396, 6596)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(6496, c="#666666", zorder=-1)
        string="Ba6496"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Ba6496_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 6490, 6502)
    
        p0 = np.array([-0.1, 6496, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Ba6496_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Ba6496_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)


#Europium

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Eu4129_amperr[index] ==0:
        continue
    if np.abs(Eu4129_amp[index]/Eu4129_amperr[index]) > 3 and chi[index] < 3 and Eu4129_amp[index] < -0.20 and np.abs(Eu4129_wavelength[index] - 4129) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4029, 4229)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4129, c="#666666", zorder=-1)
        string="Eu4129"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4129_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4123, 4135)
    
        p0 = np.array([-0.1, 4129, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Eu4129_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Eu4129_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Eu4205_amperr[index] ==0:
        continue
    if np.abs(Eu4205_amp[index]/Eu4205_amperr[index]) > 3 and chi[index] < 3 and Eu4205_amp[index] < -0.025 and np.abs(Eu4205_wavelength[index] - 4205) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4105, 4305)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4205, c="#666666", zorder=-1)
        string="Eu4205"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4205_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4199, 4211)
    
        p0 = np.array([-0.1, 4205, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Eu4205_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Eu4205_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Eu4435_amperr[index] ==0:
        continue
    if np.abs(Eu4435_amp[index]/Eu4435_amperr[index]) > 3 and chi[index] < 3 and Eu4435_amp[index] < -0.025 and np.abs(Eu4435_wavelength[index] - 4435) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4335, 4535)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4435, c="#666666", zorder=-1)
        string="Eu4435"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4435_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4429, 4441)
    
        p0 = np.array([-0.1, 4435, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Eu4435_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Eu4435_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Eu4522_amperr[index] ==0:
        continue
    if np.abs(Eu4522_amp[index]/Eu4522_amperr[index]) > 3 and chi[index] < 3 and Eu4522_amp[index] < -0.025 and np.abs(Eu4522_wavelength[index] - 4522) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4422, 4622)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4522, c="#666666", zorder=-1)
        string="Eu4522"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu4522_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4516, 4528)
    
        p0 = np.array([-0.1, 4522, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Eu4522_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Eu4522_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)
"""
for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if Eu6645_amperr[index] ==0:
        continue
    if np.abs(Eu6645_amp[index]/Eu6645_amperr[index]) > 3 and chi[index] < 3 and Eu6645_amp[index] < -0.025 and np.abs(Eu6645_wavelength[index] - 6645) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(6545, 6745)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(6645, c="#666666", zorder=-1)
        string="Eu6645"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/Eu6645_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 6639, 6651)
    
        p0 = np.array([-0.1, 6645, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(Eu6645_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/Eu6645_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)
"""
#Strontium

for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if St4077_amperr[index] ==0:
        continue
    if np.abs(St4077_amp[index]/St4077_amperr[index]) > 3 and chi[index] < 3 and St4077_amp[index] < -0.025 and np.abs(St4077_wavelength[index] - 4077) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(3977, 4177)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4077, c="#666666", zorder=-1)
        string="St4077"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/St4077_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4071, 4083)
    
        p0 = np.array([-0.1, 4077, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(St4077_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/St4077_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)


for index, star in enumerate(csv[0:450000]):
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    

    if St4215_amperr[index] ==0:
        continue
    if np.abs(St4215_amp[index]/St4215_amperr[index]) > 3 and chi[index] < 3 and St4215_amp[index] < -0.025 and np.abs(St4215_wavelength[index] - 4215) < 2:
        fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)

        fig.axes[0].set_xlim(4115, 4315)
        fig.axes[1].set_ylim(0.7, 1.2)
        fig.axes[1].axvline(4215, c="#666666", zorder=-1)
        string="St4215"+starID[index] + "Index: " + str(index)
        fig.suptitle(string)
        fig.savefig('figures/St4215_pics/' +str(index) +'.png')
        plt.close("all")
        del fig
        print(index)

        x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux, observed_ivar, model_flux, 4209, 4221)
    
        p0 = np.array([-0.1, 4215, 2])
        p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)
        p_opt_error = np.sqrt(np.diag(p_cov))
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, facecolor="k")
        ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)
        try:
            draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
            model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
            draw_between = np.percentile(model_draws, [16, 84], axis=0)

            ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
            ax.fill_between(x, draw_between[0], draw_between[1],
                        facecolor="r", alpha=0.5)
        
            gaussstring= " Index: " + str(index) + "Gaussian Amp=" + str(St4215_amp[index])
            fig.suptitle(gaussstring)
            fig.savefig('figures/St4215_pics/'+str(index)+'GAUSS'+'.png')
        except ValueError :
            print("Error: Problem with drawing gaussian")
            continue
        plt.close("all")
        del fig
        print(index)
"""


#Barium Numbers

print("Number matched by Ba 4554 filter: {}".format(sum(N_Ba4554)))

print("Number matched by Ba 4934 filter: {}".format(sum(N_Ba4934)))

print("Number matched by Ba 5853 filter: {}".format(sum(N_Ba5853)))

print("Number matched by Ba 6141 filter: {}".format(sum(N_Ba6141)))

print("Number matched by Ba 6496 filter: {}".format(sum(N_Ba6496)))
#Europium Lines
print("Number matched by Eu 4129 filter: {}".format(sum(N_Eu4129)))

print("Number matched by Eu 4205 filter: {}".format(sum(N_Eu4205)))

print("Number matched by Eu 4435 filter: {}".format(sum(N_Eu4435)))

print("Number matched by Eu 4522 filter: {}".format(sum(N_Eu4522)))

print("Number matched by Eu 6645 filter: {}".format(sum(N_Eu6645)))

#Stontium Lines
print("Number matched by St 4077 filter: {}".format(sum(N_St4077)))

print("Number matched by St 4215 filter: {}".format(sum(N_St4215)))

