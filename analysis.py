import numpy as np
import os
import matplotlib.pyplot as plt

import lamost
import utils


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


# Plot a special star.
star_index = 93

observed_flux = all_observed_flux[star_index]
observed_ivar = all_observed_ivar[star_index]
model_flux = all_model_flux[star_index]

fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)





# Fit a Gaussian to the line at 6707 Angstroms.
x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
    observed_ivar, model_flux, 4534, 4574)

p0 = np.array([-0.1,  4554, 2])
p_opt, p_cov = lamost.fit_gaussian(x, y, y_err, p0)


fig, ax = plt.subplots()
ax.scatter(x, y, facecolor="k")
ax.errorbar(x, y, yerr=y_err, fmt=None, ecolor="k", alpha=0.5)

draws = np.random.multivariate_normal(p_opt, p_cov, size=30)
model_draws = np.array([lamost.gaussian(x, *draw) for draw  in draws])
draw_between = np.percentile(model_draws, [16, 84], axis=0)

ax.plot(x, lamost.gaussian(x, *p_opt), c='r')
ax.fill_between(x, draw_between[0], draw_between[1],
    facecolor="r", alpha=0.5)

# Record p_opt and the uncertainty in each parameter.
p_opt_error = np.diag(np.sqrt(p_cov))





L = len(p_opt)
LI_p_optsLi = np.zeros((N, L))
LI_p_covsLi = np.zeros((N, L, L))
LI_p_optsBa = np.zeros((N, L))
LI_p_covsBa = np.zeros((N, L, L))
LI_p_optsEu = np.zeros((N, L))
LI_p_covsEu = np.zeros((N, L, L))


for index, star in enumerate(catalog[0:10000]):

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]

#Lithium
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 6700, 6720)

    p0 = np.array([-0.1, 6707, 2])
    p_optLi, p_covLi = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorLi = np.sqrt(np.diag(p_cov))

#Barium 4554
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4534, 4574)

    p0 = np.array([-0.15, 4554, 2])
    p_optBa, p_covBa = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa = np.sqrt(np.diag(p_cov))
   

#Europium 4129 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4109, 4149)

    p0 = np.array([-0.15, 4129, 2])
    p_optEu, p_covEu = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu = np.sqrt(np.diag(p_cov))




    
    # try >2 sigma and without abs(p_opt[0]) constraint
    #if np.abs(p_opt[0]/p_opt_error[0]) > 5 and abs(p_opt[0]) > 0.025:
        #fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            #model_flux)

        #fig.axes[0].set_xlim(4000, 5000)
        #fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    
        

    LI_p_optsLi[index] = p_optLi
    LI_p_covsLi[index] = p_covLi
    LI_p_optsBa[index] = p_optBa
    LI_p_covsBa[index] = p_covBa
    LI_p_optsEu[index] = p_optEu
    LI_p_covsEu[index] = p_covEu

    print(index)


LI_p_opt_errLi = np.array([np.sqrt(np.diag(each)) for each in LI_p_covsLi])
LI_p_opt_errBa = np.array([np.sqrt(np.diag(each)) for each in LI_p_covsBa])
LI_p_opt_errEu = np.array([np.sqrt(np.diag(each)) for each in LI_p_covsEu])

# Save these.
for index, label_name in enumerate(("amplitude", "wavelength", "sigma")[0:10000]):
    catalog["LI_6707_{}".format(label_name)] = LI_p_optsLi[:, index]
    catalog["LI_6707_{}_err".format(label_name)] = LI_p_opt_errLi[:, index]
    catalog["BA_4554_{}".format(label_name)] = LI_p_optsBa[:, index]
    catalog["BA_4554_{}_err".format(label_name)] = LI_p_opt_errBa[:, index]
    catalog["EU_4129_{}".format(label_name)] = LI_p_optsEu[:, index]
    catalog["EU_4129_{}_err".format(label_name)] = LI_p_opt_errEu[:, index]


catalog.write("my-catalog.csv")



