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
Li_p_opts = np.zeros((N, L))
Li_p_covs = np.zeros((N, L, L))

Ba4554_p_opts = np.zeros((N, L))
Ba4554_p_covs = np.zeros((N, L, L))
Ba4934_p_opts = np.zeros((N, L))
Ba4934_p_covs = np.zeros((N, L, L))
Ba5853_p_opts = np.zeros((N, L))
Ba5853_p_covs = np.zeros((N, L, L))
Ba6141_p_opts = np.zeros((N, L))
Ba6141_p_covs = np.zeros((N, L, L))
Ba6496_p_opts = np.zeros((N, L))
Ba6496_p_covs = np.zeros((N, L, L))

Eu4129_p_opts = np.zeros((N, L))
Eu4129_p_covs = np.zeros((N, L, L))
Eu4205_p_opts = np.zeros((N, L))
Eu4205_p_covs = np.zeros((N, L, L))
Eu4435_p_opts = np.zeros((N, L))
Eu4435_p_covs = np.zeros((N, L, L))
Eu4522_p_opts = np.zeros((N, L))
Eu4522_p_covs = np.zeros((N, L, L))
Eu6645_p_opts = np.zeros((N, L))
Eu6645_p_covs = np.zeros((N, L, L))

St4077_p_opts = np.zeros((N, L))
St4077_p_covs = np.zeros((N, L, L))
St4215_p_opts = np.zeros((N, L))
St4215_p_covs = np.zeros((N, L, L))


for index, star in enumerate(catalog[0:1000]):

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
        observed_ivar, model_flux, 4544, 4564)

    p0 = np.array([-0.15, 4554, 2])
    p_optBa4554, p_covBa4554 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa4554 = np.sqrt(np.diag(p_cov))

#Barium 4934
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4924, 4944)

    p0 = np.array([-0.15, 4934, 2])
    p_optBa4934, p_covBa4934 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa4934 = np.sqrt(np.diag(p_cov))

#Barium 5853
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 5843, 5863)

    p0 = np.array([-0.15, 5853, 2])
    p_optBa5853, p_covBa5853 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa5853 = np.sqrt(np.diag(p_cov))

#Barium 6141
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 6131, 6151)

    p0 = np.array([-0.15, 6141, 2])
    p_optBa6141, p_covBa6141 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa = np.sqrt(np.diag(p_cov))

#Barium 6496.9
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 6486.9, 6506.9)

    p0 = np.array([-0.15, 6496.9, 2])
    p_optBa6496, p_covBa6496 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorBa6496 = np.sqrt(np.diag(p_cov))
   

#Europium 4129 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4109, 4149)

    p0 = np.array([-0.15, 4129, 2])
    p_optEu4129, p_covEu4129 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu4129 = np.sqrt(np.diag(p_cov))

#Europium 4205 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4195, 4215)

    p0 = np.array([-0.15, 4205, 2])
    p_optEu4205, p_covEu4205 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu4205 = np.sqrt(np.diag(p_cov))

#Europium 4435 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4425, 4445)

    p0 = np.array([-0.15, 4435, 2])
    p_optEu4435, p_covEu4435 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu4435 = np.sqrt(np.diag(p_cov))

#Europium 4522 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4512, 4532)

    p0 = np.array([-0.15, 4522, 2])
    p_optEu4522, p_covEu4522 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu4522 = np.sqrt(np.diag(p_cov))

#Europium 6645 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 6635, 6655)

    p0 = np.array([-0.15, 6645, 2])
    p_optEu6645, p_covEu6645 = lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorEu6645 = np.sqrt(np.diag(p_cov))

#Strontium 4077.71 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4067, 4087)

    p0 = np.array([-0.15, 4077.71, 2])
    p_optSt4077, p_covSt4077= lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorSt4077 = np.sqrt(np.diag(p_cov))

#Strontium 4215.52 
    x, y, y_err, indices = lamost.get_data_to_fit(wavelengths, observed_flux,
        observed_ivar, model_flux, 4205, 4225)

    p0 = np.array([-0.15, 4215.52 , 2])
    p_optSt4215, p_covSt4215= lamost.fit_gaussian(x, y, y_err, p0)
    p_opt_errorSt4215 = np.sqrt(np.diag(p_cov))



    
    # try >2 sigma and without abs(p_opt[0]) constraint
    #if np.abs(p_opt[0]/p_opt_error[0]) > 5 and abs(p_opt[0]) > 0.025:
        #fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar,
            #model_flux)

        #fig.axes[0].set_xlim(4000, 5000)
        #fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    
        

    Li_p_opts[index] = p_optLi
    Li_p_covs[index] = p_covLi

    Ba4554_p_opts[index] = p_optBa4554
    Ba4554_p_covs[index] = p_covBa4554
    Ba4934_p_opts[index] = p_optBa4934
    Ba4934_p_covs[index] = p_covBa4934
    Ba5853_p_opts[index] = p_optBa5853
    Ba5853_p_covs[index] = p_covBa5853
    Ba6141_p_opts[index] = p_optBa6141
    Ba6141_p_covs[index] = p_covBa6141
    Ba6496_p_opts[index] = p_optBa6496
    Ba6496_p_covs[index] = p_covBa6496

    Eu4129_p_opts[index] = p_optEu4129
    Eu4129_p_covs[index] = p_covEu4129
    Eu4205_p_opts[index] = p_optEu4205
    Eu4205_p_covs[index] = p_covEu4205
    Eu4435_p_opts[index] = p_optEu4435
    Eu4435_p_covs[index] = p_covEu4435
    Eu4522_p_opts[index] = p_optEu4522
    Eu4522_p_covs[index] = p_covEu4522
    Eu6645_p_opts[index] = p_optEu6645
    Eu6645_p_covs[index] = p_covEu6645
    Eu6645_p_opts[index] = p_optEu6645
    Eu6645_p_covs[index] = p_covEu6645

    St4077_p_opts[index] = p_optSt4077
    St4077_p_covs[index] = p_covSt4077
    St4077_p_opts[index] = p_optSt4077
    St4077_p_covs[index] = p_covSt4077
    St4215_p_opts[index] = p_optSt4215
    St4215_p_covs[index] = p_covSt4215




    print(index)


Li_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Li_p_covs])

Ba4554_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Ba4554_p_covs])
Ba4934_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Ba4934_p_covs])
Ba5853_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Ba5853_p_covs])
Ba6141_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Ba6141_p_covs])
Ba6496_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Ba6496_p_covs])

Eu4129_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Eu4129_p_covs])
Eu4205_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Eu4205_p_covs])
Eu4435_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Eu4435_p_covs])
Eu4522_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Eu4522_p_covs])
Eu6645_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in Eu6645_p_covs])

St4077_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in St4077_p_covs])
St4215_p_opt_err = np.array([np.sqrt(np.diag(each)) for each in St4215_p_covs])


# Save these.
for index, label_name in enumerate(("amplitude", "wavelength", "sigma")[0:1000]):
    catalog["LI_6707_{}".format(label_name)] = Li_p_opts[:, index]
    catalog["LI_6707_{}_err".format(label_name)] = Li_p_opt_err[:, index]

    catalog["BA_4554_{}".format(label_name)] = Ba4554_p_opts[:, index]
    catalog["BA_4554_{}_err".format(label_name)] = Ba4554_p_opt_err[:, index]
    catalog["BA_4934_{}".format(label_name)] = Ba4934_p_opts[:, index]
    catalog["BA_4934_{}_err".format(label_name)] = Ba4934_p_opt_err[:, index]
    catalog["BA_5853_{}".format(label_name)] = Ba5853_p_opts[:, index]
    catalog["BA_5853_{}_err".format(label_name)] = Ba5853_p_opt_err[:, index]
    catalog["BA_6141_{}".format(label_name)] = Ba6141_p_opts[:, index]
    catalog["BA_6141_{}_err".format(label_name)] = Ba6141_p_opt_err[:, index]
    catalog["BA_6496_{}".format(label_name)] = Ba6496_p_opts[:, index]
    catalog["BA_6496_{}_err".format(label_name)] = Ba6496_p_opt_err[:, index]

    catalog["EU_4129_{}".format(label_name)] = Eu4129_p_opts[:, index]
    catalog["EU_4129_{}_err".format(label_name)] = Eu4129_p_opt_err[:, index]
    catalog["EU_4205_{}".format(label_name)] = Eu4205_p_opts[:, index]
    catalog["EU_4205_{}_err".format(label_name)] = Eu4205_p_opt_err[:, index]
    catalog["EU_4435_{}".format(label_name)] = Eu4435_p_opts[:, index]
    catalog["EU_4435_{}_err".format(label_name)] = Eu4435_p_opt_err[:, index]
    catalog["EU_4522_{}".format(label_name)] = Eu4522_p_opts[:, index]
    catalog["EU_4522_{}_err".format(label_name)] = Eu4522_p_opt_err[:, index]
    catalog["EU_6645_{}".format(label_name)] = Eu6645_p_opts[:, index]
    catalog["EU_6645_{}_err".format(label_name)] = Eu6645_p_opt_err[:, index]

    catalog["ST_4077_{}".format(label_name)] = St4077_p_opts[:, index]
    catalog["ST_4077_{}_err".format(label_name)] = St4077_p_opt_err[:, index]
    catalog["ST_4215_{}".format(label_name)] = St4215_p_opts[:, index]
    catalog["ST_4215_{}_err".format(label_name)] = St4215_p_opt_err[:, index]


catalog.write("my-catalog.csv")



