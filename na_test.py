import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib.ticker as mtick

import lamost
import utils


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

candidates = Table.read("Tc_catalog.fits")

#Barium
Ba4554_wavelength = candidates['Ba4554_wavelength']
Ba4554_amp = candidates['Ba4554_amplitude']
Ba4554_amperr = candidates['Ba4554_amplitude_err']

#Strontium
St4077_wavelength = candidates['St4077_wavelength']
St4077_amp = candidates['St4077_amplitude']
St4077_amperr = candidates['St4077_amplitude_err']

#Sodium
Na5889_wavelength = candidates['Na5889_wavelength']
Na5889_amp = candidates['Na5889_amplitude']
Na5889_amperr = candidates['Na5889_amplitude_err']


snrg = candidates['snrg']
chi = candidates['cannon_red_chisq']
starID = candidates['id']
teff = candidates['cannon_teff']
surfg = candidates['cannon_logg']
met = candidates['cannon_alpha_m']

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


#Sodium
N_Na5889 = (np.abs(Na5889_amp/(Na5889_amperr.astype(float) + 1e-10)) > 3) \
       * (chi < 3) \
       * (Na5889_amp < -0.05) \
       * (np.abs(Na5889_wavelength - 5889) < 2) \
       * (Na5889_amperr.astype(float) > 0)\
       * (snrg > 30)

N_Na = N_St4077*N_Ba4554*N_Na5889

print("Barium star Na Matches: {}".format(sum(N_Na)))

'''
for index, star in enumerate(candidates):
    if not N_Na[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)   
    fig.axes[0].set_xlim(5789, 5989)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(5889, c="#666666", zorder=-1)
    string="Na5889"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Na_figure/' +str(index) +'Na' +'.png')
    
    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)   
    fig.axes[0].set_xlim(4000, 4600)
    fig.axes[1].set_ylim(0.7, 1.2)
    fig.axes[1].axvline(4077, c="#666666", zorder=-1)
    fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    string="Ba4554 and St4077"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    fig.suptitle(string)
    fig.savefig('Na_figure/' +str(index) +'StBa' +'.png')

    fig = utils.plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux)
    fig.axes[1].axvline(4077, c="#666666", zorder=-1)
    fig.axes[1].axvline(4554, c="#666666", zorder=-1)
    fig.axes[1].axvline(5889, c="#666666", zorder=-1)
    string="Barium star full"+starID[index] + "Index: " + str(index)
    fig.suptitle(string)
    fig.savefig('Na_figure/'+str(index)+'zfull'+'.png')
    
    plt.close("all")
    
    del fig
    print(index)
    print(star['id'])
'''
ID = ['GAC116N20B', 'VB061N34V2', 'HD074049N040418V01', 
      'HD145704N052550B', 'HD172251S022317B01']
extinction = [0.0546, 0.2898, 0.0465, 0.0343, 0.4382]
counter = 0
fig, ax = plt.subplots(figsize=(5,12), nrows=5,ncols=1, sharex=True)
for index, star in enumerate(candidates):
    if not N_Na[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]

    ax[counter] = utils.plot_spectrummod(wavelengths, observed_flux, observed_ivar, model_flux, passed_ax=ax[counter])   
    #ax[counter].set_title('Index: '+str(index)+' E(B-V)='+
                                #str(extinction[counter]))
    ax[counter].annotate(str(ID[counter]),
                (5985, 1.15), fontsize=10, 
                horizontalalignment='right', verticalalignment='top')
    ax[counter].annotate('E(B-V)='+ str(extinction[counter]),
                (5985, 0.75), fontsize=10, 
                horizontalalignment='right', verticalalignment='bottom')
    #string="Na5889"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    #fig.suptitle(string)
    #fig.savefig('Na_figure/' +str(index) +'Natest' +'.png')
    #ax.set_xlim(5789, 5989)
    #ax.set_ylim(0.7, 1.2)
    #ax.axvline(5889, c="#666666", zorder=-1)
    
    #plt.close("all")
    
    #del fig
    print(index)
    print(counter)
    #print(starID[index])
    counter = counter + 1
    plt.close()
fig.text(0.5, 0.04, 'Wavelengths (Angstroms)', ha='center', va='center')
fig.text(0.04, 0.5, 'Normalized flux', ha='center', va='center', rotation='vertical')
#ax.set_xlabel("Wavelengths (Angstroms)")
#ax.set_ylabel("Normalized flux")
fig.savefig("Na_figure/na_figure.png", dpi=300, overwrite=True)
fig.savefig("Na_figure/na_figure.pdf", dpi=300, overwrite=True)
fig.show()


'''
for index, star in enumerate(candidates):
    if not N_Na[index]: continue

    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig, ax = utils.plot_spectrummod(wavelengths, observed_flux, observed_ivar, model_flux)   
    ax.set_xlim(5789, 5989)
    ax.set_ylim(0.7, 1.2)
    ax.axvline(5889, c="#666666", zorder=-1)
    string="Na5889"+starID[index] + "Index: " + str(index)+" Teff: "+str(teff[index])+" SurfG: "+str(surfg[index])+" M: "+str(met[index])
    #fig.suptitle(string)
    fig.savefig('Na_figure/' +str(index) +'Natest' +'.png')
    
    
    plt.close("all")
    
    del fig
    print(index)
extinction = [0.0546, 0.2898, 0.0465, 0.0343, 0.4382]
counter = 0
figure, axes = plt.subplots(5, 1)
for index, star in enumerate(candidates):
    if not N_Na[index]: continue
    
    observed_flux = all_observed_flux[index]
    observed_ivar = all_observed_ivar[index]
    model_flux = all_model_flux[index]
    
    fig = utils.plot_spectrummod(wavelengths, observed_flux, observed_ivar, model_flux)   
    #fig.axes[0].set_xlim(5789, 5989)
    fig.axes[counter].set_ylim(0.7, 1.2)
    fig.axes[counter].axvline(5889, c="#666666", zorder=-1)
    fig.axes[counter].set_title('Index: '+str(index)+' E(B-V)='+
                                str(extinction[counter]))
    plt.close("all")
    
    del fig
    print(counter)
    counter = counter +1
    


string="Na enhanced candidates"
figure.suptitle(string)
figure.savefig('Na_figure/na_fig.png')
'''





