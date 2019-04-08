
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

import lamost

def fill_between_steps(ax, x, y1, y2=0, h_align='mid', **kwargs):
    """
    Fill between for step plots in matplotlib.

    **kwargs will be passed to the matplotlib fill_between() function.
    """

    # If no Axes opject given, grab the current one:

    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = np.repeat((x[1:] - x[:-1]), 2)
    xstep = np.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = np.append(xx, xx.max() + xstep[-1])

    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)#[:-1]
    if type(y2) == np.ndarray:
        y2 = y2.repeat(2)#[:-1]

    # now to the plotting part:
    return ax.fill_between(xx, y1, y2=y2, **kwargs)



def plot_spectrum(wavelengths, observed_flux, observed_ivar, model_flux):

    gs = GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(13, 4))
    ax_residual = plt.subplot(gs[0])
    ax_spectrum = plt.subplot(gs[1], sharex=ax_residual)

    observed_flux_error = observed_ivar**-2
    ax_spectrum.plot(wavelengths, observed_flux, c="k", drawstyle="steps-mid")
    fill_between_steps(ax_spectrum, wavelengths, 
        observed_flux - observed_flux_error, observed_flux + observed_flux_error,
        facecolor="k", alpha=0.5)

    model_flux_error = lamost.scatters
    ax_spectrum.plot(wavelengths, model_flux, c="r", drawstyle="steps-mid")
    fill_between_steps(ax_spectrum, wavelengths, model_flux - model_flux_error,
        model_flux + model_flux_error, facecolor="r", alpha=0.5)

    residual_flux = observed_flux - model_flux
    residual_flux_error = np.sqrt(model_flux_error**2 + observed_flux_error**2)

    ax_residual.plot(wavelengths, residual_flux, c="k", drawstyle="steps-mid")
    fill_between_steps(ax_residual, wavelengths,
        residual_flux - residual_flux_error, residual_flux + residual_flux_error,
        facecolor="k", alpha=0.5)

    for ax in (ax_spectrum, ax_residual):
        ax.set_xlim(wavelengths[0], wavelengths[-1])

    ax_spectrum.set_ylim(0, 1.2)

    value = np.mean(np.abs(np.percentile(residual_flux, [1, 99])))
    ax_residual.set_ylim(-value, +value)

    # Hide x-label ticks on the residual axes
    plt.setp(ax_residual.get_xticklabels(), visible=False)

    ax_spectrum.set_xlabel("Wavelengths (Angstroms)")
    ax_spectrum.set_ylabel("Normalized flux")
    \
    fig.tight_layout()

    return fig

def plot_spectrummod(wavelengths, observed_flux, observed_ivar, model_flux, passed_ax=None):

    #gs = GridSpec(2, 1, height_ratios=[1, 4])
    #fig, ax = plt.subplots(figsize=(13, 4))

    if passed_ax:
      ax = passed_ax
    else:
      fig = plt.figure()
      ax = fig.add_subplot(111)
    #fig = plt.subplots(figsize=(13, 4))
    #ax = fig.add_subplot(111)
    #ax_residual = plt.subplot(gs[0])
    #ax_spectrum = plt.subplot(gs[1], sharex=ax_residual)

    observed_flux_error = observed_ivar**-2
    ax.plot(wavelengths, observed_flux, c="k", drawstyle="steps-mid")
    fill_between_steps(ax, wavelengths, 
        observed_flux - observed_flux_error, observed_flux + observed_flux_error,
        facecolor="k", alpha=0.5)

    model_flux_error = lamost.scatters
    ax.plot(wavelengths, model_flux, c="r", drawstyle="steps-mid")
    fill_between_steps(ax, wavelengths, model_flux - model_flux_error,
        model_flux + model_flux_error, facecolor="r", alpha=0.5)

    #residual_flux = observed_flux - model_flux
    #residual_flux_error = np.sqrt(model_flux_error**2 + observed_flux_error**2)

    #ax_residual.plot(wavelengths, residual_flux, c="k", drawstyle="steps-mid")
    #fill_between_steps(ax_residual, wavelengths,
        #residual_flux - residual_flux_error, residual_flux + residual_flux_error,
        #facecolor="k", alpha=0.5)

    #for ax in (ax_spectrum, ax_residual):
        #ax.set_xlim(wavelengths[0], wavelengths[-1])
    
    #ax.set_xlim(wavelengths[0], wavelengths[-1])
    ax.set_xlim(5789, 5989)
    ax.set_ylim(0.7, 1.2)
    ax.axvline(5889, c="#666666", zorder=-1)


    #value = np.mean(np.abs(np.percentile(residual_flux, [1, 99])))
    #ax_residual.set_ylim(-value, +value)

    # Hide x-label ticks on the residual axes
    #plt.setp(ax.get_xticklabels(), visible=False)


    #\
    #ax.tight_layout()

    if passed_ax:
      return ax
    else:
     fig.show()
