#Greatly inspired by python-fsps demo specbymass.py
#Adapted by Lauren Taylor: lauren.y.taylor96@gmail.com (please email with Qs)

import numpy as np
import matplotlib.pyplot as plt
import fsps

def spectrum_masscuts(sps, masslims, age, plot_filters={}, plot=1, save=0):
    """Plots the spectral contributions from different mass cuts of the SPS plus the total spectum.
        
        Parameters
        __________
        sps : StellarPopulation
            A Stellar Population object from the fsps module. All parameters should be set.
        masslims : list
            A list of masses at which to make the various mass cuts.
        age : int, float
            The age of the stellar population at which to plot the spectra.
        plot_filters : dict
            A dict of pattern {fsps_name : name_to_plot}. Available filters can be found with fsps.list_filters().
        plot: bool
            Show plot of spectra at different mass cute (and filters if provided).
        save: bool
            Save the plot as a .png file.
            
    """
    # Set up stellar mass ranges and arrays to hold output spectra
    mspec = np.zeros([len(sps.wavelengths), len(masslims)])
    dmspec = np.zeros_like(mspec)
    wave, spec_tot = sps.get_spectrum(tage=age)

    dfig, dax = plt.subplots(figsize=(12, 7))
    
    for i,mc in enumerate(masslims):
        sps.params['masscut'] = mc #Change the mass cut to mc
        w, spec = sps.get_spectrum(tage=age) #Spectrum due to new mass cut
        mspec[:,i] = spec.copy()
        if i == 0:
            # Plot the lowest mass bin spectrum
            dax.plot(w, spec, label=f"0.08 < M < {mc}")
            # skip the rest of the loop
            continue

        # Subtract the last upper limit from the spectrum to get the
        #  spectrum due to stars within the bin
        dmspec[:, i] = spec - mspec[:, i-1]
        # Plot it
        label = f"{masslims[i-1]} < M < {mc}"
        dax.plot(w, dmspec[:, i], label=label)
        # Plot the total spectrum
        if mc == masslims[-1]:
            dax.plot(w, spec, label="Total", color='black')

# prettify the axes, titles, etc
    dax.legend(loc=0)
    dax.set_xlabel("$\lambda (\AA)$")
    dax.set_ylabel("$F_\lambda$")
    dax.set_yscale("log")
    dax.set_ylim(1e-30,1e-10)
    dax.set_xlim(0,10000)
    fstring = f"age = {age} Gyr, log$Z/Z_\odot$={sps.params['logzsol']}"
    dax.set_title(fstring)
    
    try:
        filters = [fsps.filters.get_filter(f) for f in plot_filters]
        for f in filters:
            wavelength, transmission = f.transmission
            tmax = transmission.max()
            wmax = wavelength[transmission.argmax()]
            dax.fill_between(wavelength, .5e-12,
                            alpha=0.3, color='grey')
            dax.text(wmax, 1e-12, plot_filters[f.name], fontdict={'size': 16})
    except TypeError:
        print('Filters must be a dict. Instead of filters="sdss_g" try filters={"sdss_g": "g"}.')

    
    if plot:
        plt.show()
    if save:
        plt.savefig(f"Age={age}_logzsol={sps.params['logzsol']}.png")

    return dfig

#TBC
def add_spectra(wavelengths_grids, spectra):
    #1. interpolate spectra onto the same wavelength grid if necessary
    #2. Add spectra in spectra list together
    #3. Plot against wavelength grid
    return

def main():
    #########################THICK DISK PARMETERS#################################
    age=11
    #The basic stellar population with constant metallicity from GALAH (Sharma et al. 2019), using Milky Way extinction law (dust_type=1)
    sps = fsps.StellarPopulation(zcontinuous=1, logzsol=-0.162, dust_type=1, masscut=120)

    #IMF parameters given in Sharma et al. 2019, Table 4
    sps.params['imf_type'] = 2
    sps.params['imf1'] = -0.5
    sps.params['imf2'] = -0.5
    sps.params['imf3'] = -0.5

    #SFH of a single 2 Gyr long burst of star-forming (Sharma et al. 2019, section 4.2.3), see also Robin et al. 2003 Table 1
    sps.params['sfh'] = 1
    sps.params['sf_start'] = 0
    sps.params['sf_trunc'] = 2
    sps.params['fburst'] = 1
    ##############################################################################
    
    # Get a few filter transmission curves, in case you want to call plot_filters
    filterlist = {'galex_fuv': 'FUV', 'galex_nuv': 'NUV', 'sdss_u': 'u',
        'sdss_g': 'g', 'sdss_i': 'i', '2mass_j': 'J'}

    #Chosen mass cuts for observing their spectral contributions
    masslims=[1.0, 3.0, 15.0, 30.0, 120.0]

    #Default values: show=1, save=0 - change save to 1 if you want save the image
    fig = spectrum_masscuts(sps, masslims, age)

if __name__ == "__main__":
    main()
