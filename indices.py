import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.table as Table
import LoadData as ld

#input median stacked fluxes,
wavs = np.arange(2400, 4200, 1.5)

def Dn4000(fluxes):
    #n is for narrow - less sensitive to reddening effects, used more in literature - 4000A to 4100A instead of 4050 - 4250
    #indices for Dn4000 for rest frame wavelenghts of spectra
    #measure the strength of the 4000A break (Dn4000)
    mask_3850 = (wavs <= 3950) & (wavs >= 3850) #blue continuum
    flux_3850 = np.nanmean(fluxes[mask_3850])
    #print(flux_3850)
    mask_4000 = (wavs >= 4000) & (wavs <= 4100) #red continuum
    flux_4000 = np.nanmean(fluxes[mask_4000])
    #print(flux_4000)
    Dn4000_index = flux_4000/flux_3850


    return Dn4000_index


def C29_33(fluxes):
    c = 3*10**8
    F_nu = fluxes*(wavs**2)/c
    #indices for the C(29-33) line for rest frame wavelenghts of spectra
    mask_2900 = (wavs <= 3100) & (wavs >= 2700)
    flux_2900 = np.nanmean(F_nu[mask_2900])
    mask_3300 = (wavs <= 3500) & (wavs >= 3100)
    flux_3300 = np.nanmean(F_nu[mask_3300])

    flux_ratio = flux_2900/flux_3300
    C_ind = -2.5*np.log10(flux_ratio)


    return C_ind


def H_delta(wavs, fluxes):
    #from Bolagh.M, 1999
    #changed to make either side 40Angstroms so no weighting change needed
    mask_blue = (wavs > 4042) & (wavs<4082) #blue continuum

    flux_blue = fluxes[mask_blue]
    mask_red = (wavs > 4122) & (wavs<4162) #red continuum
    flux_red = fluxes[mask_red]
    #flux_red = fluxes[mask_red]
    line = (wavs>4082) & (wavs < 4122)
    flux_line = fluxes[line]
    flux_continuum = np.mean(np.mean(flux_blue)+np.mean(flux_red)) #one number
    H_delta_EW = ((flux_continuum - np.mean(flux_line))*40)/flux_continuum
    #multiply by 40 angstroms / continuum flux
    # EW found by summing fluxes in region, and fitting with gaussian

    return H_delta_EW #negative if emission, postive if absorption
#mask_2825 = (wavs > 2725) & (wavs < 2825)
#print(wavs[mask_2825])


def Mg_UV(flux):
    mask_2625 = (wavs<2625)& (wavs > 2525)
    mask_2725 = (wavs > 2625) & (wavs < 2725)
    mask_2825 = (wavs > 2725) & (wavs < 2825)
    int_flux_2725 = np.trapz(flux[mask_2725], x = wavs[mask_2725])
    int_flux_2625 = np.trapz(flux[mask_2625], x = wavs[mask_2625])
    int_flux_2825 = np.trapz(flux[mask_2825], x = wavs[mask_2825])

    Mg_UV_index = (2 * int_flux_2725)/(int_flux_2625 + int_flux_2825)

    return Mg_UV_index
