import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.table as Table
import LoadData as ld

#input median stacked fluxes,
wavs = np.arange(2400, 4200, 1.25)

def Dn4000(fluxes):

    #indices for Dn4000 for rest frame wavelenghts of spectra
    #measure the strength of the 4000A break (Dn4000)
    mask_3850 = (wavs <= 3950) & (wavs >= 3850) #blue continuum
    flux_3850 = np.nanmean(fluxes[mask_3850])
    #print(flux_3850)
    mask_4000 = (wavs >= 4000) & (wavs <= 4100) #red continuum
    flux_4000 = np.nanmean(fluxes[mask_4000])
    #print(flux_4000)
    D4000_index = flux_4000/flux_3850


    return D4000_index


def C29_33(fluxes):
    #indices for the C(29-33) line for rest frame wavelenghts of spectra
    mask_2900 = (wavs <= 3100) & (wavs >= 2700)
    flux_2900 = np.nanmean(fluxes[mask_2900])
    mask_3300 = (wavs <= 3500) & (wavs >= 3100)
    flux_3300 = np.nanmean(fluxes[mask_3300])

    flux_ratio = flux_2900/flux_3300
    C_ind = -2.5*np.log10(flux_ratio)


    return C_ind


def H_delta():
    #from Brolagh.M, 1999
    mask_blue = (wavs >4030 )& (wavs <4082) #blue continuum
    mask_red = (wavs > 4122)& (wavs<4170) #red continuum
    line = (wavs>4082) & (wavs< 4122)
    return H_delt_ind #negative if emission, postive if absorption

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
