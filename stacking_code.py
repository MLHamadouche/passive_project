import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import spectres
import LoadData as ld
import matplotlib.pyplot as plt

"""
Outside of the function, define the pandas list of IDs that is needed to index the redshifts on line 49 (after for loop).
The function takes the argument 'objects_list' whcih is just a list of strings of the IDs of spectra to be loaded
by the LoadData function (imported as ld.load_vandels_spectra), or whichever loading data function as long as it takes the same
ID format (e.g. 'CDFS-HST034930SELECT') and outputs a 3d array, where spectrum[:,0]= wavelengths, spectrum[:,1] = fluxes,
spectrum[:,2] = flux errors. new_wavs is defined inside the function but can be modified.
"""
#both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()
#df1 = pd.DataFrame(both_xray)
#passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
#df2 = pd.DataFrame(passive_cut)
#new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)

####ID list for indexing redshifts########
#ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))

####ist of objects - creating the list of IDs in correct format from the pandas dataframe######
#objects1  = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))
#all_obs = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))
#objects_list = list(set(all_obs).difference(objects1)) #<- list of IDs####


concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df = pd.DataFrame(concat_3dhst)

ID_ = df.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))


def stacks(objects_list): #input array of redshifts for given list of objects
    new_wavs = np.arange(2400, 4200, 1.25)
    new_spec = np.zeros(len(new_wavs))
    new_errs = np.zeros(len(new_wavs))

    old_spec =[]
    old_errs = []

    spec = []
    spec_err = []
    med_norm = []
    med_spec_units  = []

    median_spec = np.zeros(len(new_wavs))
    errs = np.zeros(len(new_wavs))
    spec_ = np.zeros(len(new_wavs))
    spec_errs = np.zeros(len(new_wavs))
    med_spec_units=[]
    med_spectrum =np.zeros(len(new_wavs))
    new_redshifts = []
    all_masses = []
    for ID in objects_list:
        masses = ID_.loc[ID, 'log10(M*)']
        all_masses.append(masses)
        z = ID_.loc[ID, 'zspec']
        new_redshifts.append(z)
        spectrum = ld.load_vandels_spectra(ID)
        wav_mask = (spectrum[:,0]>5200) & (spectrum[:,0]<9250)

        flux = spectrum[:,1][wav_mask]
        flux_errs = spectrum[:,2][wav_mask]
        wavs = spectrum[:,0][wav_mask]
        #plt.plot(wavs, flux)
        #plt.ylim(-2*10**-18, 2*10**-18)
        #plt.savefig(str(ID)+'.pdf')
        #plt.close()
        zeros_mask = (flux == 0.)|(flux_errs == 0.)
        flux[zeros_mask] = np.nan
        flux_errs[zeros_mask] = np.nan

        rest_wavs = wavs/(1.0 + z)
        mask =  (rest_wavs > 3000) & (rest_wavs < 3500) # fairly featureless region of spectrum
        old_spec = flux/np.nanmedian(flux[mask]) #normalisation median from that region
        old_errs = flux_errs/np.nanmedian(flux[mask])

        med_norm.append(np.nanmedian(flux[mask]))
        new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)
        spec.append(new_spec)
        spec_err.append(new_errs)



    spec = np.transpose(spec)
    spec_err = np.transpose(spec_err)
    standev_err = []
    med_new = np.nanmedian(med_norm)
    for m in range(len(new_wavs)):
        spec_ = spec[m,:]
        spec_errs = spec_err[m,:]
        #standev_err[m] = np.std(spec_, axis=0)
        median_spec[m]=np.nanmedian(spec_)

    med_spec_units = median_spec*med_new #test removing normalisations to see fluxes nh

    return median_spec #returns an array of the new median stacked fluxes



####example:#####
#med_stack_missing = stacks(objects_list)
