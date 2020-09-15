import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import os
import spectres
import LoadData as ld
import matplotlib.pyplot as plt
from matplotlib import gridspec
from glob import glob

"""
globpath = os.path.join('new_vandels_spec/', '*.fits')
filelist = glob(globpath)
hdulist = fits.open(filelist[-1])
delt_wav = hdulist[0].header['CDELT1']
wav_first_pixel = hdulist[0].header['CRVAL1']
"""

#passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
#df = pd.DataFrame(passive_cut)
#ID_list = df.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8"))
#redshifts = passive_cut['zspec']
passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
df = pd.DataFrame(passive_cut)
ID_list = df.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))


def stacks(new_wavs, objects):
    new_spec = np.zeros(len(new_wavs))
    new_errs = np.zeros(len(new_wavs))

    old_spec =[]
    old_errs = []

    new_objs = objects
    #resampling spectra onto new wavelength grid with spectres
    #z_mask = (redshifts <= 1.5) & (redshifts>= 0.9) # 89 objects
    #new_objs = objects[z_mask] #uncomment if want mask
    #print('lens', len(new_wavs), len(new_objs))
    specc = np.zeros((len(new_wavs),len(new_objs)))
    specerrs = np.zeros((len(new_wavs),len(new_objs)))
    spec = []
    spec_err = []
    med_norm = []
    for ID in new_objs:
        z = ID_list.loc[ID, 'zspec']
        spectrum = ld.load_vandels_spectra(ID)
        rest_wavs=spectrum[:,0]/(1.+z)
        mask = (rest_wavs < 3500) & (rest_wavs > 3000) # fairly featureless region of spectrum
        old_spec = spectrum[:,1]/np.nanmedian(spectrum[:,1][mask]) #normalisation median from that region
        old_errs = spectrum[:,2]/np.nanmedian(spectrum[:,2][mask])
        med_norm.append(np.nanmedian(spectrum[:,1][mask]))
        new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)
        spec.append(new_spec)
        spec_err.append(new_errs)

    spec = np.transpose(spec)
    spec_err = np.transpose(spec_err)

    median_spec = np.zeros(len(new_wavs))
    errs = np.zeros(len(new_wavs))
    spec_ = np.zeros(len(new_wavs))
    spec_errs = np.zeros(len(new_wavs))

    for m in range(len(new_wavs)):
        spec_ = spec[m,:]
        spec_errs = spec_err[m,:]
        median_spec[m]=np.nanmedian(spec_)
        errs[m] = np.nanmedian(spec_errs) #nopes
    #errs[m] = np.sqrt(1/np.sum(spec_errs**2)) #??? weighted errors
    #print(len(median_spec))
    med_new = np.nanmedian(med_norm) # multiply stacked median spectrum by this value to get back the units for flux
    med_spec_units = median_spec*med_new

    return med_spec_units

# see UVJ_plots.py for code to stack up UVJ diagram
"""
fig, (ax1, ax2) = plt.subplots(2, figsize=(15,7))
grid = gridspec.GridSpec(2, 1, height_ratios=[6,1])
ax1 = plt.subplot(grid[0])
#ax2 = plt.subplot(grid[1])
ax1.plot(new_wavs, med_spec_units*10**18, color="black", lw=1.5)#, label = '1.25 $\mathrm{\AA}\ $ Sampling')
#ax2.plot(new_wavs, errs, color="red", lw=1.5)
ax1.set_ylabel("Flux $(10^{-18}\ \mathrm{W/m^2/\\AA)}$", size=15) # $(\mathrm{W/m^2/\\AA)}$
#ax2.set_ylabel("Flux Error", size=15)
ax1.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
ax1.set_xlim(2300, 4300)
#ax2.set_xlim(2300, 4300)
ax1.set_title('Median Stacked Spectra of galaxies 1 < z < 1.5')
plt.savefig('stack_plot_all.pdf')
"""
"""
new_waves = np.arange(2400, 4200, 1.5)
objects1 =  np.array(passive_cut_new['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut_new['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut_new['CAT'].str.decode("utf-8"))

med_stack = stacks(new_waves, objects1)

plt.figure(figsize=(15,7))
plt.plot(new_waves, med_stack, color="black", lw=1.5 )
plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
plt.xlim(2300, 4250)
plt.title('Median Stacked Spectra of galaxies 1 < z < 1.5', size=17)
plt.savefig('stack_plot_all.pdf')
"""
