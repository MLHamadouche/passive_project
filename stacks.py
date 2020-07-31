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

globpath = os.path.join('vandelsspec/', '*.fits')
filelist = glob(globpath)
hdulist = fits.open(filelist[-1])
delt_wav = hdulist[0].header['CDELT1']
wav_first_pixel = hdulist[0].header['CRVAL1']

passive_cut = Table.read('FirstProjectCatalogs/xmatch_spec_derived237objs.fits').to_pandas()
df = pd.DataFrame(passive_cut)
ID_list = df.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8"))

redshifts = passive_cut['z_spec']
objects = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8"))

new_wavs = np.arange(2400, 4200, 1.5) #new wavelength grid

new_spec = np.zeros(len(new_wavs))
new_errs = np.zeros(len(new_wavs))
rest_wavs = []
old_spec =[]
old_errs = []
spectrum = []
#resampling spectra onto new wavelength grid with spectres
z_mask = (redshifts <= 1.5) & (redshifts>= 1.0) # 89 objects
new_objs = objects[z_mask]
specc = np.zeros((len(new_wavs),len(new_objs)))
specerrs = np.zeros((len(new_wavs),len(new_objs)))
spec = []
spec_err = []
for ID in new_objs:
    z = ID_list.loc[ID, 'z_spec']
    spectrum=ld.load_vandels_spectra(ID)
    rest_wavs = spectrum[:,0]/(1.+z)
    mask = (rest_wavs < 3500) & (rest_wavs > 3000)
    old_spec = spectrum[:,1]/np.median(spectrum[:,1][mask]) #normalisation ????
    old_errs = spectrum[:,2]/np.median(spectrum[:,2][mask])

    new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)
    spec.append(new_spec)
    spec_err.append(new_errs)

spec = np.transpose(spec)
spec_err = np.transpose(spec_err)

#for n in range(len(new_objs)):
#    specc[:,n] = spec
#    specerrs[:,n] = spec_err

#print(specc, specerrs)
#print(specc.shape)
#print(specc[0:10])

median_spec = np.zeros(len(new_wavs))
errs = np.zeros(len(new_wavs))
spec_ = np.zeros(len(new_wavs))
spec_errs = np.zeros(len(new_wavs))

for m in range(len(new_wavs)):
    spec_ = spec[m,:]
    spec_errs = spec_err[m,:]
    median_spec[m] = np.median(spec_)
    errs[m] = np.median(spec_errs) #nopes
    #errs[m] = np.sqrt(1/np.sum(spec_errs**2)) #??? weighted errors


fig, (ax1, ax2) = plt.subplots(2, figsize=(12,7))
grid = gridspec.GridSpec(2, 1, height_ratios=[4,1])
ax1 = plt.subplot(grid[0])
ax2 = plt.subplot(grid[1])
ax1.plot(new_wavs, median_spec, color="black", lw=1.5)#, label = '1.25 $\mathrm{\AA}\ $ Sampling')
ax2.plot(new_wavs, errs, color="red", lw=1.5)
ax1.set_ylabel("Flux (arbitrary units)", size=15) # $(\mathrm{W/m^2/\\AA)}$
ax2.set_ylabel("Flux Error", size=15)
ax2.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
ax1.set_xlim(2300, 4300)
ax2.set_xlim(2300, 4300)
ax1.set_title('Median Stacked Spectra of galaxies 1 < z < 1.5')
plt.savefig('stack_test_plot_all.pdf')
