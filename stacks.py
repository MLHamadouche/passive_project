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
#print(ID_list)

new_wavs = np.arange(2500, 4200,5) #new wavelength grid

new_spec = np.zeros(len(new_wavs))
new_errs = np.zeros(len(new_wavs))
rest_wavs = []
old_spec =[]
old_errs = []
spectrum = []


#resampling spectra onto new wavelength grid with spectres

for ID in objects:
    z = ID_list.loc[ID, 'z_spec']
    spectrum=ld.load_vandels_spectra(ID)
    rest_wavs = spectrum[:,0]/(1.+z)
    mask = (rest_wavs < 3500) & (rest_wavs > 3000)
    old_spec = spectrum[:,1]/np.max(spectrum[:,1][mask]) #normalisation ????
    old_errs = spectrum[:,2]/np.max(spectrum[:,1][mask])

    new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)

#plt.plot(new_wavs, new_spec)
#plt.show()
specc = np.zeros((len(new_wavs),len(objects)))
specerrs = np.zeros((len(new_wavs),len(objects)))

for n in range(len(objects)):
    specc[:,n] = new_spec
    specerrs[:,n] = new_errs

#print(specc, specerrs)
print(specc.shape)

median_spec = np.zeros(len(new_wavs))
errs = np.zeros(len(new_wavs))
spec_ = np.zeros(len(new_wavs))
spec_errs = np.zeros(len(new_wavs))

for m in range(len(new_wavs)):
    spec_ = specc[m,:]
    spec_errs = specerrs[m,:]
    median_spec[m] = np.nanmedian(spec_)
    errs[m] = np.sqrt(1/np.sum(spec_errs**2)) #???

print(median_spec.shape)
print(new_wavs.shape)

fig, (ax1, ax2) = plt.subplots(2, figsize=(12,7))
grid = gridspec.GridSpec(2, 1, height_ratios=[3,1])
ax1 = plt.subplot(grid[0])
ax2 = plt.subplot(grid[1])
ax1.plot(new_wavs, median_spec, color="black", lw=1.5, )
ax2.plot(new_wavs, errs, color="red", lw=1.5)
ax1.set_ylabel("Flux $( \mathrm{W/m^2/\\AA)}$", size=15)
ax2.set_ylabel("Flux Error", size=15)
ax2.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
ax1.set_xlim(2400, 4300)
ax2.set_xlim(2400, 4300)
plt.show()
