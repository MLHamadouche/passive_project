import numpy as np
from astropy.io import fits
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt
from glob import glob
import os
from astropy import units as u
from astropy.coordinates import SkyCoord

globpath = os.path.join('vandelsspec/', '*.fits')
print(globpath)
# glob searches through directories similar to the Unix shell
filelist = glob(globpath)
objID = []
objRA = []
objDEC = []
redshift = []
flag = []

hdulist = fits.open(filelist[-1])
print(hdulist.info())
#print(hdulist[0].columns)
#flux = hdulust[0].data
flux = hdulist[0].data
flux_err = hdulist[3].data
print(flux, flux_err)
#print(np.arange(0,30000))
redshift = hdulist[0].header['HIERARCH PND Z']

wav_first_pixel = hdulist[0].header['CRVAL1']
delt_wav = hdulist[0].header['CDELT1']

#wave = np.arange(wav_first_pixel, , delt_wav )
#waves = np.array(wave)
#print(wave)

wa_end = wav_first_pixel + (2154*delt_wav)
wa = np.arange(wav_first_pixel, wa_end, delt_wav)
print(len(wa))
spectrum  = np.c_[wa, flux, flux_err]

#wave=np.arange(0,30000,30000/2154)
#waves = np.array(wave)
#print(sp1[0].header)
plt.plot(wa*(1+float(redshift)), flux, linewidth = 0.5)
ax = plt.gca()
# recompute the ax.dataLim
ax.relim()
# update ax.viewLim using the new dataLim
ax.autoscale_view()
plt.draw()
plt.show()


"""
for i in range(len(filelist)):
    sp = fits.open(filelist[i])
    header = sp[0].header
    objID.append(header['HIERARCH PND OBJID'])
    objRA.append(header['HIERARCH PND OBJRA'])
    objDEC.append(header['HIERARCH PND OBJDEC'])
    redshift.append(header['HIERARCH PND Z'])
    flag.append(header['HIERARCH PND ZFLAGS'])

#print(header['HIERARCH PND OBJID'])

ra = np.array(objRA)
dec = np.array(objDEC)
col1 = fits.Column(name='ID', format='10A', array=objID)
col2 = fits.Column(name='z_spec', format='E', array=redshift)
col3 = fits.Column(name='RA', format='E',  array=objRA)
col4 = fits.Column(name='DEC', format='E', array=objDEC)
col5 = fits.Column(name='ZFLAGS', format='E', array=flag)
"""
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
#file =  "all_vandels_spec.fits"
#hdu.writeto(file)


"""
derived_cat = Table.read('catalogs/VANDELS_DERIVED_v1.1.fits').to_pandas()
ra2 = np.array(derived_cat['RA'])
dec2 = np.array(derived_cat['DEC'])
c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)

#print(idx)
#hdulist = fits.open("vandels_spectra_june20/sc_CDFS002669_P1M3Q4_P1M4Q4_022_1.fits")

matched_cat = Table.read('all_vandels_spec.fits').to_pandas()
df = pd.DataFrame(matched_cat)
mask1 = (matched_cat['z_spec']>= 1.0) & (matched_cat['z_spec'] <= 1.5)
mask2 = (matched_cat['ZFLAGS'] >= 3.0) & (matched_cat['ZFLAGS'] <=4.0)
#print(hdulist.info)
mask = mask1 & mask2
filtered_Data = df[mask]
fits = Table.from_pandas(filtered_Data)
table = fits.write('passive_cut_to_xmatch.fits')

matched_cat = Table.read('matched_spec_derived.fits').to_pandas()
df = pd.DataFrame(matched_cat)

mask1 = (matched_cat['z_spec'] >= 1.0) & (matched_cat['z_spec'] <= 1.5)
mask2 = (matched_cat['ZFLAGS'] >= 3.0) & (matched_cat['ZFLAGS'] <=4.0)

mask = mask1 & mask2

filtered_data = df[mask]
print(filtered_data)
fits  = Table.from_pandas(filtered_data)
file = 'z1&1.5filtered_matched_spec.fits'
tab = fits.write(file)
"""

#df = pd.DataFrame(df.values[flag3], df.index[mask], df.columns).astype(df.dtypes)
#flags = matched_cat['ZFLAGS'] == 3. and  matched_cat['ZFLAGS'] == 4.

#print(flags)
