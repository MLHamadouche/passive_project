import numpy as np
import bagpipes as pipes
#import LoadData as ld
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import os
import re
from collections import OrderedDict
from glob import glob

# Define Ha equivalent width index (these numbers are just an example)
ha_ew = {}
ha_ew["name"] = "Ha_ew"
ha_ew["type"] = "EW"
ha_ew["feature"] = [6553., 6573.] # The spectral region being measured
ha_ew["continuum"] = [[6533., 6553.], [6573., 6593.]] # The side regions to be compared against

hd_ew = {}
hd_ew["name"] = "Hd_ew"
hd_ew["type"] = "EW"
hd_ew["feature"] = [4082., 4122.]
hd_ew["continuum"] = [[4042., 4082.], [4122,4162]]


Mg_UV_red = {}
Mg_UV_red["name"] = "Mg_UV"
Mg_UV_red["type"] = "break"
#Mg_UV["feature"] = [2625., 2725.]
Mg_UV_red["continuum"] = [[2525., 2625.], [2725.,2825.]] #red

Mg_UV_blue = {}
Mg_UV_blue["name"] = "Mg_UV"
Mg_UV_blue["type"] = "break"
#Mg_UV["feature"] = [2625., 2725.]
Mg_UV_blue["continuum"] = [[2625., 2725.], [2725.,2825.]] #blue

D4000 = {}
D4000["name"] = "D4000"
D4000["type"] = "break"
D4000["continuum"] = [[3850.,3950.],[4000., 4100.]]

C29_33 = {}
C29_33["name"] = "C29_33"
C29_33["type"] = "break"
C29_33["continuum"] = [[2700.,3100.],[3100., 3500.]]

# Define index list variable, the idea here is similar to filt_list
index_list = [hd_ew, D4000]
index_list2 = [Mg_UV_red,Mg_UV_blue]
index_list3 = [C29_33]
# Set up simple model_components dictionary
nebular = {}
nebular["logU"] = -3.

const = {}
const["age_min"] = 0.
const["age_max"] = 1.
const["metallicity"] = 1.
const["massformed"] = 10.

model_comp = {}
model_comp["redshift"] = 1.5660
model_comp["constant"] = const
model_comp["nebular"] = nebular


### Example of making a model_galaxy to demonstrate index calculation
model = pipes.model_galaxy(model_comp, index_list=index_list)

print(model.indices) # EW specified above (absorption is positive)
cat = Table.read("Re_cat.fits").to_pandas()
df = pd.DataFrame(cat)
df = df.groupby(df['log10(M*/Msun)']>10.4).get_group(True)
R_e = df["Re_kpc"]
redshift = df['redshifts'].values
#print(redshifts)
ids = np.array(df['IDs'].str.decode("utf-8").str.rstrip().values)
print(len(ids))
ID = []
for i in ids:
    ID.append(i)
### Example of loading in index data to a galaxy object

# Define load_indices function to return equivalent width value
def load_ha_ew_value(ID):
    indices = np.zeros((1, 2))
    indices[0, 0] = -15. #index value
    indices[0, 1] = 0.1 #index uncertainty
    return indices

# Example load_data function for photometry
#def load_phot_data(ID):
#    photometry = np.random.randn(12, 2)
#    return photometry
inds = np.zeros((len(ID),2))
inds_Mg_red = np.zeros((len(ID),2))
inds_Mg_blue = np.zeros((len(ID),2))

d4000 = np.zeros((len(ID),2))


colour_ind = np.zeros((len(ID), 2))

def load_vandels_spectra(ID_no):
    pre = ID_no.split('-')[0]
    #print(pre)
    new_ID = re.search('\d+', ID_no).group()
    #ID = ID.lstrip('0')
    ID_new = str(pre) + str(new_ID)
    globpath = os.path.join('new_vandels_spec/', '*.fits')
    filelist = glob(globpath)
    #print(filelist[0])
    for i in range(len(filelist)):
        if ID_new in str(filelist[i]):
            hdulist = fits.open(filelist[i])
            flux = hdulist[0].data
            flux_err = hdulist[3].data
            redshift = hdulist[0].header['HIERARCH PND Z']
            wav_first_pixel = hdulist[0].header['CRVAL1']
            delt_wav = hdulist[0].header['CDELT1']
            wa_end = wav_first_pixel + (2154*delt_wav)
            wave = np.arange(wav_first_pixel, wa_end, delt_wav)
            wav_mask = (wave>5200) & (wave<9250)
            spectrum=np.c_[wave[wav_mask], flux[wav_mask], flux_err[wav_mask]]

    return spectrum
def load_vandels_spectra_freq(ID_no):
    pre = ID_no.split('-')[0]
    #print(pre)
    new_ID = re.search('\d+', ID_no).group()
    #ID = ID.lstrip('0')
    ID_new = str(pre) + str(new_ID)
    globpath = os.path.join('new_vandels_spec/', '*.fits')
    filelist = glob(globpath)
    #print(filelist[0])
    for i in range(len(filelist)):
        if ID_new in str(filelist[i]):
            hdulist = fits.open(filelist[i])
            flux = hdulist[0].data
            flux_err = hdulist[3].data
            redshift = hdulist[0].header['HIERARCH PND Z']
            wav_first_pixel = hdulist[0].header['CRVAL1']
            delt_wav = hdulist[0].header['CDELT1']
            wa_end = wav_first_pixel + (2154*delt_wav)
            wave = np.arange(wav_first_pixel, wa_end, delt_wav)
            wav_mask = (wave>5200) & (wave<9250)
            spectrum2=np.c_[wave[wav_mask], (wave[wav_mask]**2/3*10**8)*flux[wav_mask], (wave[wav_mask]**2/3*10**8)*flux_err[wav_mask]]

    return spectrum2

for i in range(len(ID)):
    spectrum = load_vandels_spectra(ID[i])
    galaxy = pipes.galaxy(ID[i], load_data = load_vandels_spectra, index_list=index_list2,
    spectrum_exists = True, photometry_exists=False, index_redshift = redshift[i], load_indices = 'from_spectrum')

    galaxy3 = pipes.galaxy(ID[i], load_data = load_vandels_spectra_freq, index_list=index_list3,
    spectrum_exists = True, photometry_exists=False, index_redshift = redshift[i], load_indices = 'from_spectrum')

    try:
        galaxy2 = pipes.galaxy(ID[i], load_data = load_vandels_spectra, index_list=index_list,
        spectrum_exists = True, photometry_exists=False, index_redshift = redshift[i], load_indices = 'from_spectrum')
    except:
         continue
    colour_ind[i] = galaxy3.indices
    wavs = spectrum[:,0][i]
    d4000[i] = galaxy2.indices[1]
    inds[i] = galaxy2.indices[0]
    inds_Mg_red[i] = galaxy.indices[0]
    inds_Mg_blue[i] = galaxy.indices[1]


#print(inds[:,0])
#print(d4000[:,0])
#print(len(inds))
av_MgUV = (inds_Mg_red + inds_Mg_blue)/2

#print(av_MgUV)
#wavs = np.arange(2400, 4200, 1.25)
colour_ind_a = 2.5*np.log10((colour_ind[:,0]))
colour_ind[:,1] = colour_ind[:,1]#/colour_ind[:,0]
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('title')
plt.rc('xlabel',)
plt.rc('ylabel', )
fig1, ax1 = plt.subplots(figsize=[10,6.5])
ax1.errorbar(d4000[:,0],inds[:,0],xerr = d4000[:,1], yerr = inds[:,1], ecolor='k',marker='s', ms = 5, color ='k', linestyle =" ", elinewidth = 0.2, capsize = 4, mew = 0.2)
#cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'Redshift, z', size=12)
#ax1.scatter(colour_index, Mg_UV_index)
ax1.set_ylabel("EW($H\delta$)", size=13)
ax1.set_xlabel("$D_{n}4000$", size=13)
ax1.set_xlim(0.98, 1.9)
#ax1.set_ylim(-5., 12.)
plt.title('$D_{n}4000 \ v \ EW(H\delta)$', size =14)# excluding possible AGN (CDFS + UDS)')
#plt.savefig('EW(Hdelta)vDn4000_bagpipes.pdf')
plt.close()

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('title',)
plt.rc('xlabel',)
plt.rc('ylabel', )
fig1, ax1 = plt.subplots(figsize=[10,6.5])
#im1 = ax1.scatter(inds_Mg, colour_ind, s=130, c=redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
ax1.errorbar(colour_ind_a, av_MgUV[:,0],marker='s', xerr = 0.434*2.5*colour_ind[:,1]/colour_ind[:,0], yerr = av_MgUV[:,1], ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3)
#cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'Redshift, z', size=12)
#ax1.scatter(colour_index, Mg_UV_index)
ax1.set_ylabel("$\mathrm{Mg_{UV}}$", size=13)
ax1.set_xlabel("C(29-33)", size=13)
#ax1.set_yticks([1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
#ax1.set_xticks([0.4, 0.6, 0.8, 1.0])
ax1.set_ylim(0.6, 2.3)
ax1.set_xlim(0.0, 1.2)
plt.title('$\mathrm{C(29 - 33) \ v \ Mg_{UV}}$', size =14)# excluding possible AGN (CDFS + UDS)')
plt.savefig('colourVmgUV_bagpipes.pdf')
plt.close()



mg_uv_massi, colour_massi = np.loadtxt('MgUV_colour.dat', delimiter = ' ', unpack=True)
hdelta_massi, d4000_massi = np.loadtxt('hdelta_d4000.dat',  delimiter = ' ', unpack=True)

one = np.linspace(-2,7, 100)
plt.scatter(d4000_massi, d4000[:,0])
plt.plot(one, one)
plt.xlim(1.1, 1.8)
plt.ylim(1.1,1.8)
plt.xlabel('d4000_massi')
plt.ylabel('d4000_bagpipes')
plt.show()

plt.scatter(hdelta_massi, inds[:,0])
plt.plot(one, one)
plt.xlim(-1., 6.5)
plt.ylim(-1.,6.5)
plt.xlabel('hdelta_massi')
plt.ylabel('hdelta_bagpipes')
plt.show()

plt.scatter(mg_uv_massi, av_MgUV[:,0])
plt.plot(one, one)
plt.xlim(0.8, 2.2)
plt.ylim(0.8,2.2)
plt.xlabel('mg_UV_massi')
plt.ylabel('mg_UV_bagpipes')
plt.show()


plt.scatter(colour_massi, colour_ind_a)
plt.plot(one, one)
plt.xlim(0.0, 1.5)
plt.ylim(0.0,1.2)
plt.xlabel('c(29-33)_massi')
plt.ylabel('c(29-33)_bagpipes')
plt.show()

"""
# Set up example filt_list
filt_list = np.loadtxt("filters/list_of_filters.txt", dtype="str")


# Create example galaxy object to demonstrate loading index values
galaxy = pipes.galaxy("test", load_phot_data, filt_list=filt_list,
                      load_indices=load_ha_ew_value, index_list=index_list,
                      spectrum_exists=False)

print galaxy.indices # Similar idea to galaxy.photometry


### Example of fitting indices

# Set up simple fit_instructions dictionary
nebular = {}
nebular["logU"] = -3.

const = {}
const["age_min"] = (0., 0.1)
const["age_max"] = 1.
const["metallicity"] = (0.1, 2.5)
const["massformed"] = (0., 13.)

fit_info = {}
fit_info["redshift"] = 0.001
fit_info["constant"] = const
fit_info["nebular"] = nebular


# Set up example fit object
fit = pipes.fit(galaxy, fit_info)


# Fit the photometry + index value
fit.fit(verbose=True)
"""
