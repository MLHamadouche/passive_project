import numpy as np
import bagpipes as pipes
import LoadData as ld
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 4})
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


Mg_UV = {}
Mg_UV["name"] = "Mg_UV"
Mg_UV["type"] = "break"
#Mg_UV["feature"] = [2625., 2725.]
Mg_UV["continuum"] = [[2525., 2625.], [2725.,2825.]]

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
index_list2 = [Mg_UV,C29_33]

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
inds_Mg = np.zeros((len(ID),2))

d4000 = np.zeros((len(ID),2))

#print(inds_Mg)
colour_ind = np.zeros((len(ID), 2))
for i in range(len(ID)):

    galaxy = pipes.galaxy(ID[i], load_data = ld.load_vandels_spectra, index_list=index_list2,
    spectrum_exists = True, photometry_exists=False, index_redshift = redshift[i], load_indices = 'from_spectrum')

    try:
        galaxy2 = pipes.galaxy(ID[i], load_data = ld.load_vandels_spectra, index_list=index_list,
        spectrum_exists = True, photometry_exists=False, index_redshift = redshift[i], load_indices = 'from_spectrum')
    except:
         continue
#load_indices=load_hd_ew_value,

    #print(galaxy.indices[i])
    #inds_Mggalaxy.indices[i].values
    #print(inds_Mg[i,:])
    d4000[i] = galaxy2.indices[1]
    inds[i] = galaxy2.indices[0]
    inds_Mg[i] = galaxy.indices[0]
    colour_ind[i] = galaxy.indices[1]

print(inds[:,0])
print(d4000[:,0])
#print(len(inds))


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
fig1, ax1 = plt.subplots(figsize=[10,6.5])
ax1.errorbar(d4000[:,0],inds[:,0],xerr = d4000[:,1], yerr = inds[:,1], ecolor='k',marker='s', ms = 5, color ='k', linestyle =" ", elinewidth = 0.2, capsize = 4, mew = 0.2)
#cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'Redshift, z', size=12)
#ax1.scatter(colour_index, Mg_UV_index)
ax1.set_ylabel("EW($H\delta$)", size=13)
ax1.set_xlabel("$D_{n}4000$", size=13)
#ax1.set_xlim(0.98, 1.9)
#ax1.set_ylim(-5., 12.)
plt.title('$D_{n}4000 \ v \ EW(H\delta)$', size =14)# excluding possible AGN (CDFS + UDS)')
#plt.savefig('EW(Hdelta)vDn4000_bagpipes.pdf')
plt.show()

"""
fig1, ax1 = plt.subplots(figsize=[10,6.5])

#im1 = ax1.scatter(inds_Mg, colour_ind, s=130, c=redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
ax1.errorbar(colour_ind[:,0],inds_Mg[:,0],marker='s', xerr = colour_ind[:,1], yerr = inds_Mg[:,1], ecolor='k', color ='pink', linestyle =" ")
#cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'Redshift, z', size=12)
#ax1.scatter(colour_index, Mg_UV_index)
ax1.set_ylabel("Mg_UV", size=13)
ax1.set_xlabel("C(29-33)", size=13)
plt.title('colour index v Mg_UV', size =14)# excluding possible AGN (CDFS + UDS)')
#plt.savefig('EW(Hdelta)vDn4000.pdf')
plt.show()

"""
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
