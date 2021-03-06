import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import spectres
import LoadData as ld
from glob import glob
import os
from astropy.stats import sigma_clip
from indices import C29_33, Dn4000, Mg_UV, H_delta
#from stacks import stacks
plt.rc('text', usetex=True)

#uds_xray = Table.read('FirstProjectCatalogs/uds_possible_xray_matches_massi.fits').to_pandas()
#cdfs_xray = Table.read('FirstProjectCatalogs/cdfs_possible_xray_matches_massi.fits').to_pandas()
both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()
df1 = pd.DataFrame(both_xray)#, index = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8")))
#df1 = pd.DataFrame(uds_xray)
#ID_uds = np.array(uds_xray['FIELD'].str.decode("utf-8").str.rstrip() + uds_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + "SELECT")
ID_list1 = df1.set_index(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))
#ID_cdfs = np.array(cdfs_xray['FIELD'].str.decode("utf-8").str.rstrip() + cdfs_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + "SELECT")
#z_uds = uds_xray['zspec'].values
#z_cdfs = cdfs_xray['zspec'].values
objects1  = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))
#redshifts = np.concatenate([z_uds, z_cdfs])
#redshifts = both_xray['zspec']
#print(redshifts)
#ID_uds = uds_xray['ID_1']
#IDs_uds = df1.set_index(s.decode('utf-8') for s in uds_xray['ID_1'])
#all_IDs = np.concatenate([ID_uds, ID_cdfs])


#objects1 =  np.array(passive_cut_new['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut_new['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut_new['CAT'].str.decode("utf-8"))

passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
df2 = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )
ID_list2 = df2.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))

#redshifts = passive_cut['zspec']

all_obs = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))
#print(ID_list.type)
objects_list = list(set(all_obs).difference(objects1))

#df = df2.merge(df1, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only']
#print(df)
#ID_list = df.set_index(df['FIELD'].str.decode("utf-8").str.rstrip() + df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + df['CAT'].str.decode("utf-8"))
#new_df= df1.index.difference(df2.index)

#new_df = df2.merge(df1, indicator = True, how='left').loc[lambda x : x['_merge']!='both'].drop_duplicates(keep=False)
new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)
#print(new_df)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))
#print(ID_list.index.values)
#print(passive_cut)

vandels_cat_new = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()
agn = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))

ID_pipes = vandels_cat_new['#ID'].values
df_bp = pd.DataFrame(vandels_cat_new)
IDs = df_bp.set_index(s.decode('utf-8') for s in vandels_cat_new['#ID'])
print(len(IDs))


new_wavs = np.arange(2400, 4200, 1.25)

def stacks(objects_list): #input array of redshifts for given list of objects
    new_spec = np.zeros(len(new_wavs))
    new_errs = np.zeros(len(new_wavs))

    old_spec =[]
    old_errs = []

    spec = []
    spec_err = []
    med_norm = []
    med_spec_units  = []

    colour_index = []
    D4000_index =[]
    Mg_UV_index = []
    H_delta_EW = []
    all_ages = []
    ages = []
    masses = []
    median_spec = np.zeros(len(new_wavs))
    errs = np.zeros(len(new_wavs))
    spec_ = np.zeros(len(new_wavs))
    spec_errs = np.zeros(len(new_wavs))
    med_spec_units=[]
    ID_plz = []
    med_spectrum =np.zeros(len(new_wavs))
    for ID in objects_list:
        ID_plz.append(ID)
        z = ID_list.loc[ID, 'zspec']
        ages_all = IDs.loc[ID, "mass_weighted_age_50"]
        ages.append(ages_all)
        masses.append(IDs.loc[ID, 'stellar_mass_50'])
        #define the list of objects outside of function- must be an array of IDs as strings
        #in format e.g. CDFS_HST034930SELECT
        spectrum = ld.load_vandels_spectra(ID)
        wav_mask = (spectrum[:,0]>5200) & (spectrum[:,0]<9250)

        flux = spectrum[:,1][wav_mask]
        flux_errs = spectrum[:,2][wav_mask]
        wavs = spectrum[:,0][wav_mask]

        zeros_mask = (flux == 0.)|(flux_errs == 0.)
        flux[zeros_mask] = np.nan
        flux_errs[zeros_mask] = np.nan

        rest_wavs = wavs/(1.0 + z)
        mask =  (rest_wavs > 3000) & (rest_wavs < 3500) # fairly featureless region of spectrum
        old_spec = flux/np.nanmedian(flux[mask]) #normalisation median from that region
        old_errs = flux_errs/np.nanmedian(flux[mask])

        med_norm.append(np.nanmedian(flux[mask]))
        new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)
        colour_index.append(C29_33(new_spec))
        D4000_index.append(Dn4000(new_spec))
        Mg_UV_index.append(Mg_UV(new_spec))
        H_delta_EW.append(H_delta(new_spec, new_wavs))
        spec.append(new_spec)
        spec_err.append(new_errs)

    spec = np.transpose(spec)

    spec_err = np.transpose(spec_err)
    standev_err = []
    med_new = np.nanmedian(med_norm)
    for m in range(len(new_wavs)):
        spec_ = spec[m,:]
        #print(spec_.shape)
        spec_errs = spec_err[m,:]
        #standev_err[m] = np.std(spec_, axis=0)
        median_spec[m]=np.nanmedian(spec_)

    med_spec_units = median_spec*med_new

    return med_spec_units, colour_index, D4000_index, Mg_UV_index, H_delta_EW, ages, masses, ID_plz

#med_stack = stacks(objects1)
#print(med_stack)
"""
med_spec_units, colour_index, D4000_index, Mg_UV_index, H_delta_EW, ages, masses, ID_plz = stacks(ID_list.index.values)

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
print(len(D4000_index))

fig1, ax1 = plt.subplots(figsize=[12,8.5])
im1 = ax1.scatter(colour_index, Mg_UV_index, s=180, c=masses, cmap=plt.cm.magma, marker='o', linewidth=0.5 )#edgecolors='black',
#ax1.errorbar(colour_index, Mg_UV_index ,marker='s', xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
cbar = fig1.colorbar(im1, ax=ax1)
cbar.set_label(r'log$_{10}$(M*/M$_{\odot}$)', size=12)
#ax1.scatter(colour_index, Mg_UV_index ,marker='o',color = 'k',zorder = 1)
ax1.set_ylim(0.8, 2.3)
ax1.set_xlim(0.2, 1.2)
ax1.set_xlabel("C(29-33)", size=13)
ax1.set_ylabel("Mg$_{UV}$", size=13)
#plt.title('', size =14)# excluding possible AGN (CDFS + UDS)')
ax1.set_title(f'Full sample ({len(D4000_index)} objects) coloured by mass')
plt.savefig('all_obs228_mgUVvC29_33_cbar_masses.pdf')
plt.close()

fig2, ax2 = plt.subplots(figsize=[12,8.5])
im2 = ax2.scatter(D4000_index, H_delta_EW, s=180, c=masses, cmap=plt.cm.magma, marker='o', linewidth=0.5 )#edgecolors='black',
#ax1.errorbar(colour_index, Mg_UV_index ,marker='s', xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
cbar2 = fig2.colorbar(im2, ax=ax2)
cbar2.set_label(r'log$_{10}$(M*/M$_{\odot}$)', size=12)
#ax2.scatter(D4000_index, H_delta_EW, marker='o',color = 'k',zorder = 1)
ax2.set_xlabel("D$_{n}4000$", size=13)
ax2.set_ylabel("EW(H$\delta$)", size=13)
ax2.set_title(f'Full sample ({len(D4000_index)} objects) coloured by mass')
plt.savefig('all_obs228_d4000vHdelta_cbar_mass.pdf')
plt.close()
"""


def plot_stackssingle(stack):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, stack*10**18, color="black", lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    #plt.xlim(2300, 4250)
    plt.ylim(0 ,2.0)
    plt.title('Median Stacked Spectra excluding possible AGN (CDFS + UDS)')# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stack_plot_exc_x_ray_check_nans.pdf')
    plt.close()

#med_stack_missing = stacks(missing, ID_list)
#med_stack_missing = stacks(objects_list)
#plot_stackssingle(med_stack_missing)


"""
plt.figure(figsize=(15,7))
plt.plot(new_wavs, errs*10**18, color="black", lw=1.5 )
plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
#plt.xlim(2300, 4250)
plt.ylim(0 ,2.0)
plt.title('Error Spectrum (standard deviation of median stacks)')# excluding possible AGN (CDFS + UDS)')
plt.savefig('errs_stack_plot_exc_x_ray_check_nans.pdf')
plt.close()

"""
"""

#fig = plt.figure(figsize = (15,7))
def subplots(med_stack_with_agn, med_stack_missing):
    fig, (ax1, ax2) = plt.subplots(2, figsize = (15,7),sharex=True, sharey=True)
    fig.suptitle('Median Stacked Spectra')
    ax1.plot(new_wavs, med_stack_missing*10**18, 'k', lw=0.9)
    ax2.plot(new_wavs, med_stack_with_agn*10**18, 'r',lw=0.9)
    ax1.set_title('Excluding AGN (N=228)')
    ax2.set_title('Only AGN (N=36)')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)",)
    ax2.set_ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", )
    ax1.set_ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", )
    plt.savefig('stack_plot_two_figs.pdf')

#subplots(med_stack_agn, med_stack_missing)

def plot_stacks(stack1, stack2):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, stack1*10**18, color="red", lw=1.5, ls ='-', label = 'possible AGN (N=36)')
    plt.plot(new_wavs, stack2*10**18, color="black", lw=1.6, label = 'excluding possible AGN (N=228)')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    #plt.xlim(2300, 4250)
    plt.ylim(0 ,2.0)
    plt.legend()
    plt.title('Median Stacked Spectra')# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stack_plot_on_top.pdf')
    plt.close()


"""
#print(len(med_stack))
#print(med_stack)
#plot_stackssingle(med_stack_with_agn)
#plot_stacks(med_stack_agn, med_stack_missing)

#plot_stacks(med_stack_missing, med_stack_with_agn)
