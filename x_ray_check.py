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
#from stacks import stacks


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
missing = list(set(all_obs).difference(objects1))

#df = df2.merge(df1, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only']
#print(df)
#ID_list = df.set_index(df['FIELD'].str.decode("utf-8").str.rstrip() + df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + df['CAT'].str.decode("utf-8"))
#new_df= df1.index.difference(df2.index)

#new_df = df2.merge(df1, indicator = True, how='left').loc[lambda x : x['_merge']!='both'].drop_duplicates(keep=False)
new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)
#print(new_df)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))
#print(ID_list)
#print(len(ID_list))
print(len(missing))

new_wavs = np.arange(2400, 4200, 1.25)

def stacks(list_of_IDs, pd_list):
    new_spec = np.zeros(len(new_wavs))
    new_errs = np.zeros(len(new_wavs))

    old_spec =[]
    old_errs = []

    #resampling spectra onto new wavelength grid with spectres
    #z_mask = (redshifts <= 1.5) & (redshifts>= 0.9) # 89 objects
    #new_objs = objects[z_mask] #uncomment if want mask
    #print('lens', len(new_wavs), len(new_objs))
    #specc = np.zeros((len(new_wavs),len(new_objs)))
    #specerrs = np.zeros((len(new_wavs),len(new_objs)))
    spec = []
    spec_err = []
    med_norm = []
    med_spec_units  = []

    median_spec = np.zeros(len(new_wavs))
    errs = np.zeros(len(new_wavs))
    spec_ = np.zeros(len(new_wavs))
    spec_errs = np.zeros(len(new_wavs))
    med_spec_units=[]
    med_spectrum =[]
    for ID in list_of_IDs:
        z = pd_list.loc[ID, 'zspec']
        spectrum = ld.load_vandels_spectra(ID)
        wav = spectrum[:,0]
        flux = spectrum[:,1]
        errors = spectrum[:,0]
        #print(ID)
        #plt.plot(wav, flux)
        #plt.savefig(str(ID)+'.pdf')
        #plt.close()
        for f in flux:
            if f==0:
                f = np.nan
        for e in errors:
            if e==0:
                e = np.nan
        plt.plot(wav, flux)
        plt.savefig(str(ID)+'.pdf')
        plt.close()

        rest_wavs = spectrum[:,0]/(1.0 + z)
        #print(f'rest_wavs={rest_wavs}')
        #print(max(rest_wavs), min(rest_wavs))
        #print(f'flux={spectrum[:,1]}')
        mask =  (rest_wavs > 3000) & (rest_wavs < 3500) # fairly featureless region of spectrum
        old_spec = spectrum[:,1]/np.nanmedian(spectrum[:,1][mask]) #normalisation median from that region
        old_errs = spectrum[:,2]/np.nanmedian(spectrum[:,2][mask])
        med_norm.append(np.nanmedian(spectrum[:,1][mask]))
        #print(f'old_spec:\n {old_spec}')
        #input()
        new_spec, new_errs = spectres.spectres(new_wavs, rest_wavs, old_spec, spec_errs=old_errs)
        #print(new_spec)
        for i in new_spec:
            i  = float(i)
        for j in new_errs:
            j = float(j)
        spec.append(new_spec)
        spec_err.append(new_errs)

    spec = np.transpose(spec)
    #print(spec)
    spec_err = np.transpose(spec_err)

    med_new = np.nanmedian(med_norm)
    for m in range(len(new_wavs)):
        spec_ = spec[m,:]
        #print(spec_)
        #input()
        spec_errs = spec_err[m,:]
        median_spectrum = median_spec[m] #newline
        #print(np.median(spec_))
        #input()
        median_spec[m]=np.nanmedian(spec_)

        median_spectrum = median_spec[m]
        #median_spectrum_units = median_spectrum * np.median(spec_)

    med_spec_units = median_spec*med_new

    return med_spec_units

#med_stack = stacks(objects1)
#print(med_stack)

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

#new_waves = np.arange(2400, 4300, 1.25)
#med_stack_agn = stacks(objects1, ID_list1)
med_stack_missing = stacks(missing, ID_list)
plot_stackssingle(med_stack_missing)
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



#print(len(med_stack))
#print(med_stack)
#plot_stackssingle(med_stack_with_agn)
#plot_stacks(med_stack_agn, med_stack_missing)

#plot_stacks(med_stack_missing, med_stack_with_agn)
