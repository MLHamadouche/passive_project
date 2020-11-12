import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import spectres
import LoadData as ld
import matplotlib.pyplot as plt
from indices import C29_33, Dn4000, Mg_UV, H_delta
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
    colour_index = []
    D4000_index =[]
    Mg_UV_index = []
    H_delta_EW = []
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
        #spectres resamples spectra to the new wavelengths
        colour_index.append(C29_33(new_spec))
        D4000_index.append(Dn4000(new_spec))
        Mg_UV_index.append(Mg_UV(new_spec))
        H_delta_EW.append(H_delta(new_spec, new_wavs))
        np.savetxt('hdelta_d4000.dat', np.c_[H_delta_EW, D4000_index], delimiter=' ')
        np.savetxt('MgUV_colour.dat', np.c_[Mg_UV_index, colour_index], delimiter=' ')
        #print(len(H_delta_EW))
        spec.append(new_spec)
        spec_err.append(new_errs)



    spec = np.transpose(spec)
    #colour_index.append(C29_33(spec))
    #D4000_index.append(Dn4000(spec))
    #Mg_UV_index.append(Mg_UV(spec))
    from astropy.io import ascii
    table = ascii.read('indices_for_stacks.dat',  data_start=0).to_pandas()
    df = pd.DataFrame(table)
    print(df)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig1, ax1 = plt.subplots(figsize=[12,8.5])
    #im1 = ax1.scatter(D4000_index, H_delta_EW, s=130, c=new_redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
    #ax1.errorbar(colour_index, Mg_UV_index ,marker='s', """xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
    #cbar = fig1.colorbar(im1, ax=ax1)
    ax1.scatter(colour_index, Mg_UV_index ,marker='o',color = 'k',zorder = 1)
    #cbar.set_label(r'Redshift, z', size=12)
    ax1.scatter([df['C(29-33)'][0],df['C(29-33)'][2],df['C(29-33)'][4]],[df['Mg_UV'][0],df['Mg_UV'][2],df['Mg_UV'][4]], marker = '*', color = 'r', s=350, zorder = 2, edgecolors = 'k', label = 'above relation')
    ax1.scatter([df['C(29-33)'][1],df['C(29-33)'][3],df['C(29-33)'][5]],[df['Mg_UV'][1],df['Mg_UV'][3],df['Mg_UV'][5]], marker = '*', color = 'g', s=350, zorder = 2, edgecolors = 'k', label = 'below relation')
    ax1.set_ylim(0.8, 2.3)
    ax1.set_xlim(0.2, 1.2)
    #ax1.scatter(colour_index, Mg_UV_index)

    ax1.set_xlabel("C(29-33)", size=13)
    ax1.set_ylabel("MgUV", size=13)
    #plt.title('', size =14)# excluding possible AGN (CDFS + UDS)')
    plt.legend()
    plt.savefig('massi_colurvmgUV_blackcircles_plusstacks.pdf')
    plt.close()
    from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    fig2, ax2 = plt.subplots(figsize=[12,8.5])
    ax2.scatter(D4000_index, H_delta_EW, marker='o',color = 'k',zorder = 1)
    #cbar.set_label(r'Redshift, z', size=12)
    f1 = ax2.scatter(df['Dn4000'][0],df['EW(Hdelta)'][0], marker = '*', color = 'r', s=300, zorder = 2, edgecolors = 'k',) #label = '10.5 - 10.75')
    f2 = ax2.scatter(df['Dn4000'][2], df['EW(Hdelta)'][2], marker = 'd', color = 'r', s=300, zorder = 2, edgecolors = 'k',) #label = '10.75 - 11.0')
    f3 = ax2.scatter(df['Dn4000'][4], df['EW(Hdelta)'][4], marker = '^', color = 'r', s=300, zorder = 2, edgecolors = 'k',) #label = '11.0 - 11.3')
    s1 = ax2.scatter(df['Dn4000'][1],df['EW(Hdelta)'][1], marker = '*', color = 'g', s=300, zorder = 2, edgecolors = 'k', label = '10.5 - 10.75')
    s2 = ax2.scatter(df['Dn4000'][3], df['EW(Hdelta)'][3], marker = 'd', color = 'g', s=300, zorder = 2, edgecolors = 'k', label = '10.75 - 11.0')
    s3 = ax2.scatter(df['Dn4000'][5], df['EW(Hdelta)'][5], marker = '^', color = 'g', s=300, zorder = 2, edgecolors = 'k', label = '11.0 - 11.3')
    #legend1 = ax2.legend(*scatter.legend_elements(),
    #                loc="lower left", title="Above or below Van Der Wel Relation")
    ax2.set_xlabel("Dn4000)", size=13)
    ax2.set_ylabel("EW(Hdelta)", size=13)
    first_legend = plt.legend([s1, s2, s3], ['10.5 - 10.75', '10.75 - 11.0','11.0 - 11.3'], loc = 'lower right')
    ax2 = plt.gca().add_artist(first_legend)
    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color='red', label='above relation')
    green_patch = mpatches.Patch(color='green', label='below relation')
    plt.legend(handles=[red_patch, green_patch], loc = 'upper right')
    #ax2.legend(['red patch', 'green patch'], ['above relation', 'below relation'],
    #      loc='upper right', frameon=True)
    #from matplotlib.legend import Legend
    #leg = Legend(ax2, [(s1,s2,s3)], [('10.5 - 10.75', '10.75 - 11.0', '11.0 - 11.3')], loc='lower right')
    #ax2.add_artist(leg)
    plt.savefig('massi_hdeltad4000_blackcircles_plusstacks.pdf')
    plt.close()

    spec_err = np.transpose(spec_err)
    standev_err = []
    med_new = np.nanmedian(med_norm)
    for m in range(len(new_wavs)):
        spec_ = spec[m,:]
        spec_errs = spec_err[m,:]
        #standev_err[m] = np.std(spec_, axis=0)
        median_spec[m]=np.nanmedian(spec_)

    med_spec_units = median_spec*med_new #test removing normalisations to see fluxes nh

    return med_spec_units #returns an array of the new median stacked fluxes



####example:#####
#med_stack_missing = stacks(objects_list)
