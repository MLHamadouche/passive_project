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
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import scipy
from stacking_code import stacks
plt.rc('text', usetex=True)
both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()
df1 = pd.DataFrame(both_xray)#, index = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8")))
passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
df2 = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )
new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)
ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))
#candels = Table.read('FirstProjectCatalogs/concat_candels_passive_match.fits').to_pandas()
concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df = pd.DataFrame(concat_3dhst)
ID_list1 = np.array(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
ID_ = df.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
agn = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))

vandels_cat_pipes = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()
df1 = pd.DataFrame(vandels_cat_pipes)
ID_pipes = vandels_cat_pipes['#ID'].values
IDs = df1.set_index(s.decode('utf-8') for s in vandels_cat_pipes['#ID'])

#df2 = pd.DataFrame(vandels_cat_pipes)
#df2 = df2.groupby(df2['stellar_mass_50']>10.4).get_group(True)
#ssfr2 = df2["ssfr_50"]
UV = vandels_cat_pipes['UV_colour_50']
VJ = vandels_cat_pipes['VJ_colour_50']
age1 = vandels_cat_pipes['mass_weighted_age_50']
ssfr = vandels_cat_pipes['ssfr_50']
re = concat_3dhst['re']
redshifts = concat_3dhst['zspec']

mass = passive_cut["log10(M*)"]

#print(len(ssfr))

#input()
list_IDs = []
ages = []
masses = []
re = []
re_errs =[]

cosmo  = FlatLambdaCDM(H0=70, Om0=0.3)
Mpc_to_kpc = 1000
for i in ID_list1:
    for j in ID_pipes:
        j = j.decode("utf-8")
        if i == j:
            list_IDs.append(i)
new_IDs = []
new_redshifts = []
q_ratio = []
ssfr2 = []
for ID in list_IDs:
    if ID not in agn:
        new_redshifts.append(ID_.loc[ID, "zspec"])
        ssfr2.append(IDs.loc[ID, "ssfr_50"])
        ages.append(IDs.loc[ID, "mass_weighted_age_50"])
        #print(all_ages)
        masses.append(ID_list.loc[ID, "log10(M*)"])
        #print(all_masses)
        re.append(ID_.loc[ID, 're'])
        re_errs.append(ID_.loc[ID, 'dre'])
        q_ratio.append(ID_.loc[ID, 'q'])
        new_IDs.append(ID)
#print(all_sizes)
#new IDs has passive galaxies excluding possible AGN found in the 3DHST catalog/CANDELS catalog
masses = np.array(masses)
arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(np.array(new_redshifts))
print('max=', max(re))

Re_kpc = (np.array(re)*u.arcsec)/arcsec_per_kpc
Re_kpc_errs = (np.array(re_errs)*u.arcsec)/arcsec_per_kpc


Rc = (np.sqrt(np.array(q_ratio))*Re_kpc) /u.kpc
Rc_errs = (np.sqrt(np.array(q_ratio))*Re_kpc_errs) /u.kpc

col1 = fits.Column(name='IDs', format='30A', array=new_IDs)
col2 = fits.Column(name='redshifts', format='E', array=new_redshifts)
col3 = fits.Column(name='age', format='E',  array=ages)
col4 = fits.Column(name='log10(M*/Msun)', format='E', array=masses)
col5 = fits.Column(name='Re_kpc', format='E', array=Re_kpc/u.kpc)
col6 = fits.Column(name='Re_kpc_errs', format='E', array=Re_kpc_errs/u.kpc)
col7 = fits.Column(name ='SSFR', format = 'E', array = ssfr2)
hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
#file =  "Re_cat.fits"
#hdu.writeto(file)

#print(f'Rc : {Rc}')

x = np.linspace(9.5, 11.5, len(Rc))
shen = 0.56*x + (-5.54)

"""
plt.plot(x, shen, 'k', lw = 0.7)
plt.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 8)
plt.ylabel(r'$\mathrm{log_{10}{(R_{c})}}$', size = 8)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('3D-HST log10(Rc) versus log10(M*/Msun)', size = 8)
plt.xlim(9.75, 11.5)
#plt.show()
plt.close()
#plt.savefig('Rc_mass_relation_3DHST_density.pdf')
"""
x_model = np.linspace(10.4, 11.2, len(Rc))

c_model = np.arange(-8.0, -3.0, 0.01)
model_mask = (masses > 10.4)
best_chi  = np.inf

errs = 0.434*(Rc_errs/Rc)

for cvals in range(len(c_model)):
    c_vals = c_model[cvals]

    y_model = 0.56*masses + c_vals
    diffs = y_model - np.log10(Rc)
    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq = np.sum((diffs**2)/((10*errs)**2))

    if chisq < best_chi:
        best_c = c_vals
        best_chi = chisq

#print(f'best_c {best_c} \n best_chi {best_chi}')

cat2 = Table.read("Re_cat.fits").to_pandas()
df2 = pd.DataFrame(cat2)
df2 = df2.groupby(df2['log10(M*/Msun)']>10.4).get_group(True)
R_e = df2["Re_kpc"]
R_e_errs = df2["Re_kpc_errs"]
mass = df2["log10(M*/Msun)"]
ssfr_50 = df2["SSFR"]
#print(len(df2))
redshifts = df2['redshifts']
alpha = 0.76
log_A = 0.22
x = np.linspace(9.5, 11.5, len(R_e))
log_Reff = log_A + alpha*np.log10((10**x)/(5*10**10))


def vdw_relation(logA, alpha, x_values):
    logR_eff = logA + alpha*(np.log10((10**x_values)/(5*10**10)))
    return logR_eff

#print(shen[20])
#print(f'std mass = {np.std(masses)}')
#print(f'std Re = {np.std(Re_kpc_errs/u.kpc)}')
print(10**0.22)
log_A_model = np.arange(-0.3, 0.3, 0.01)
best_chi2 = np.inf

for cvals2 in range(len(log_A_model)):
    c_vals2 =log_A_model[cvals2]

    vdw_model = vdw_relation(c_vals2, 0.76, mass)
    diffs2 = vdw_model - np.log10(R_e)
    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq2 = np.sum((diffs2**2)/((4.34*(R_e_errs/R_e))))

    if chisq2 < best_chi2:
        best_chi_vdw = chisq2
        best_c_vdw = c_vals2

#print(f'best_c_vdw: {best_c_vdw} \n best_chi_vdw: {best_chi_vdw}')
#print(vdw_model)

vdw_norm_model = vdw_relation(best_c_vdw, 0.76, x)

y_model = 0.56*x + (best_c)

#print(best_c_vdw - 0.22)
#print((best_c_vdw - 0.22))
#print(y_model[20])

#print(y_model - shen)
#diff_shen_me = abs(y_model[0]) - abs(shen[0])
#print(10**diff_shen_me)



fig1, ax1 = plt.subplots(figsize=[12,8.5])
im1 = ax1.scatter(mass, np.log10(R_e), s=130, c=redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
cbar1 = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar1.set_label(r'Redshift (z)', size=12)
#ax1.plot(x, y_model, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
ax1.plot(x, vdw_norm_model, 'k', lw = 1.2, label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
ax1.plot(x, log_Reff, 'r' ,lw = 1.2, label= 'v. der Wel z = 1.25 ETG relation (not normalised)')
#ax1.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax1.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(prop={'size': 10})
plt.title('3D-HST log10(Re) versus log10(M*/Msun)', size = 13)
plt.xlim(10.3, 11.4)
#plt.show()
#plt.savefig('Re_v_M*_test_with_Arjens_relationandbestfit_TEST.pdf')
plt.close()


mg_uv_massi, colour_massi = np.loadtxt('MgUV_colour.dat', delimiter = ' ', unpack=True)
hdelta_massi, d4000_massi = np.loadtxt('hdelta_d4000.dat',  delimiter = ' ', unpack=True)


print(len(d4000_massi))
print(np.median(np.log10(R_e)))
fig2, ax2 = plt.subplots(figsize=[12,8.5])
im2 = ax2.scatter(mass, np.log10(R_e), s=130, c=hdelta_massi, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
cbar2 = fig2.colorbar(im2, ax=ax2)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar2.set_label(r'$EW(H\delta)$', size=12)
#ax1.plot(x, y_model, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
ax2.plot(x, vdw_norm_model, 'k', lw = 1.2, label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
#ax1.plot(x, log_Reff, 'r' ,lw = 1.2, label= 'v. der Wel z = 1.25 ETG relation (not normalised)')
#ax1.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
ax2.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax2.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(prop={'size': 10})
plt.title('3D-HST log10(Re) versus log10(M*/Msun)', size = 13)
plt.xlim(10.3, 11.4)
#plt.show()
#plt.savefig('Re_v_M*_cbar_Hdelta.pdf')
plt.close()

fig, (ax, ax3) = plt.subplots(1, 2, figsize=(14,6))
im = ax.scatter(mass, np.log10(R_e), s=130, c=d4000_massi, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'$D_{n}4000$', size=12)
ax.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)
ax.set_xlim(10.3, 11.3)
ax.set_ylim(-0.4, 1.1)
im3 = ax3.scatter(mass, np.log10(R_e), s=130, c=hdelta_massi, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.set_label(r'$EW(H\delta)$', size=12)
ax3.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)# label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
ax3.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax3.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax3.set_xlim(10.3, 11.3)
ax3.set_ylim(-0.4, 1.1)
#plt.savefig("Re_v_M_cbar_both.pdf")
plt.close()

figure, (axs, axs1) = plt.subplots(1, 2, figsize=(14,6))
ims = axs.scatter(mass, np.log10(R_e), s=130, c=redshifts, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbars = figure.colorbar(ims, ax=axs)
cbars.set_label(r'Redshift ($z$)', size=12)
axs.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
axs.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
axs.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)
axs.set_xlim(10.3, 11.3)
axs.set_ylim(-0.4, 1.1)
ims1 = axs1.scatter(mass, np.log10(R_e), s=130, c=ssfr_50, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )
cbars1 = fig.colorbar(ims1, ax=axs1)
cbars1.set_label(r'sSFR', size=12)
axs1.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)# label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
axs1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
axs1.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
axs1.set_xlim(10.3, 11.3)
axs1.set_ylim(-0.4, 1.1)
#plt.savefig("Re_v_M_cbar_redshiftandssfr.pdf")
plt.close()



fig3, (ax4, ax5) = plt.subplots(1, 2, figsize=(14,6))
im4 = ax4.scatter((np.log10(R_e)- np.median(np.log10(R_e))), hdelta_massi, s=130, c=d4000_massi, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbar4 = fig.colorbar(im4, ax=ax4)
cbar4.set_label(r'$D_{n}4000$', size=12)
#ax.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax4.set_xlabel(r'$\Delta \mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax4.set_ylabel(r'$EW(H\delta)$', size = 12)
#ax4.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)
#ax4.set_xlim(10.3, 11.3)
#ax4.set_ylim(-0.4, 1.1)
im5 = ax5.scatter((np.log10(R_e)- np.median(np.log10(R_e))), d4000_massi, s=130, c=hdelta_massi, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbar5 = fig3.colorbar(im5, ax=ax5)
cbar5.set_label(r'$EW(H\delta)$', size=12)
#ax5.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)# label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
#ax5.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax5.set_xlabel(r'$\Delta \mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax5.set_ylabel(r'$D_{n}4000$', size = 12)
#ax5.set_xlim(10.3, 11.3)
#ax5.set_ylim(-0.4, 1.1)
#plt.savefig("deltaRe_v_d4000andhdelta.pdf")
plt.close()


#cat = Table.read('TEST_CAT.fits').to_pandas()
cat = Table.read("Re_cat.fits").to_pandas()
df = pd.DataFrame(cat)
df = df.groupby(df['log10(M*/Msun)']>10.4).get_group(True)
R_e = df["Re_kpc"]
#print(np.log10(R_c))
#print(y_model)
index = np.log10(R_e).index.to_list()
#print(index)
IDs = cat['IDs']
IDs_above = []
IDs_below = []
mass_ = df['log10(M*/Msun)']
#size = np.array(np.log10(R_c))
size  = np.array(np.log10(R_e))
#y_model = np.array(y_model)
vdw_norm_model = np.array(vdw_norm_model)

#mask = (size>y_model)
#mask1 = (size<y_model)
mask = (size>vdw_norm_model)
mask1 = (size<vdw_norm_model)
#index_masked = np.log10(R_c)[mask].index.to_list()
#index_masked2 = np.log10(R_c)[mask1].index.to_list()
index_masked = np.log10(R_e)[mask].index.to_list()
index_masked2 = np.log10(R_e)[mask1].index.to_list()
#print(size.shape)
#print('size_mask1:',size[mask1])
#print(len(size[mask1]))

#print(index_masked)
IDs_above = IDs[index_masked].str.decode("utf-8").str.rstrip().values
IDs_below= IDs[index_masked2].str.decode("utf-8").str.rstrip().values

#print(IDs_above, IDs_below)
all_IDs = np.concatenate((IDs_above, IDs_below), axis = None)

print(len(all_IDs))

""""
for names in all_IDs:
    spectrum = ld.load_vandels_spectra(names)
    wav_mask = (spectrum[:,0]>5200) & (spectrum[:,0]<9250)
    flux = spectrum[:,1]#[wav_mask]

    colour_index.append(C29_33(flux))
    D4000_index.append(Dn4000(flux))
    Mg_UV_index.append(Mg_UV(flux))

plt.figure(figsize=(12,8.5))
plt.scatter(colour_index, Mg_UV_index, color=color, lw=1.5 )
plt.xlabel("C(29-33)", size=17)
plt.ylabel("Mg_UV index", size=17)
#plt.xlim(2350, 4250)
#plt.ylim(0 ,1.75)
plt.title('colour index versus Mg_UV index', size =18)# excluding possible AGN (CDFS + UDS)')
plt.savefig('C_2933vMgUV_allobjects.pdf')
plt.close()
"""

new_wavs = np.arange(2400, 4200, 1.25)
#print((0.56*x + (best_c)).shape)
#print(mask)
#for i in size[mask]:

#    IDs_above.append(df.loc[(np.log10(R_c)[mask]==i),'IDs'].values.tolist())

#print(len(IDs_above), len(IDs_below))
from x_ray_check import stacks

stacking_all ,colour_index, D4000_index, Mg_UV_index, H_delta_EW, ages, masses = stacks(all_IDs)
print(len(all_IDs))
input()

data_indices = {'C(29-33)': colour_index, 'Mg_UV': Mg_UV_index, 'H_delta_EW': H_delta_EW, 'D4000': D4000_index, 'ages': ages, 'masses':masses}
df_new = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
D4000_ind= df_new['D4000']
H_delta_EW_ind = df_new['H_delta_EW']
ages = df_new['ages']
mass = df_new['masses']
colour_ind = df_new['C(29-33)']
mgUV_ind = df_new['Mg_UV']


#data_indices = {'C(29-33)': colour_index, 'Mg_UV': Mg_UV_index, 'H_delta_EW': H_delta_EW, 'D4000': D4000_index, 'ages': ages, 'masses':masses}
df_new1 = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_new1 = df_new1.groupby((df_new1['masses']>10.5) & (df_new1['masses']<=10.75)).get_group(True)

D4000_ind1= df_new1['D4000']
H_delta_EW_ind1 = df_new1['H_delta_EW']
ages1 = df_new1['ages']
mass1 = df_new1['masses']
colour_ind1 = df_new1['C(29-33)']
mgUV_ind1 = df_new1['Mg_UV']

df_new2 = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_new2 = df_new2.groupby((df_new2['masses']>10.75) & (df_new2['masses']<=11.0)).get_group(True)
D4000_ind2= df_new2['D4000']
H_delta_EW_ind2 = df_new2['H_delta_EW']
ages2 = df_new2['ages']
mass2 = df_new2['masses']
colour_ind2 = df_new2['C(29-33)']
mgUV_ind2 = df_new2['Mg_UV']

df_new3 = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_new3 = df_new3.groupby((df_new3['masses']>11.0) & (df_new3['masses']<=11.3)).get_group(True)
D4000_ind3= df_new3['D4000']
H_delta_EW_ind3 = df_new3['H_delta_EW']
ages3 = df_new3['ages']
mass3 = df_new3['masses']
colour_ind3 = df_new3['C(29-33)']
mgUV_ind3 = df_new3['Mg_UV']


av_mg1, av_d40001, av_hdelta1, av_colour1 = np.nanmedian(mgUV_ind1), np.nanmedian(D4000_ind1), np.nanmedian(H_delta_EW_ind1),np.nanmedian(colour_ind1)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')
av_mg2, av_d40002, av_hdelta2, av_colour2 = np.nanmedian(mgUV_ind2), np.nanmedian(D4000_ind2), np.nanmedian(H_delta_EW_ind2),np.nanmedian(colour_ind2)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')

av_mg3, av_d40003, av_hdelta3, av_colour3 = np.nanmedian(mgUV_ind3), np.nanmedian(D4000_ind3), np.nanmedian(H_delta_EW_ind3),np.nanmedian(colour_ind3)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')


from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
figtree, axtree = plt.subplots(figsize=[12,8.5])
imtree = axtree.scatter(D4000_ind, H_delta_EW_ind, s=200, c=ages, cmap=plt.cm.magma, marker='o',linewidth=0.5)
#print(len(D4000_ind_all))
#ax2.scatter(D4000_index, H_delta_EW, marker='o',color = 'k',zorder = 1)
cbartree = figtree.colorbar(imtree, ax=axtree)
#cbartree.set_label(r'log$_{10}$(M*/M$_{\odot}$)', size=12)
cbartree.set_label(r'Age (Gyr)', size=12)
#im1 = axtree.scatter(av_d40001, av_hdelta1, s=200, c=ages[ages1], cmap=plt.cm.magma, marker='*',linewidth=0.5, edgecolors = 'k', label = '10.5 - 10.75')
#im2 = axtree.scatter(av_d40002, av_hdelta2, s=200, c=ages[ages2], cmap=plt.cm.magma, marker='d',linewidth=0.5, edgecolors = 'k', label = '10.75 - 11.0')
#im3 = axtree.scatter(av_d40003, av_hdelta3, s=200, c=ages[ages3], cmap=plt.cm.magma, marker='^',linewidth=0.5, edgecolors = 'k', label = '11.0 - 11.3')
s1 = axtree.scatter(av_d40001, av_hdelta1, marker = '*', color = 'teal', s=300, zorder = 2, edgecolors = 'k', label = '10.5 - 10.75')
s2 = axtree.scatter(av_d40001, av_hdelta2, marker = 'd', color = 'teal', s=250, zorder = 2, edgecolors = 'k', label = '10.75 - 11.0')
s3 = axtree.scatter(av_d40002, av_hdelta3, marker = '^', color = 'teal', s=250, zorder = 2, edgecolors = 'k', label = '11.0 - 11.3')
axtree.set_xlabel("D$_{n}$4000", size=13)
axtree.set_ylabel("EW(H$\delta$)", size=13)
plt.legend()
plt.savefig('NO_SIZE_corrected_hdeltad4000cbar_medianmassbinned_agecbar_only78.pdf')
plt.close()

input()

stacking_above, colour_index_above, D4000_index_above, Mg_UV_index_above, H_delta_EW_above, ages_above, masses_above = stacks(IDs_above)
data_indices_above = {'C(29-33)': colour_index_above, 'Mg_UV': Mg_UV_index_above, 'H_delta_EW': H_delta_EW_above, 'D4000': D4000_index_above, 'ages': ages_above, 'masses':masses_above}
df_over = pd.DataFrame(data_indices_above, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_over = df_over.groupby(df_over['masses']<=10.5).get_group(True)
D4000_above = df_over['D4000']
H_delta_EW_above = df_over['H_delta_EW']
ages_above = df_over['ages']
masses_over = df_over['masses']
colour_index_o = df_over['C(29-33)']
mgUV_over = df_over['Mg_UV']

av_mg_1, av_d4000_1, av_hdelta_1, av_colour_1 = np.nanmedian(mgUV_over), np.nanmedian(D4000_above), np.nanmedian(H_delta_EW_above),np.nanmedian(colour_index_o)
print(f'For 10.4 <=log10(M*) <=10.5: \n MgUV = {av_mg_1} \n D4000 = {av_d4000_1} \n EW(Hdelta) = {av_hdelta_1} \n C29-33 = {av_colour_1}')


df_2 = pd.DataFrame(data_indices_above, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_2 = df_2.groupby((df_2['masses']>10.5) & (df_2['masses']<=10.75)).get_group(True)
D4000_2 = df_2['D4000']
H_delta_EW_2 = df_2['H_delta_EW']
ages_2 = df_2['ages']
masses_2 = df_2['masses']
colour_index_2 = df_2['C(29-33)']
mgUV_2 = df_2['Mg_UV']
av_mg_2, av_d4000_2, av_hdelta_2, av_colour_2 = np.nanmedian(mgUV_2), np.nanmedian(D4000_2), np.nanmedian(H_delta_EW_2),np.nanmedian(colour_index_2)
print(f'For 10.5 <=log10(M*) <=10.75: \n MgUV = {av_mg_2} \n D4000 = {av_d4000_2} \n EW(Hdelta) = {av_hdelta_2} \n C29-33 = {av_colour_2}')


df_3 = pd.DataFrame(data_indices_above, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_3 = df_3.groupby((df_3['masses']>10.75) & (df_3['masses']<=11.0)).get_group(True)
D4000_3 = df_3['D4000']
H_delta_EW_3 = df_3['H_delta_EW']
ages_3 = df_3['ages']
masses_3= df_3['masses']
colour_index_3 = df_3['C(29-33)']
mgUV_3 = df_3['Mg_UV']

av_mg_3, av_d4000_3, av_hdelta_3, av_colour_3 = np.nanmedian(mgUV_3), np.nanmedian(D4000_3), np.nanmedian(H_delta_EW_3),np.nanmedian(colour_index_3)
print(f'For 10.75 <log10(M*) <=11.0: \n MgUV = {av_mg_3} \n D4000 = {av_d4000_3} \n EW(Hdelta) = {av_hdelta_3} \n C29-33 = {av_colour_3}')



df_4 = pd.DataFrame(data_indices_above, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_4 = df_4.groupby((df_4['masses']>11.0) & (df_4['masses']<=11.3)).get_group(True)

D4000_4 = df_4['D4000']
H_delta_EW_4 = df_4['H_delta_EW']
ages_4 = df_4['ages']
masses_4= df_4['masses']
colour_index_4 = df_4['C(29-33)']
mgUV_4 = df_4['Mg_UV']

av_mg_4, av_d4000_4, av_hdelta_4, av_colour_4 = np.nanmedian(mgUV_4), np.nanmedian(D4000_4), np.nanmedian(H_delta_EW_4),np.nanmedian(colour_index_4)
print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_4} \n D4000 = {av_d4000_4} \n EW(Hdelta) = {av_hdelta_4} \n C29-33 = {av_colour_4}')

print(len(df_over), len(df_2), len(df_3), len(df_4))

stacking_below, colour_index_below, D4000_index_below, Mg_UV_index_below, H_delta_EW_below, ages_below, masses_below = stacks(IDs_below)

data_indices_below = {'C(29-33)': colour_index_below, 'Mg_UV': Mg_UV_index_below, 'H_delta_EW': H_delta_EW_below, 'D4000': D4000_index_below, 'ages': ages_below, 'masses':masses_below}
df_bel = pd.DataFrame(data_indices_below, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
D4000_ind_bel= df_bel['D4000']
H_delta_EW_ind_bel = df_bel['H_delta_EW']
ages_bel = df_bel['ages']
mass_bel = df_bel['masses']
colour_ind_bel = df_bel['C(29-33)']
mgUV_ind_bel = df_bel['Mg_UV']
"""
df_over_bel = pd.DataFrame(data_indices_below, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_over_bel = df_over_bel.groupby(df_over_bel['masses']<=10.5).get_group(True)
D4000_above_bel = df_over_bel['D4000']
H_delta_EW_above_bel = df_over_bel['H_delta_EW']
ages_above_bel = df_over_bel['ages']
masses_over_bel = df_over_bel['masses']
colour_index_o_bel = df_over_bel['C(29-33)']
mgUV_over_bel = df_over_bel['Mg_UV']
"""
#av_mg_1_bel, av_d4000_1_bel, av_hdelta_1_bel, av_colour_1_bel = np.nanmedian(mgUV_over_bel), np.nanmedian(D4000_above_bel), np.nanmedian(H_delta_EW_above_bel),np.nanmedian(colour_index_o_bel)
#print(f'For 10.4 <=log10(M*) <=10.5: \n MgUV = {av_mg_1_bel} \n D4000 = {av_d4000_1_bel} \n EW(Hdelta) = {av_hdelta_1_bel} \n C29-33 = {av_colour_1_bel}')


df_2_bel = pd.DataFrame(data_indices_below, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_2_bel = df_2_bel.groupby((df_2_bel['masses']>10.5) & (df_2_bel['masses']<=10.75)).get_group(True)
D4000_2_bel = df_2_bel['D4000']
H_delta_EW_2_bel = df_2_bel['H_delta_EW']
ages_2_bel = df_2_bel['ages']
masses_2_bel = df_2_bel['masses']
colour_index_2_bel = df_2_bel['C(29-33)']
mgUV_2_bel = df_2_bel['Mg_UV']
av_mg_2_bel, av_d4000_2_bel, av_hdelta_2_bel, av_colour_2_bel = np.nanmedian(mgUV_2_bel), np.nanmedian(D4000_2_bel), np.nanmedian(H_delta_EW_2_bel),np.nanmedian(colour_index_2_bel)
print(f'For 10.5 <=log10(M*) <=10.75: \n MgUV = {av_mg_2_bel} \n D4000 = {av_d4000_2_bel} \n EW(Hdelta) = {av_hdelta_2_bel} \n C29-33 = {av_colour_2_bel}')


df_3_bel = pd.DataFrame(data_indices_below, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_3_bel = df_3_bel.groupby((df_3_bel['masses']>10.75) & (df_3_bel['masses']<=11.0)).get_group(True)
D4000_3_bel = df_3_bel['D4000']
H_delta_EW_3_bel = df_3_bel['H_delta_EW']
ages_3_bel = df_3_bel['ages']
masses_3_bel= df_3_bel['masses']
colour_index_3_bel = df_3_bel['C(29-33)']
mgUV_3_bel = df_3_bel['Mg_UV']

av_mg_3_bel, av_d4000_3_bel, av_hdelta_3_bel, av_colour_3_bel = np.nanmedian(mgUV_3_bel), np.nanmedian(D4000_3_bel), np.nanmedian(H_delta_EW_3_bel),np.nanmedian(colour_index_3_bel)
print(f'For 10.75 <log10(M*) <=11.0: \n MgUV = {av_mg_3_bel} \n D4000 = {av_d4000_3_bel} \n EW(Hdelta) = {av_hdelta_3_bel} \n C29-33 = {av_colour_3_bel}')



df_4_bel = pd.DataFrame(data_indices_below, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_4_bel = df_4_bel.groupby((df_4_bel['masses']>11.0) & (df_4_bel['masses']<=11.3)).get_group(True)

D4000_4_bel = df_4_bel['D4000']
H_delta_EW_4_bel = df_4_bel['H_delta_EW']
ages_4_bel = df_4_bel['ages']
masses_4_bel= df_4_bel['masses']
colour_index_4_bel = df_4_bel['C(29-33)']
mgUV_4_bel = df_4_bel['Mg_UV']

av_mg_4_bel, av_d4000_4_bel, av_hdelta_4_bel, av_colour_4_bel = np.nanmedian(mgUV_4_bel), np.nanmedian(D4000_4_bel), np.nanmedian(H_delta_EW_4_bel),np.nanmedian(colour_index_4_bel)
print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_4_bel} \n D4000 = {av_d4000_4_bel} \n EW(Hdelta) = {av_hdelta_4_bel} \n C29-33 = {av_colour_4_bel}')



plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
"""
fig1, ax1 = plt.subplots(figsize=[12,8.5])
#im1 = ax1.scatter(D4000_index, H_delta_EW, s=130, c=new_redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
#ax1.errorbar(colour_index, Mg_UV_index ,marker='s',xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
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
"""

both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()

dfx = pd.DataFrame(both_xray)#, index = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8")))

passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()

dfp = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )

new_df=pd.concat([dfx,dfp]).drop_duplicates(subset = 'ID_1', keep=False)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))

IDs_vals = ID_list.index.values
stack_all_all, colour_index_all, D4000_index_all, Mg_UV_index_all, H_delta_EW_all, ages_all, masses_all = stacks(IDs_vals)
print(len(IDs_vals))

data_indices_all = {'C(29-33)': colour_index_all, 'Mg_UV': Mg_UV_index_all, 'H_delta_EW': H_delta_EW_all, 'D4000': D4000_index_all, 'ages': ages_all, 'masses':masses_all}
df_xp = pd.DataFrame(data_indices_all, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])

D4000_ind_all= df_xp['D4000']
H_delta_EW_ind_all = df_xp['H_delta_EW']
ages_all = df_xp['ages']
mass_all = df_xp['masses']
colour_ind_all = df_xp['C(29-33)']
mgUV_ind_all = df_xp['Mg_UV']
print(len(D4000_ind_all))

df_xp = df_xp.groupby((df_xp['masses']>10.4) & (df_xp['masses']<=10.5)).get_group(True)

D4000_ind_all1= df_xp['D4000']
H_delta_EW_ind_all1 = df_xp['H_delta_EW']
ages_all1 = df_xp['ages']
mass_all1 = df_xp['masses']
colour_ind_all1 = df_xp['C(29-33)']
mgUV_ind_all1 = df_xp['Mg_UV']
print(len(D4000_ind_all1))


df_xp3 = pd.DataFrame(data_indices_all, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_xp3 = df_xp3.groupby((df_xp3['masses']>10.5) & (df_xp3['masses']<=10.75)).get_group(True)
D4000_ind_all3= df_xp3['D4000']
H_delta_EW_ind_all3 = df_xp['H_delta_EW']
ages_all3 = df_xp3['ages']
mass_all3 = df_xp3['masses']
colour_ind_all3 = df_xp3['C(29-33)']
mgUV_ind_all3 = df_xp3['Mg_UV']
print(len(D4000_ind_all3))

df_xp4 = pd.DataFrame(data_indices_all, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_xp4 = df_xp4.groupby((df_xp4['masses']>10.75) & (df_xp4['masses']<=11.0)).get_group(True)
D4000_ind_all4= df_xp4['D4000']
H_delta_EW_ind_all4 = df_xp4['H_delta_EW']
ages_all4 = df_xp4['ages']
mass_all4 = df_xp4['masses']
colour_ind_all4 = df_xp4['C(29-33)']
mgUV_ind_all4 = df_xp4['Mg_UV']
print(len(D4000_ind_all4))

df_xp5 = pd.DataFrame(data_indices_all, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses'])
df_xp5 = df_xp5.groupby((df_xp5['masses']>11.0) & (df_xp5['masses']<=11.3)).get_group(True)
D4000_ind_all5= df_xp5['D4000']
H_delta_EW_ind_all5 = df_xp5['H_delta_EW']
ages_all5 = df_xp5['ages']
mass_all5 = df_xp5['masses']
colour_ind_all5 = df_xp5['C(29-33)']
mgUV_ind_all5 = df_xp5['Mg_UV']
print(len(D4000_ind_all5))

av_mg_all_1, av_d4000_all_1, av_hdelta_all_1, av_colour_all_1 = np.nanmedian(mgUV_ind_all1), np.nanmedian(D4000_ind_all1), np.nanmedian(H_delta_EW_ind_all1),np.nanmedian(colour_ind_all1)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')
av_mg_all_2, av_d4000_all_2, av_hdelta_all_2, av_colour_all_2 = np.nanmedian(mgUV_ind_all3), np.nanmedian(D4000_ind_all3), np.nanmedian(H_delta_EW_ind_all3),np.nanmedian(colour_ind_all3)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')

av_mg_all_3, av_d4000_all_3, av_hdelta_all_3, av_colour_all_3 = np.nanmedian(mgUV_ind_all4), np.nanmedian(D4000_ind_all4), np.nanmedian(H_delta_EW_ind_all4),np.nanmedian(colour_ind_all4)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')

av_mg_all_4, av_d4000_all_4, av_hdelta_all_4, av_colour_all_4 = np.nanmedian(mgUV_ind_all5), np.nanmedian(D4000_ind_all5), np.nanmedian(H_delta_EW_ind_all5),np.nanmedian(colour_ind_all5)
#print(f'For 11.0 <log10(M*) <=11.3: \n MgUV = {av_mg_all_1} \n D4000 = {av_d4000_all_1} \n EW(Hdelta) = {av_hdelta_all_1} \n C29-33 = {av_colour_all_1}')


from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
figtree2, axtree2 = plt.subplots(figsize=[12,8.5])
imtree2 = axtree2.scatter(D4000_ind, H_delta_EW_ind, s=200, c=ages, cmap=plt.cm.magma, marker='o',linewidth=0.5)
#imtree2 = axtree2.scatter(D4000_ind_all, H_delta_EW_ind_all, s=200, c=ages_all, cmap=plt.cm.magma, marker='o',linewidth=0.5)
#print(len(D4000_ind_all))
#ax2.scatter(D4000_index, H_delta_EW, marker='o',color = 'k',zorder = 1)
cbartree = figtree2.colorbar(imtree2, ax=axtree2)
#cbartree.set_label(r'log$_{10}$(M*/M$_{\odot}$)', size=12)
cbartree.set_label(r'Age (Gyr)', size=12)
#f1 = ax2.scatter(av_d4000_1, av_hdelta_1, marker = '*', color = 'r', s=300, zorder = 2, edgecolors = 'k',) #label = '10.5 - 10.75')
f2 = axtree2.scatter(av_d4000_2, av_hdelta_2, marker = '*', color = 'grey', s=250, zorder = 2, edgecolors = 'k') #label = '10.75 - 11.0')
f3 = axtree2.scatter(av_d4000_3, av_hdelta_3, marker = 'd', color = 'grey', s=250, zorder = 2,edgecolors = 'k' ) #label = '11.0 - 11.3')
f4 = axtree2.scatter(av_d4000_4, av_hdelta_4, marker = '^', color = 'grey', s=250, zorder = 2,edgecolors = 'k')
#s1 = axtree2.scatter(av_d4000_all_1, av_hdelta_all_1, marker = 'o', color = 'teal', s=300, zorder = 2, edgecolors = 'k', label = '10.4 - 10.5')
s2 = axtree2.scatter(av_d4000_2_bel, av_hdelta_2_bel, marker = '*', color = 'teal', s=250, zorder = 2, edgecolors = 'k', label = '10.5 - 10.75')
s3 = axtree2.scatter(av_d4000_3_bel, av_hdelta_3_bel, marker = 'd', color = 'teal', s=250, zorder = 2, edgecolors = 'k', label = '10.75 - 11.0')
s4 = axtree2.scatter(av_d4000_4_bel, av_hdelta_4_bel, marker = '^', color = 'teal', s=250, zorder = 2,edgecolors = 'k',label = '11.0 - 11.3')

axtree2.set_xlabel("D$_{n}$4000", size=13)
axtree2.set_ylabel("EW(H$\delta$)", size=13)
first_legend = plt.legend([s2, s3, s4], ['10.5 - 10.75', '10.75 - 11.0','11.0 - 11.3'], loc = 'lower right')
axtree2 = plt.gca().add_artist(first_legend)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='grey', label='Above A. van der Wel., 2014 \n (normalised) M* v R$_{e}$ relation')
green_patch = mpatches.Patch(color='teal', label='Below A. van der Wel., 2014 \n (normalised) M* v R$_{e}$ relation')
plt.legend(handles=[red_patch, green_patch], loc = 'upper right')
#plt.legend()
plt.savefig('SIZE_corrected_hdeltad4000cbar_medianmassbinned_agecbar78.pdf')
plt.close()

print("END")

input()

def plot_stackssingle(stack, name, color):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,1.75)
    plt.title('median stacks '+ str(name) +' line \n 9.7 < (log(M*) <= 11.4)', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_vdw_relation_ALL_TEST.pdf')
    plt.close()
    #plt.show()

#plot_stackssingle(stacking_all, 'all_newHdeltaEW', 'r')

#stacking_below = stacks(IDs_below)
from indices import C29_33, Dn4000, Mg_UV, H_delta
from stacking_code import stacks

def plot_stacks(stack1, stack2):
    plt.figure(figsize=(20,8))
    #plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below)})')
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.5, ls ='-', label = f'above relation (N = {len(IDs_above)})')
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below)})')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux ", size=17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title('Median stacks above and below normalised van der Wel ETG relation', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacks_abovebelow_vdw_relation_stacknormalisationremoved.pdf')
    plt.close()
    #plt.show()

#plot_stacks(stacking_above, stacking_below)


new_wavs = np.arange(2400, 4200, 1.25)
low_lim, upp_lim = 10.75, 11.0
low_lim2, upp_lim2 = 11.0, 11.3
low_lim3, upp_lim3 = 10.5, 10.75

cat1 = Table.read("Re_cat.fits").to_pandas()

def stack_lims(lower_lim, higher_lim):
    mass_mask = (masses>lower_lim) & (masses < higher_lim)
    df0 = pd.DataFrame(cat1[mass_mask])
    R_e = df0["Re_kpc"]
    index = np.log10(R_e).index.to_list()
    size  = np.array(np.log10(R_e))
    mass_ = df0['log10(M*/Msun)']
    x_array = np.linspace(lower_lim, higher_lim, len(size))
    vdw_norm_model = best_c_vdw + alpha*np.log10((10**x_array)/(5*10**10))
    mask = (size>vdw_norm_model)
    mask1 = (size<vdw_norm_model)
    index_masked_mass = np.log10(R_e)[mask].index.to_list()
    index_masked2_mass = np.log10(R_e)[mask1].index.to_list()

    IDs_above = IDs[index_masked_mass].str.decode("utf-8").str.rstrip().values
    IDs_below = IDs[index_masked2_mass].str.decode("utf-8").str.rstrip().values

    stacking_above = stacks(IDs_above)
    stacking_below = stacks(IDs_below)
    stacking_both = stacking_above, stacking_below

    return stacking_both


stacking_both_1075_11 = stack_lims(low_lim, upp_lim)
stacks_above_1075_11, stacks_below_1075_11 = stacking_both_1075_11[0], stacking_both_1075_11[1]
stacking_both_105_1075 = stack_lims(low_lim3, upp_lim3)
stacks_above_105_1075, stacks_below_105_1075 = stacking_both_105_1075[0], stacking_both_105_1075[1]
stacking_both_11_113= stack_lims(low_lim2, upp_lim2)
stacks_above_11_113, stacks_below_11_113 = stacking_both_11_113[0], stacking_both_11_113[1]

c_ind_below1075 = C29_33(stacks_below_1075_11)
c_ind_above1075 = C29_33(stacks_above_1075_11)
dn4000_above1075 = Dn4000(stacks_above_1075_11)
dn4000_below1075 = Dn4000(stacks_below_1075_11)
Mg_UV_index_below1075 = Mg_UV(stacks_below_1075_11)
Mg_UV_index_above1075 = Mg_UV(stacks_above_1075_11)
H_delta_EW_above1075 = H_delta(stacks_above_1075_11, new_wavs)
H_delta_EW_below1075 = H_delta(stacks_below_1075_11, new_wavs)

c_ind_below105 = C29_33(stacks_below_105_1075)
c_ind_above105 = C29_33(stacks_above_105_1075)
dn4000_above105 = Dn4000(stacks_above_105_1075)
dn4000_below105 = Dn4000(stacks_below_105_1075)
Mg_UV_index_below105 = Mg_UV(stacks_below_105_1075)
Mg_UV_index_above105 = Mg_UV(stacks_above_105_1075)
H_delta_EW_above105 = H_delta(stacks_above_105_1075, new_wavs)
H_delta_EW_below105 = H_delta(stacks_below_105_1075, new_wavs)

c_ind_below11 = C29_33(stacks_below_11_113)
c_ind_above11 = C29_33(stacks_above_11_113)
dn4000_above11 = Dn4000(stacks_above_11_113)
dn4000_below11 = Dn4000(stacks_below_11_113)
Mg_UV_index_below11 = Mg_UV(stacks_below_11_113)
Mg_UV_index_above11 = Mg_UV(stacks_above_11_113)
H_delta_EW_above11 = H_delta(stacks_above_11_113, new_wavs)
H_delta_EW_below11 = H_delta(stacks_below_11_113, new_wavs)
#print(Mg_UV_index_above11)

data = {'(log(M*/Msun))': ['10.5 - 10.75','  ','10.75 - 11.0','  ','11.0 - 11.3','  '],
        'above/below': ['above','below','above','below','above','below'],
        'EW(Hdelta)': [round(H_delta_EW_above105,5),round(H_delta_EW_below105,5),round(H_delta_EW_below1075,5),round(H_delta_EW_below1075,5),round(H_delta_EW_above11,5),round(H_delta_EW_below11,5)],
        'Mg_UV': [round(Mg_UV_index_above105,5),round(Mg_UV_index_below105,5),round(Mg_UV_index_above1075,5),round(Mg_UV_index_below1075, 5),round(Mg_UV_index_above11,5),round(Mg_UV_index_below11,5)],
        'C(29-33)': [round(c_ind_above105,5),round(c_ind_below105,5),round(c_ind_above1075,5),round(c_ind_below1075,5),round(c_ind_above11,5), round(c_ind_below11,5)],
        'Dn4000': [round(dn4000_above105,5),round(dn4000_below105,5),round(dn4000_above1075,5),round(dn4000_below1075,5),round(dn4000_above11,5),round(dn4000_below11,5)],
        }

df = pd.DataFrame(data, columns = ['(log(M*/Msun))','above/below', r'EW(Hdelta)', r'Mg_UV', 'C(29-33)', r'Dn4000'])
df['(log(M*/Msun))'] = df['(log(M*/Msun))'].astype(str)
df['Mg_UV'] = df['Mg_UV'].astype(float)
df['Dn4000'] = df['Dn4000'].astype(float)
#print (df)

table = Table.from_pandas(df)
from astropy.io import ascii
print(table)

ascii.write(table, 'indices_for_stacks.dat')

def plot_stackssingle(stack, name, color):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size=17) # $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,1.8)
    plt.title('median stacks '+ str(name) +f' line \n ({low_lim} < (log(M*) <= {upp_lim})', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_vdw_relation_'+str(low_lim)+'_'+str(upp_lim)+'_stacknormgone.pdf')
    plt.close()
    #plt.show()

def plot_stacks(stack1, stack2):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(figsize=(20,8))
    #plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below_1075_11)})')
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.5, ls ='-', label = f'above relation (N = {len(IDs_above_1075_11)})')
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below_1075_11)})')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size =17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title(f'Median stacks above and below normalised van der Wel ETG relation \n ({low_lim} < (log(M*) <= {upp_lim})', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacks_abovebelow_vdw_relation_'+str(low_lim)+'_'+str(upp_lim)+'_stacknormgone.pdf')
    plt.close()
    #plt.show()

#plot_stacks(stacking_above_1075_11, stacking_below_1075_11)

#plot_stackssingle(stacking_above_1075_11, 'above', 'r')
#plot_stackssingle(stacking_below_1075_11,'below', 'k')

input()
mask = (masses>10.4)
df2 = pd.DataFrame(cat1[mask])
R_e = df2['Re_kpc']
sigma_50 = (10**masses[mask])/(2*(R_e)**2)

redshifts = df2['redshifts']
mass_ = df2['log10(M*/Msun)']
#redshift = df0['redshifts']
#redshift = np.array(redshift[(mass_>10.4)])
print(len(redshifts))
#print(sigma_50)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig1, ax1 = plt.subplots(figsize=[12,8.5])
im1 = ax1.scatter(mass_, np.log10(sigma_50), s=130, c=redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
#cbar = fig1.colorbar(im1, ax=ax1)
y = 9.6*np.ones(len(mass_))
ax1.plot(mass_, y)
#cbar.set_label(r'Redshift (z)', size=12)
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax1.set_ylabel(r'$\mathrm{log_{10}{(\sigma_{50}/M_{\odot} kpc^{-2})}}$', size = 12)
plt.xticks(fontsize=10)
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(np.min(mass_), np.max(mass_))

plt.yticks(fontsize=10)
plt.title('Mass surface density versus log stellar mass', size = 13)
plt.savefig('sigma_50_v_mass_(masscomp).pdf')
plt.close()

sig_mask = (np.log10(sigma_50)>y)
print(sig_mask)
index_masked_sig= np.log10(R_e)[sig_mask].index.to_list()
IDs_top_half_sigma = IDs[index_masked_sig].str.decode("utf-8").str.rstrip().values
print(IDs_top_half_sigma)
print(len(IDs_top_half_sigma))

sig_mask_below = (np.log10(sigma_50)<=y)
index_masked_sig_below= np.log10(R_e)[sig_mask_below].index.to_list()
IDs_bottom_half_sigma = IDs[index_masked_sig_below].str.decode("utf-8").str.rstrip().values
print(len(IDs_bottom_half_sigma)) #32 obs below sig_50 half


stack_top_half_sig = stacks(IDs_top_half_sigma)

stack_bottom_half_sig = stacks(IDs_bottom_half_sigma)
def plot_stacks(stack1, stack2, name):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.figure(figsize=(20,8))
    #plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below_1075_11)})')
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.5, ls ='-', label = r' High $\mathrm{\sigma_{50}}$ (N = 46)')
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = r'Low $\mathrm{\sigma_{50}}$ (N = 32)')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size =17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title(r'Median stacks of low & high $\mathrm{\sigma_{50}}$', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig(str(name)+'.pdf')
    plt.close()

plot_stacks(stack_top_half_sig, stack_bottom_half_sig, 'stacks_bottom_and_top_half_sigma_50')

def plot_stacks_single(stack1, stack2, name1, name2):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.5, ls ='-', label = r' High $\mathrm{\sigma_{50}}$ (N = 46)')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size =17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    plt.legend(fontsize=14)
    plt.title(r'Median stack of high $\mathrm{\sigma_{50}}$', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig(str(name1)+'.pdf')
    plt.close()

    plt.figure(figsize=(20,8))
    #plt.plot(new_wavs, stack2*10**18, color="r", lw=1.5, ls ='-', label = r' High $\mathrm{\sigma_{50}}$ (N = 46)')
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = r'Low $\mathrm{\sigma_{50}}$ (N = 32)')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size =17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title(r'Median stack of low $\mathrm{\sigma_{50}}$', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig(str(name2)+'.pdf')
    plt.close()

plot_stacks_single(stack_top_half_sig, stack_bottom_half_sig, 'top_half_sigma_50', 'bottom_half_sigma_50')
### histogram for mass complete (ish) sample log10(M*) > 10.4
#print(np.max(masses))
fig, ax = plt.subplots(figsize = [12,8.5])
#bins = np.linspace(10.4, 11.3, 21)
#bins = [10.4, 10.45, 10.5, 10.55,10.6, 10.65, 10.7, 10.75, 10.8, 10.85, 10.9, 10.95, 11.0, 11.05, 11.1, 11.15, 11.2 , 11.25, 11.3] #,bins = bins,
ax.hist(np.log10(sigma_50),color='pink',ec="k")
ax.set_xlabel(r'$\mathrm{log_{10}{(\sigma_{50}/M_{\odot} kpc^{-2})}}$', size = 12)
ax.set_ylabel('N', size = 12)
ax.xaxis.grid(True, which='minor')
plt.savefig('histogram_sigma_50.pdf')
