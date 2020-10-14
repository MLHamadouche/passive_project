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

UV = vandels_cat_pipes['UV_colour_50']
VJ = vandels_cat_pipes['VJ_colour_50']
age1 = vandels_cat_pipes['mass_weighted_age_50']
ssfr = vandels_cat_pipes['ssfr_50']
re = concat_3dhst['re']
redshifts = concat_3dhst['zspec']

mass = passive_cut["log10(M*)"]

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
for ID in list_IDs:
    if ID not in agn:
        new_redshifts.append(ID_.loc[ID, "zspec"])

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

hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
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
x_model = np.linspace(10.4, 11.0, len(Rc))

c_model = np.arange(-8.0, -3.0, 0.01)

best_chi  = np.inf

errs = 0.434*(Rc_errs/Rc)

for cvals in range(len(c_model)):
    c_vals = c_model[cvals]

    y_model = 0.56*x_model + c_vals
    diffs = y_model - np.log10(Rc)
    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq = np.sum((diffs**2)/((10*errs)**2))

    if chisq < best_chi:
        best_c = c_vals
        best_chi = chisq

print(f'best_c {best_c} \n best_chi {best_chi}')
alpha = 0.76
log_A = 0.22
log_Reff = log_A + alpha*np.log10(10**x/(5*10**10))
#print(shen[20])
print(f'std mass = {np.std(masses)}')
print(f'std Re = {np.std(Re_kpc_errs/u.kpc)}')
log_A_model = np.arange(0.0, 0.3, 0.01)
x_model2 = np.linspace(10.3, 11.0, len(Re_kpc/u.kpc))
best_chi2 = np.inf

for cvals2 in range(len(log_A_model)):
    c_vals2 =log_A_model[cvals2]

    vdw_model = alpha*np.log10(10**x_model2/(5*10**10)) + c_vals2
    diffs2 = vdw_model - np.log10(Re_kpc/u.kpc)
    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq2 = np.sum((diffs2**2)/((4.34*(Re_kpc_errs/Re_kpc))))

    if chisq2 < best_chi2:
        best_c_vdw = c_vals2
        best_chi_vdw = chisq2

print(f'best_c_vdw: {best_c_vdw} \n best_chi_vdw: {best_chi_vdw}')
print(vdw_model)

vdw_norm_model= best_c_vdw + alpha*np.log10(10**x/(5*10**10))

y_model = 0.56*x + (best_c)

print(best_c_vdw - 0.22)
#print(y_model[20])

#print(y_model - shen)
diff_shen_me = abs(y_model[0]) - abs(shen[0])
print(10**diff_shen_me)


fig1, ax1 = plt.subplots(figsize=[12,8.5])
im1 = ax1.scatter(masses, np.log10(Re_kpc/u.kpc), s=130, c=new_redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Redshift (z)', size=12)
ax1.plot(x, y_model, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
ax1.plot(x, vdw_norm_model, 'k', lw = 1.2, label=f'local ETG relation normalised by f = {round(10**diff_shen_me,2)}')
ax1.plot(x, log_Reff, 'r' ,lw = 1.2, label= 'A. van der Wel z = 1.25 ETG relation')
#ax1.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax1.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(prop={'size': 10})
plt.title('3D-HST log10(Re) versus log10(M*/Msun)', size = 13)
plt.xlim(9.75, 11.4)
#plt.show()
plt.savefig('Re_v_M*_test_with_Arjens_relationandbestfit.pdf')
plt.close()


#cat = Table.read('TEST_CAT.fits').to_pandas()
cat = Table.read("Re_cat.fits").to_pandas()
df = pd.DataFrame(cat)
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

print(all_IDs)

from indices import C29_33, Dn4000, Mg_UV # H_delta

colour_index = []
D4000_index =[]
Mg_UV_index = []
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


#print((0.56*x + (best_c)).shape)
#print(mask)
#for i in size[mask]:

#    IDs_above.append(df.loc[(np.log10(R_c)[mask]==i),'IDs'].values.tolist())

#print(len(IDs_above), len(IDs_below))
new_wavs = np.arange(2400, 4200, 1.25)
stacking_all = stacks(all_IDs)
#stacking_above = stacks(IDs_above)


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
input()
#plot_stackssingle(stacking_above, 'above', 'r')

#stacking_below = stacks(IDs_below)
from indices import C29_33, Dn4000, Mg_UV, H_delta
"""
c_ind = C29_33(stacking_below)
print(f'C(29-33) index below: {c_ind}')
c_ind = C29_33(stacking_above)
print(f'C(29-33) index above : {c_ind}')
dn4000 = Dn4000(stacking_above)
print(f'Dn4000 index above: {dn4000}')
dn4000 = Dn4000(stacking_below)
print(f'Dn4000 index below: {dn4000}')

Mg_UV_index = Mg_UV(stacking_below)
print(f'Mg_UV index below: {Mg_UV_index}')
#plot_stackssingle(stacking_below,'below', 'k')
Mg_UV_index = Mg_UV(stacking_above)
print(f'Mg_UV index above: {Mg_UV_index}')

#H_delta_EW = H_delta(stacking_above)
#print(f'H_delta_EW above: {H_delta_EW}')

#H_delta_EW = H_delta(stacking_below)
#print(f'H_delta_EW above: {H_delta_EW}')
"""


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

"""

if np.greater(size,np.array(0.56*x + np.ones(95)*(best_c)))==True:
    #it's on the left, above the line
    IDs_above.append(df.loc[np.log10(df['R_c_size']) == round(s,8),'IDs'].values.tolist())

if np.less(size,np.array(0.56*x + np.ones(95)*(best_c)))==True:
    #it's on the right, below the line
    IDs_below.append(df.loc[np.log10(df['R_c_size']) == round(s,8),'IDs'].values.tolist()) #.str.decode("utf-8").str.rstrip()
"""
#print(f'IDs_above {IDs_above}')
#print(f'IDs_below {IDs_below}')
low_lim = 10.4
upp_lim = 11.3
mass_mask = (masses>low_lim) & (masses <= upp_lim)
cat1 = Table.read("Re_cat.fits").to_pandas()

df0 = pd.DataFrame(cat1[mass_mask])
R_e = df0["Re_kpc"]
index = np.log10(R_e).index.to_list()
size  = np.array(np.log10(R_e))
#y_model = np.array(y_model)
#vdw_norm_model = np.array(vdw_norm_model)
#cat1 = Table.read('TEST_CAT.fits').to_pandas()

#R_c = df0["R_c_size"]
mass_ = df0['log10(M*/Msun)']
#size = np.log10(R_e)
#y_model = np.array(y_model)
#x_mask = (x>10.5)&(x<=10.75)
x_array= np.linspace(low_lim, upp_lim, len(size))
#y_model  = 0.56*x_array+(best_c)
vdw_norm_model = best_c_vdw + alpha*np.log10((10**x_array)/(5*10**10))

print(len(vdw_norm_model))
print(len(R_e))
#mask = (size>y_model)
#mask1 = (size<y_model)

mask = (size>vdw_norm_model)
mask1 = (size<vdw_norm_model)

#print(IDs_above[mass_mask])



index_masked_mass= np.log10(R_e)[mask].index.to_list()
index_masked2_mass = np.log10(R_e)[mask1].index.to_list()
#print(size.shape)
#print('size_mask1:',size[mask1])
#print(len(size[mask1]))

IDs_above_1075_11 = IDs[index_masked_mass].str.decode("utf-8").str.rstrip().values
IDs_below_1075_11= IDs[index_masked2_mass].str.decode("utf-8").str.rstrip().values
#point = x0,y0
#print(IDs_above_1075_11)
#print(IDs_below_1075_11)

#if y0 > 0.56*x0 + best_c:
stacking_above_1075_11 = stacks(IDs_above_1075_11)
new_wavs = np.arange(2400, 4200, 1.25)

def plot_stackssingle(stack, name, color):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size=17) # $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,1.8)
    plt.title('median stacks '+ str(name) +f' line \n ({low_lim} < (log(M*) <= {upp_lim})', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_vdw_relation_'+str(low_lim)+'_'+str(upp_lim)+'_stacknormalisationremoved.pdf')
    plt.close()
    #plt.show()

plot_stackssingle(stacking_above_1075_11, 'above', 'r')



stacking_below_1075_11 = stacks(IDs_below_1075_11)

plot_stackssingle(stacking_below_1075_11,'below', 'k')


def plot_stacks(stack1, stack2):
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
    plt.savefig('stacks_abovebelow_vdw_relation_'+str(low_lim)+'_'+str(upp_lim)+'_stacknormalisationremoved.pdf')
    plt.close()
    #plt.show()

plot_stacks(stacking_above_1075_11, stacking_below_1075_11)

"""
c_ind = C29_33(stacking_below_1075_11)
print(f'C(29-33) index below: {c_ind}')
c_ind = C29_33(stacking_above_1075_11)
print(f'C(29-33) index above : {c_ind}')
dn4000 = Dn4000(stacking_above_1075_11)
print(f'Dn4000 index above: {dn4000}')
dn4000 = Dn4000(stacking_below_1075_11)
print(f'Dn4000 index below: {dn4000}')


Mg_UV_index = Mg_UV(stacking_below)
print(f'Mg_UV index below: {Mg_UV_index}')
#plot_stackssingle(stacking_below,'below', 'k')
Mg_UV_index = Mg_UV(stacking_above)
print(f'Mg_UV index above: {Mg_UV_index}')
#𝑦0>𝑎𝑥0+𝑏 then the point is above the line, etc.
H_delta_EW = H_delta(stacking_above)
print(f'H_delta_EW above: {H_delta_EW}')

H_delta_EW = H_delta(stacking_below)
print(f'H_delta_EW above: {H_delta_EW}')
"""
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

fig1, ax1 = plt.subplots(figsize=[12,8.5])
im1 = ax1.scatter(mass_, np.log10(sigma_50), s=130, c=redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
#cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'Redshift (z)', size=12)
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax1.set_ylabel(r'$\mathrm{log_{10}{(\sigma_{50}/M_{\odot} kpc^{-2})}}$', size = 12)
plt.xticks(fontsize=10)
#plt.yscale('log')
#plt.xscale('log')
plt.yticks(fontsize=10)
plt.title('Mass surface density versus log stellar mass', size = 13)
plt.savefig('sigma_50_v_mass_(masscomp).pdf')
plt.close()


### histogram for mass complete (ish) sample log10(M*) > 10.4
print(np.max(masses))
fig, ax = plt.subplots(figsize = [12,8.5])
bins = np.linspace(10.4, 11.3, 21)
#bins = [10.4, 10.45, 10.5, 10.55,10.6, 10.65, 10.7, 10.75, 10.8, 10.85, 10.9, 10.95, 11.0, 11.05, 11.1, 11.15, 11.2 , 11.25, 11.3] #,bins = bins,
ax.hist(np.log10(sigma_50), bins = bins, color='pink',ec="k")
plt.savefig('histogram_sigma_50.pdf')
