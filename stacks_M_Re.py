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
col5 = fits.Column(name='R_c_size', format='E', array=Rc)
col6 = fits.Column(name='Rc_errs', format='E', array=Rc_errs)

hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
#file =  "TEST_CAT.fits"
#hdu.writeto(file)

print(f'Rc : {Rc}')

x = np.linspace(9.5, 11.5, len(Rc))
shen = 0.56*x + (-5.54)


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

#print(shen[20])

y_model = 0.56*x + (best_c)

#print(y_model[20])

#print(y_model - shen)
diff_shen_me = abs(y_model[0]) - abs(shen[0])
print(10**diff_shen_me)


plt.plot(x, shen, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
plt.plot(x, y_model, 'r', lw = 0.7, label=f'local ETG relation normalised by f = {round(10**diff_shen_me,2)}')
plt.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 8)
plt.ylabel(r'$\mathrm{log_{10}{(R_{c})}}$', size = 8)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('3D-HST log10(Rc) versus log10(M*/Msun)', size = 8)
plt.xlim(9.75, 11.5)
#plt.show()
plt.savefig('TEST_NEW_SHEN_ME.pdf')

cat = Table.read('TEST_CAT.fits').to_pandas()
df = pd.DataFrame(cat)
R_c = df["R_c_size"]
print(np.log10(R_c))
print(y_model)
index = np.log10(R_c).index.to_list()
print(index)
IDs = cat['IDs']
IDs_above = []
IDs_below = []
mass_ = df['log10(M*/Msun)']
size = np.array(np.log10(R_c))
y_model = np.array(y_model)

mask = (size>y_model)
mask1 = (size<y_model)
index_masked = np.log10(R_c)[mask].index.to_list()
index_masked2 = np.log10(R_c)[mask1].index.to_list()
#print(size.shape)
#print('size_mask1:',size[mask1])
#print(len(size[mask1]))
print('rcmask =',R_c[mask])
print(index_masked)
IDs_above = IDs[index_masked].str.decode("utf-8").str.rstrip().values
IDs_below= IDs[index_masked2].str.decode("utf-8").str.rstrip().values
#print((0.56*x + (best_c)).shape)
#print(mask)
#for i in size[mask]:

#    IDs_above.append(df.loc[(np.log10(R_c)[mask]==i),'IDs'].values.tolist())

print(len(IDs_above), len(IDs_below))


#stacking_above = stacks(IDs_above)
new_wavs = np.arange(2400, 4200, 1.25)

def plot_stackssingle(stack, name, color):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.1 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,2.0)
    plt.title('median stacks '+ str(name) +' line \n 9.7 < (log(M*) <= 11.4)', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_relation_ALL_TEST.pdf')
    plt.close()
    #plt.show()

#plot_stackssingle(stacking_above, 'above', 'r')


#stacking_below = stacks(IDs_below)

#plot_stackssingle(stacking_below,'below', 'k')


def plot_stacks(stack1, stack2):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.4, label = f'below relation (N = {len(IDs_below)})')
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.3, ls ='-', label = f'above relation (N = {len(IDs_above)})')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title('Median stacks above and below normalised Shen et al. 2003 ETG relation', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacks_abovebelow_shennorm_ALL_TEST.pdf')
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

mass_mask = (masses>10.5) & (masses <= 11.0)
cat1 = Table.read('TEST_CAT.fits').to_pandas()
df0 = pd.DataFrame(cat1[mass_mask])
R_c = df0["R_c_size"]
mass_ = df['log10(M*/Msun)']
size = np.array(np.log10(R_c))
y_model = np.array(y_model)
x_mask = (x>10.5)&(x<=11)
mask = (size>y_model[mass_mask])
mask1 = (size<y_model[mass_mask])
print(np.log10(R_c))
#print(IDs_above[mass_mask])


index_masked_mass= np.log10(R_c)[mask].index.to_list()
index_masked2_mass = np.log10(R_c)[mask1].index.to_list()
#print(size.shape)
#print('size_mask1:',size[mask1])
#print(len(size[mask1]))

IDs_above_105_75 = IDs[index_masked_mass].str.decode("utf-8").str.rstrip().values
IDs_below_105_75= IDs[index_masked2_mass].str.decode("utf-8").str.rstrip().values
#point = x0,y0

#if y0 > 0.56*x0 + best_c:
stacking_above_105_75 = stacks(IDs_above_105_75)
new_wavs = np.arange(2400, 4200, 1.25)

def plot_stackssingle(stack, name, color):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,2.0)
    plt.title('median stacks '+ str(name) +' line \n (10.5 < (log(M*) <= 10.75)', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_relation_105_75_TEST.pdf')
    plt.close()
    #plt.show()

plot_stackssingle(stacking_above_105_75, 'above', 'r')


stacking_below_105_75 = stacks(IDs_below_105_75)

plot_stackssingle(stacking_below_105_75,'below', 'k')


def plot_stacks(stack1, stack2):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below_105_75)})')
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.5, ls ='-', label = f'above relation (N = {len(IDs_above_105_75)})')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title('Median stacks above and below normalised Shen et al. 2003 ETG relation \n (10.5 < (log(M*) <= 10.75)', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacks_abovebelow_shennorm_105_75_TEST.pdf')
    plt.close()
    #plt.show()

plot_stacks(stacking_above_105_75, stacking_below_105_75)



#𝑦0>𝑎𝑥0+𝑏 then the point is above the line, etc.
