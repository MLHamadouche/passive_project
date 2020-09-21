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

#concat_3dhst = Table.read('FirstProjectCatalogs/concat_candels_passive_match.fits').to_pandas()
concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df = pd.DataFrame(concat_3dhst)
ID_list1 = np.array(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
ID_ = df.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
agn = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))

vandels_cat_pipes = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()
df1 = pd.DataFrame(vandels_cat_pipes)

ID_pipes = vandels_cat_pipes['#ID'].values
IDs = df1.set_index(s.decode('utf-8') for s in vandels_cat_pipes['#ID'])

#print(IDs)

UV = vandels_cat_pipes['UV_colour_50']
VJ = vandels_cat_pipes['VJ_colour_50']
age1 = vandels_cat_pipes['mass_weighted_age_50']
ssfr = vandels_cat_pipes['ssfr_50']
re = concat_3dhst['re']
redshifts = concat_3dhst['zspec']

mass = passive_cut["log10(M*)"]

#fig, ax1 = plt.subplots(figsize=[12,7])
list_IDs = []
all_ages = []
all_masses = []
all_sizes = []
list_IDs = []
size_errs =[]

cosmo  = FlatLambdaCDM(H0=70, Om0=0.3)
Mpc_to_kpc = 1000
for i in ID_list1:
    for j in ID_pipes:
        j = j.decode("utf-8")
        if i == j:
            list_IDs.append(i)
new_IDs = []
new_redshifts = []
size_kpc = []
size_kpc = []
dA=[]
q_ratio = []
for ID in list_IDs:
    if ID not in agn:
        new_redshifts.append(ID_.loc[ID, "zspec"])

        all_ages.append(IDs.loc[ID, "mass_weighted_age_50"])
        #print(all_ages)
        all_masses.append(ID_list.loc[ID, "log10(M*)"])
        #print(all_masses)
        all_sizes.append(ID_.loc[ID, 're'])
        size_errs.append(ID_.loc[ID, 'dre'])
        q_ratio.append(ID_.loc[ID, 'q'])
        new_IDs.append(ID)
#print(all_sizes)


arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(np.array(new_redshifts))

print(np.array(new_redshifts))


print(arcsec_per_kpc)


"""
for z in new_redshifts:
    dA = cosmo.angular_diameter_distance(z)
    all_sizes = np.array(all_sizes)*u.arcsec
    size_kpc=((all_sizes*dA)).to(u.kpc, u.dimensionless_angles())
    size_kpc_errs=((np.array(size_errs)*u.arcsec)*dA).to(u.kpc,u.dimensionless_angles())
#print(len(all_ages)) #95 objects out of 228 with sizes from the 3D-HST catalog crossmatch
"""

size_kpc = (np.array(all_sizes)*u.arcsec)/arcsec_per_kpc
size_kpc_errs = (np.array(size_errs)*u.arcsec)/arcsec_per_kpc


#print(all_sizes, '\n', size_kpc, '\n', size_kpc_errs)
R_c = (np.sqrt(np.array(q_ratio))*size_kpc) /u.kpc
Rc_errs = (np.sqrt(np.array(q_ratio))*size_kpc_errs) /u.kpc

#print(R_c)

#from Shen et al 2009
#print(x)

#print(y_errs)
#actually in arcseconds

#one_sig_scatter = np.std(shen)
plt.scatter(np.array(all_masses), np.log10(10**np.array(all_masses)/(R_c)**2), marker='o', s=20, c='r', edgecolors='k')
#plt.plot(x_values, shen-1, linewidth=1.2, color='r', label=r'Shen et al. 2003 ETG relation with 1- $\mathrm{\sigma}$ scatter')
#plt.plot(x_values, (shen+one_sig_scatter)-1, linewidth=0.5, color='r', linestyle='--')
#plt.plot(x_values, (shen-one_sig_scatter)-1, linewidth=0.5, color='r', linestyle='--')
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 8)
plt.ylabel(r'$\mathrm{log_{10}{(M/R_{c}^{2}})}$', size = 8)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('95 3DHST passive galaxies mass-density ', size = 8)
plt.xlim(9.75, 11.5)
#plt.show()

#plt.savefig('Rc_mass_relation_3DHST_density.pdf')
#print(10**all_sizes)
plt.close()

def model(m,c):
    return m * x_values + c

x_values = np.linspace(9.7, 11.4, 95)

shen = model(0.56, -5.54) #-1.2
one_sig_scatter = np.std(shen)
fig, ax = plt.subplots(figsize=[10,7.5])

im = ax.scatter(np.array(all_masses), np.log10(R_c), s=50, c=new_redshifts, cmap=plt.cm.inferno, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig.colorbar(im, ax=ax)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Redshift, z', size=12)
#plt.scatter(np.array(all_masses), np.log10(R_c), marker='o', s=20, c='r', edgecolors='k')
plt.plot(x_values, shen-np.log10(2.43), linewidth=1.2, color='r', label=r'Shen et al. 2003 ETG relation with 1- $\mathrm{\sigma}$ scatter')
plt.plot(x_values, (shen+one_sig_scatter)-np.log10(2.43), linewidth=0.5, color='r', linestyle='--')
plt.plot(x_values, (shen-one_sig_scatter)-np.log10(2.43), linewidth=0.5, color='r', linestyle='--')
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 13)
plt.ylabel(r'$\mathrm{log_{10}{(R_{c}/kpc})}$', size = 13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 11})
plt.title('95 3DHST passive galaxies overplot with Shen et al. 2003 relation ', size = 13)
plt.xlim(9.75, 11.4)
#plt.show()
#plt.savefig('redshifts_R_c_log10M_kpcproper_relation_3DHST.pdf')
#print(10**all_sizes)
plt.close()

#stack galaxies above and below the lines


#rewrite this code and make sure you get the IDs and all in a table - eg. make a new table for the objects
#and the ones with sizes to make shit easier so youre not always screwing about with the crap on lines 76 - 87

#cross product for the above and below the lines


def model(m,c):
    return m * x_values + c #from Shen et al 2009


def model(c):

    mod_vals = (0.56 * x_values) + c

    return  mod_vals

#m_values = np.arange(-5.0, 5.0 ,0.01)

c_values = np.arange(-10., -1.1, 0.1)

best_chi  = np.inf
#best_m = a
#best_c = b
#print(0.434*(Rc_errs/R_c))
#print(np.log10(R_c))
log_Rc = np.log10(R_c)

errors = (np.log10(R_c + 10*Rc_errs) - np.log10(R_c - 10*Rc_errs))/2
#print(errors)

SNR = R_c/Rc_errs

#print('SNR = ', R_c/Rc_errs)


errs = 0.434*(10*Rc_errs/R_c)
#print(log_Rc)
#print(errs)

#for mvals in range(len(m_values)):
#    m_vals = m_values[mvals]
for cvals in range(len(c_values)):
    c_vals = c_values[cvals]

    y_model = model(c_vals)

    diffs = (np.log10(R_c) - y_model)

    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq = np.sum((diffs**2)/((errs)**2))

    if chisq < best_chi:
        #best_m = m_vals
        best_c = c_vals
        best_chi = chisq

print(f'best_c {best_c} \n best_chi {best_chi}')
#print('mean SNR = ', R_c/Rc_errs)

y_model2= model(best_c)



def nmodel(m,c):
    return m * x_values + c #from Shen et al 2009

topcat = nmodel(0.5192562, -6.1182485)
shen = nmodel(0.56, -5.54)#- np.log10(1.5)

#print((shen - y_model2))
#print(np.log10(2.43))
print(np.std(y_model2))
print(np.std(shen))

one_sig_scatter = np.std(shen)
fig1, ax1 = plt.subplots(figsize=[12,8.5])

im1 = ax1.scatter(np.array(all_masses), np.log10(R_c), s=50, c=new_redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Redshift, z', size=12)
#plt.scatter(np.array(all_masses), np.log10(R_c), marker='o', s=20, c='r', edgecolors='k')
#ax1.plot(x_values, y_model2, linewidth=2, color='r', label='Best fit Shen et al. normalised by $\mathrm{f_{g} = 1.82}$')
#ax1.plot(x_values, y_model2+np.std(y_model2), linewidth=1., color='r', ls = ':') #label='Best fit Shen et al. normalised by $\mathrm{f_{g} = 1.82}$')
#ax1.plot(x_values, y_model2-np.std(y_model2), linewidth=1., color='r', ls = ':')#label='Best fit Shen et al. normalised by $\mathrm{f_{g} = 1.82}$')
ax1.plot(x_values, shen, linewidth=2, color='r', label='Shen et al. 2003 ETG relation with 1-$\mathrm{\sigma}$ scatter' )#\n & McLure et al. 2013 ($\mathrm{f_{g} = 2.43}$)')
ax1.plot(x_values, (shen+one_sig_scatter), linewidth=1., color='r', linestyle=':')
ax1.plot(x_values, (shen-one_sig_scatter), linewidth=1., color='r', linestyle=':')
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 13) #-np.log10(2.43)
ax1.set_ylabel(r'$\mathrm{log_{10}{(R_{c}/kpc})}$', size = 13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 11}) #$\mathrm{R_{c}}$
plt.title('95 passive VANDELS galaxies with 3D-HST half-light radii', size = 13)
ax1.set_xlim(9.75, 11.4)
#plt.show()
#plt.savefig('redshifts_R_c_log10M_onlyshen_no_norm_3DHST.pdf')
plt.close()

stack_obs_above = []
stack_obs_below=[]
objID = vandels_cat_pipes["#ID"]

#print(np.log10(R_c).value)#, '\n', np.transpose(np.log10(R_c)))

size_mass = np.transpose(np.array([all_masses, np.log10(R_c).value]))
#shen_x = np.transpose(np.array([x_values, shen]))
y_model2_x = np.transpose(np.array([x_values, y_model2]))

#print(shen_x)
cross = np.cross(size_mass, y_model2_x)

#print(cross)
#print(new_IDs)

#print(len(np.log10(R_c).value[mask]))
#ind = all_masses.value[mask]

#print(ind)

#print(new_IDs.index(np.log10(R_c).value[mask]))

#print(np.log10(R_c).value[mask].index(all_sizes))
#for ob in np.log10(R_c):
#    ind.append(list(np.log10(R_c).value[mask]).index(ob.value))
#print(ind)
    #stack_obs_above.append(new_IDs[ind])



col1 = fits.Column(name='ID', format='30A', array=new_IDs)
col2 = fits.Column(name='redshifts', format='E', array=new_redshifts)
col3 = fits.Column(name='age', format='E',  array=all_ages)
col4 = fits.Column(name='log10(M*)', format='E', array=all_masses)
col5 = fits.Column(name='R_c', format='E', array=R_c)
col6 = fits.Column(name='Rc_errs', format='E', array=Rc_errs)


hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
#file =  "R_c_v_mass_passive_cat.fits"
#hdu.writeto(file)
#print(stack_obs_above)

cat = Table.read("R_c_v_mass_passive_cat.fits").to_pandas()
masses = cat['log10(M*)']

mass_mask = (masses > 10.5) & (masses < 11)
dataframe = pd.DataFrame(cat)
mask = cross[mass_mask] > 0
cat_mask = pd.DataFrame(cat[mass_mask])
circ_radii = cat_mask['R_c']
print(len(circ_radii))
ind = cat[mass_mask].set_index(cat_mask['ID'].str.decode("utf-8").str.rstrip())

index = np.log10(circ_radii)[mask].index.to_list()

mask2 = cross[mass_mask] < 0
index2 = np.log10(circ_radii)[mask2].index.to_list()


#print(len(cat_mask))

for i in index:
    stack_obs_above.append(cat_mask['ID'].str.decode("utf-8").str.rstrip()[i])

for id in index2:
    stack_obs_below.append(cat_mask['ID'].str.decode("utf-8").str.rstrip()[id])


print(len(stack_obs_below))
print(len(stack_obs_above))
#for ID in stack_obs_above:
#    z = ID_.loc[ID, 'zspec']

#below = 32 objects being stacked
#above = 28 objects being stacked

stacking_above = stacks(stack_obs_above)
new_wavs = np.arange(2400, 4200, 1.25)

def plot_stackssingle(stack, name):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, stack*10**18, color="black", lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    #plt.xlim(2300, 4250)
    plt.ylim(0 ,2.0)
    plt.title('median stacks '+ str(name) +' line \n (10.5 < log(M*) < 11.0)')# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_relation.pdf')
    plt.close()

plot_stackssingle(stacking_above, 'above')


stacking_below = stacks(stack_obs_below)

plot_stackssingle(stacking_below,'below')

#for ID in list_IDs:
    #new_re.append(ID_.loc[ID, 're'])
    #new_UV.append(IDs.loc[ID, 'UV_colour_50'])
    #new_VJ.append(IDs.loc[ID, 'VJ_colour_50'])
