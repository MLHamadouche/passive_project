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

both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()

df1 = pd.DataFrame(both_xray)#, index = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8")))

passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()

df2 = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )

new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))

#candels  = Table.read('FirstProjectCatalogs/concat_candels_passive_match.fits').to_pandas()
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
for z in new_redshifts:
    dA = cosmo.angular_diameter_distance(z)
    all_sizes = np.array(all_sizes)*u.arcsec
    size_kpc=((all_sizes*dA)).to(u.kpc, u.dimensionless_angles())
    size_kpc_errs=((np.array(size_errs)*u.arcsec)*dA).to(u.kpc,u.dimensionless_angles())
#print(len(all_ages)) #95 objects out of 228 with sizes from the 3D-HST catalog crossmatch
#print(all_sizes, '\n', size_kpc, '\n', size_kpc_errs)
R_c = (np.sqrt(np.array(q_ratio))*size_kpc) /u.kpc
Rc_errs = (np.sqrt(np.array(q_ratio))*size_kpc_errs) /u.kpc

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
plt.ylabel(r'$\mathrm{log_{10}{(M/R_{e}^{2}})}$', size = 8)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('95 3D-HST passive galaxies mass-density ', size = 8)
plt.xlim(9.75, 11.5)
plt.savefig('Rc_mass_relation_3DHST_density.pdf')
#print(10**all_sizes)
plt.close()
"""
shen = model(0.56, -5.54) #-1.2
one_sig_scatter = np.std(shen)
plt.scatter(np.array(all_masses), np.log10(all_sizes, where=0<all_sizes, out=np.nan*all_sizes), marker='o', s=20, c='r', edgecolors='k')
plt.plot(x_values, shen-1, linewidth=1.2, color='r', label=r'Shen et al. 2003 ETG relation with 1- $\mathrm{\sigma}$ scatter')
plt.plot(x_values, (shen+one_sig_scatter)-1, linewidth=0.5, color='r', linestyle='--')
plt.plot(x_values, (shen-one_sig_scatter)-1, linewidth=0.5, color='r', linestyle='--')
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 8)
plt.ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc})}$', size = 8)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('90 CANDELS passive galaxies overplot with Shen et al. 2003 relation ', size = 8)
plt.xlim(9.75, 11.5)
plt.savefig('Size_mass_relation_CANDELS_offset-1.pdf')
#print(10**all_sizes)
def model(m,c):
    return m * x_values + c #from Shen et al 2009

"""

x_values = np.linspace(9.75, 11.3, 95)

def model(m,c):

    mod_vals = (m * x_values) + c

    return  mod_vals

m_values = np.arange(0.1, 1.0 ,0.01)

c_values = np.arange(-0.85, 0.85, 0.001)

best_chi  = np.inf
#best_m = a
#best_c = b
#print(0.434*(Rc_errs/R_c))
#print(np.log10(R_c))
log_Rc = np.log10(R_c)
errs = 0.434*(Rc_errs/R_c)
print(log_Rc)
#print(errs)

for mvals in range(len(m_values)):
    m_vals = m_values[mvals]
    for cvals in range(len(c_values)):
        c_vals = c_values[cvals]

        y_model = model(m_vals, c_vals)

        diffs = (np.log10(R_c) - y_model)

        #print(diffs)
        #print(y_model, '\n', np.log10(R_c))
        chisq = np.sum(diffs**2/((0.434*(Rc_errs/R_c))**2))

        if chisq < best_chi:
            best_m = m_vals
            best_c = c_vals
            best_chi = chisq

print(f'best_m {best_m} \n best_c {best_c} \n best_chi {best_chi}')

#def model(m,c):
#    return m * np.linspace(9.8, 11.5, 95) + c #from Shen et al 2009
#ymodel = log_model(best_m, best_c)
ymodel=model(best_m, best_c)
topcat = model(0.5192562, -6.1182485)
shen = (model(0.56, -5.54)
one_sig_scatter = np.std(shen)

plt.scatter(all_masses, np.log10(R_c), marker='o', s=20, c='r', edgecolors='k')
plt.plot(x_values, ymodel, linewidth=1.2, color='k', label='Best fit line', )
plt.plot(x_values, shen- np.log10(1.5)), linewidth=1.2, color='r', label=r'Shen et al. 2003 ETG relation with 1- $\mathrm{\sigma}$ scatter')
plt.plot(x_values, (shen+one_sig_scatter)- np.log10(1.5)), linewidth=0.5, color='r', linestyle='--')
plt.plot(x_values, (shen-one_sig_scatter)- np.log10(1.5)), linewidth=0.5, color='r', linestyle='--')
plt.errorbar(np.array(all_masses), np.log10(R_c), yerr = 0.434*(Rc_errs/R_c), linestyle=' ')
#plt.plot(x_values, topcat, linewidth=1.0, color='g'  (1/2.303)*(Rc_errs/R_c)
plt.xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 8)
plt.ylabel(r'$\mathrm{log_{10}{(R_{c}/kpc})}$', size = 8)

plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(prop={'size': 7})
plt.title('95 passive galaxies overplot with Shen et al. 2003 relation and best fit line', size = 8)
plt.xlim(9.8, 11.3)
plt.ylim(-0.75, 1.5)
plt.savefig('Rc_mass_relation_3dhst_shen+bestfit.pdf')
plt.close()
stack_obs_above = []
stack_obs_below=[]
objID = vandels_cat_pipes["#ID"]







"""
cross = np.cross(all_sizes, shen)
for s in cross:
    if s > 0:
        stack_obs_above.append(new_IDs[s])


print(stack_obs_above)
#stacks1 = all_sizes[all_sizes<(np.array(shen)-1)]

#print(len(stacks1))


#print(shen-1)
for vals in np.array(shen):
    for obj in all_sizes:
        if obj > vals:
            stack_obs_above.append(new_IDs.index(objs))#.decode("utf-8")
        else:
            stack_obs_below.append(new_IDs.index(objs))#.decode("utf-8"))

print(stack_obs_above, '\n', stack_obs_below)

"""


#for ID in list_IDs:
    #new_re.append(ID_.loc[ID, 're'])
    #new_UV.append(IDs.loc[ID, 'UV_colour_50'])
    #new_VJ.append(IDs.loc[ID, 'VJ_colour_50'])
