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
#new IDs has passive galaxies excluding possible AGN found in the 3DHST catalog/CANDELS catalog

arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(np.array(new_redshifts))
print('max=', max(all_sizes))
#print(np.array(new_redshifts))


#print(arcsec_per_kpc)


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

#x_values = np.linspace(10.3, 10.8, 95)
x_values0 = np.linspace(10., 11., 95)

def model(m,c):
    return m * x_values0 + c


shen = model(0.56, -5.54) #-1.2
one_sig_scatter = np.std(shen)
fig, ax = plt.subplots(figsize=[10,7.5])

im = ax.scatter(np.array(all_masses), np.log10(R_c), s=50, c=new_redshifts, cmap=plt.cm.inferno, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig.colorbar(im, ax=ax)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Redshift, z', size=12)
#plt.scatter(np.array(all_masses), np.log10(R_c), marker='o', s=20, c='r', edgecolors='k')
plt.plot(x_values0, shen-np.log10(2.43), linewidth=1.2, color='r', label=r'Shen et al. 2003 ETG relation with 1- $\mathrm{\sigma}$ scatter')
plt.plot(x_values0, (shen+one_sig_scatter)-np.log10(2.43), linewidth=0.5, color='r', linestyle='--')
plt.plot(x_values0, (shen-one_sig_scatter)-np.log10(2.43), linewidth=0.5, color='r', linestyle='--')
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



def model(c):

    mod_vals = (0.56 * x_values1) + c

    return  mod_vals

m_values = np.arange(0.5, 0.6 ,0.01)

c_values = np.arange(-6., -5., 0.01)

best_chi  = np.inf
#best_m = a
#best_c = b
#print(0.434*(Rc_errs/R_c))
#print(np.log10(R_c))
log_Rc = np.log10(R_c)

errors = (np.log10(R_c + 10*Rc_errs) - np.log10(R_c - 10*Rc_errs))/2
#print(errors)

SNR = R_c/Rc_errs

print('SNR = ', R_c/Rc_errs)


errs = 0.434*(10*Rc_errs/R_c)
#print(log_Rc)
#print(errs)
#errs = 0.33*R_c
#for mvals in range(len(m_values)):
#    m_vals = m_values[mvals]
for cvals in range(len(c_values)):
    c_vals = c_values[cvals]

    y_model = 0.56*x_values0 + c_vals

    diffs = (np.log10(R_c) - y_model)

    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq = np.sum((diffs**2)/((errs)**2))

    if chisq < best_chi:
        #best_m = m_vals
        best_c = c_vals
        best_chi = chisq

print(f'best_c {best_c} \n best_chi {best_chi}')# \n best_m {best_m}')
#print('mean SNR = ', R_c/Rc_errs)

#y_model2 = model(best_c)

#new_x = np.linspace(9.8, 11.4, 95)
x_values1 = np.linspace(9.7, 11.5, 95)

def nmodel(m,c):
    return m * x_values1 + c #from Shen et al 2009

#y_model2 = 0.56*x_values1 + best_c
#y_model2 = best_m*x_values1 + best_c
y_model2 = 0.56*x_values1 + best_c
topcat = nmodel(0.5192562, -6.1182485)
shen = nmodel(0.56, -5.54)#- np.log10(1.5)

print((shen - y_model2))
print(10**(shen - y_model2))
#print(np.log10(2.43))
#print('YMODEL2:', y_model2)
#print(np.log10(R_c).value)

#print(np.std(y_model2))
#print(np.std(shen))
#x_values1 = np.linspace(9.7, 11.5, 95)
from matplotlib.offsetbox import AnchoredText
one_sig_scatter = np.std(shen)
fig1, ax1 = plt.subplots(figsize=[12,8.5])

im1 = ax1.scatter(np.array(all_masses), np.log10(R_c), s=50, c=new_redshifts, cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig1.colorbar(im1, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Redshift, z', size=12)
#plt.scatter(np.array(all_masses), np.log10(R_c), marker='o', s=20, c='r', edgecolors='k')
ax1.plot(x_values1, y_model2, linewidth=2, color='r', label=f'Best fit Shen et al. normalised by f = {round(np.mean(10**(shen - y_model2)),2)}')
ax1.plot(x_values1, y_model2+np.std(y_model2), linewidth=1., color='r', ls = ':')#label='Best fit Shen et al. normalised by $\mathrm{f_{g} = 1.82}$')
ax1.plot(x_values1, y_model2-np.std(y_model2), linewidth=1., color='r', ls = ':')#label='Best fit Shen et al. normalised by $\mathrm{f_{g} = 1.82}$')
ax1.plot(x_values1, shen, linewidth=2, color='k', label='Shen et al. 2003 local ETG relation with 1-$\mathrm{\sigma}$ scatter' )#\n & McLure et al. 2013 ($\mathrm{f_{g} = 2.43}$)')
ax1.plot(x_values1, (shen+one_sig_scatter), linewidth=1., color='k', linestyle=':')
ax1.plot(x_values1, (shen-one_sig_scatter), linewidth=1., color='k', linestyle=':')
ax1.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 13) #-np.log10(2.43)
ax1.set_ylabel(r'$\mathrm{log_{10}{(R_{c}/kpc})}$', size = 13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
equation = 'y = ' + str(0.56) + 'x' ' + ' + str(round(best_c,2))
at = AnchoredText(equation, frameon=True,
                  loc='upper right', prop=dict(size=8))

at.patch.set_boxstyle("round,pad=0.,rounding_size=0.1")
ax1.add_artist(at)
plt.legend(prop={'size': 11}) #$\mathrm{R_{c}}$
plt.title('95 passive VANDELS galaxies with 3D-HST half-light radii \n (1.0 < z < 1.6)', size = 13)
#ax1.set_xlim(9.75, 11.4)
plt.grid()
plt.savefig('SHENlocal+my_normalisation_fit_z_R_log10M_3DHST.pdf')
plt.close()
#plt.show()


col1 = fits.Column(name='IDs', format='30A', array=new_IDs)
col2 = fits.Column(name='redshifts', format='E', array=new_redshifts)
col3 = fits.Column(name='age', format='E',  array=all_ages)
col4 = fits.Column(name='log10(M*/Msun)', format='E', array=all_masses)
col5 = fits.Column(name='R_c_size', format='E', array=R_c)
col6 = fits.Column(name='Rc_errs', format='E', array=Rc_errs)

hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
#file =  "R_c_v_mass_passive_cat.fits"
#hdu.writeto(file)
cat = Table.read("R_c_v_mass_passive_cat.fits").to_pandas()
R_c_size = cat['R_c_size']
masses = cat['log10(M*/Msun)']

size_mass = np.transpose(np.array([all_masses, np.log10(R_c_size)]))
#print(size_mass)
y_model2_x = np.transpose(np.array([x_values1, y_model2]))
#print(y_model2_x)

print('MAX:', max(y_model2_x[:,0]), max(y_model2_x[:,1]), 'MIN:', min(y_model2_x[:,0]), min(y_model2_x[:,1]))


stack_obs_above = []
stack_obs_below=[]
#max = [11.4,0.5839999999999863]
max = [max(y_model2_x[:,0]), max(y_model2_x[:,1])]
#min= [9.7 ,-0.36800000000001454]
min = [min(y_model2_x[:,0]), min(y_model2_x[:,1])]
is_above = lambda p,a,b: np.cross(p-a, b-a) < 0
a = np.array(min)
b = np.array(max)
p1 = np.array([1.1130000e+01 ,2.8398392e-01])
p2 = np.array([1.1010000e+01, 5.6337827e-01])
p3 = np.array([1.1070000e+01, 2.4578544e-01])
p4 =np.array([1.1140000e+01 , 4.9349055e-01])
p5 = np.array([1.1200000e+01,  6.3046807e-01])
p6=np.array([1.1100000e+01 , 5.2439910e-01])
p7 = np.array([1.1130000e+01, 2.3757902e-01])
p8=np.array([1.1070000e+01 ,8.0749393e-01])
p9 = np.array([1.1080000e+01, 7.2835326e-01])
p10=np.array([1.1270000e+01, 4.4071025e-01])

(fig, ax) = plt.subplots()
model_points = np.array([a,b]) # Add points: (1,2) , (3,5)
model_points_x = model_points[:,0] # For every point, get 1st value, which is x.
model_points_y = model_points[:,1] # For every point, get 2nd value, which is y.
ax.plot(model_points_x, model_points_y, marker="o", color="k")
ax.scatter(all_masses, np.log10(R_c))
if is_above(p1,a,b):
    ax.scatter(p1[0], p1[1], color='green')
else:
    ax.scatter(p1[0], p1[1], color='red')

if is_above(p2,a,b):
    ax.scatter(p2[0], p2[1], color='green')
else:
    ax.scatter(p2[0], p2[1], color='red')

if is_above(p3,a,b):
    ax.scatter(p3[0], p3[1], color='green')
else:
    ax.scatter(p3[0], p3[1], color='red')

if is_above(p4,a,b):
    ax.scatter(p4[0], p4[1], color='green')
else:
    ax.scatter(p4[0], p4[1], color='red')

if is_above(p5,a,b):
    ax.scatter(p5[0], p5[1], color='green')
else:
    ax.scatter(p5[0], p5[1], color='red')

# Point 2:
if is_above(p6,a,b):
    ax.scatter(p6[0], p6[1], color='green')
else:
    ax.scatter(p6[0], p6[1], color='red')
if is_above(p7,a,b):
    ax.scatter(p7[0], p7[1], color='green')
else:
    ax.scatter(p7[0], p7[1], color='red')

# Point 2:
if is_above(p8,a,b):
    ax.scatter(p8[0], p8[1], color='green')
else:
    ax.scatter(p8[0], p8[1], color='red')
if is_above(p9,a,b):
    ax.scatter(p9[0], p9[1], color='green')
else:
    ax.scatter(p9[0], p9[1], color='red')

# Point 2:10
if is_above(p10,a,b):
    ax.scatter(p10[0], p10[1], color='green')
else:
    ax.scatter(p10[0], p10[1], color='red')

plt.savefig('is-point-above-below-line.png')
plt.close()
#input()
#1.1130000e+01  2.8398392e-01
#1.1010000e+01  5.6337827e-01
#1.1070000e+01  2.4578544e-01
#1.1140000e+01  4.9349055e-01
#1.1200000e+01  6.3046807e-01
#1.1100000e+01  5.2439910e-01
#1.1130000e+01  2.3757902e-01
#1.1070000e+01  8.0749393e-01
#1.1080000e+01  7.2835326e-01
#1.1270000e+01  4.4071025e-01

#cross_product = (y_model2_x[:,0]*size_mass[:,1]) - (y_model2_x[:,1]*size_mass[:,0])
cross_product = np.cross(np.array(size_mass), np.array(y_model2_x))


(fig, ax) = plt.subplots()
model_points = np.array([a,b]) # Add points: (1,2) , (3,5)
model_points_x = model_points[:,0] # For every point, get 1st value, which is x.
model_points_y = model_points[:,1] # For every point, get 2nd value, which is y.
ax.plot(model_points_x, model_points_y, marker="o", color="k")
ax.scatter(all_masses, np.log10(R_c))

#for i in size_mass[:,0]:
    #for j in size_mass[:,1]:
points_above = []
points_below = []
df0 = pd.DataFrame(cat)
size = df0['R_c_size']
index = df0.set_index('IDs')
mass = df0['log10(M*/Msun)']
for i in size_mass:
    #print('printing i in size_mass', i)
    if is_above(i, a, b):
        ax.scatter(i[0], i[1], color='green')

        #print('coords above', i)
        points_above.append((i[0],10**i[1]))
        #points_above.append(cat.loc[(mass == i[0])&(size ==10**i[1]), 'IDs'].to_list())
    else:
        ax.scatter(i[0],i[1], color='red')
        #points_below.append(cat.loc[(mass == i[0])&(size ==10**i[1]), 'IDs'].to_list())
        points_below.append((i[0],10**i[1]))
#plt.scatter(masses, np.log10(R_c_size), )
plt.plot(x_values1, y_model2)
plt.grid()
plt.savefig('test_all_points.png')
plt.close()

#print(points_above ,'\n', points_below)
Idnumber = cat['IDs']

#mass_mask0 = (masses > 11.00)
df0 = pd.DataFrame(cat)
size = np.log10(df0['R_c_size'])
#print(size)
index = df0.set_index('IDs')
mass = df0['log10(M*/Msun)']
ID_above1 = []
ID_below1 =[]
ID_above2 = []
ID_below2 =[]
ID_above3 = []
ID_below3 =[]
ID_above4 = []
ID_below4 =[]
ID_above = []
ID_below =[]
for i in points_above:
    if (i[0] > 9.7) and (i[0]<=12.0):
        #print(i)
        massi_ab = np.float(i[0])
        size_ab  = np.log10(i[1])
        ID_above.append(df0.loc[((mass == massi_ab) &(size ==size_ab)) , 'IDs'].str.decode('utf-8').values.tolist())

for j in points_below:
    if (j[0] > 9.7 ) and (j[0]<=12.0):
        #print(j)
        massi_be = np.float(j[0])
        size_be  = np.log10(j[1])
        ID_below.append(df0.loc[((mass==massi_be)&(size==size_be) ), 'IDs'].str.decode('utf-8').values.tolist())#&(masses==mass)








#& (masses==mass)
for i in points_above:
    if (i[0] > 11.0) and (i[0]<=12.0):
        #print(i)
        massi_ab4 = np.float(i[0])
        size_ab4  = np.log10(i[1])
        ID_above4.append(df0.loc[((mass == massi_ab4) &(size ==size_ab4)) , 'IDs'].str.decode('utf-8').values.tolist())
    elif (i[0] > 10.75) and (i[0]<=11.0):
        massi_ab3 = np.float(i[0])
        size_ab3 = np.log10(i[1])
        ID_above3.append(df0.loc[((mass == massi_ab3) &(size ==size_ab3)) , 'IDs'].str.decode('utf-8').values.tolist())
    elif (i[0] > 10.5) and (i[0]<=10.75):
        massi_ab2 = np.float(i[0])
        size_ab2 = np.log10(i[1])
        ID_above2.append(df0.loc[((mass == massi_ab2) &(size ==size_ab2)) , 'IDs'].str.decode('utf-8').values.tolist())
    elif (i[0] > 9.7) and (i[0]<=10.5):
        massi_ab1 = np.float(i[0])
        size_ab1 = np.log10(i[1])
        ID_above1.append(df0.loc[((mass == massi_ab1) &(size ==size_ab1)) , 'IDs'].str.decode('utf-8').values.tolist())




for j in points_below:
    if (j[0] > 11.0 ) and (j[0]<=12.0):
        #print(j)
        massi_be4 = np.float(j[0])
        size_be4  = np.log10(j[1])
        ID_below4.append(df0.loc[((mass==massi_be4)&(size==size_be4) ), 'IDs'].str.decode('utf-8').values.tolist())#&(masses==mass)
    elif (j[0] > 10.75) and (j[0]<=11.0):
        massi_be3 = np.float(j[0])
        size_be3 = np.log10(j[1])
        ID_below3.append(df0.loc[((mass == massi_be3) &(size ==size_be3)) , 'IDs'].str.decode('utf-8').values.tolist())
    elif (j[0] > 10.5) and (j[0]<=10.75):
        massi_be2 = np.float(j[0])
        size_be2 = np.log10(j[1])
        ID_below2.append(df0.loc[((mass == massi_be2) &(size ==size_be2)) , 'IDs'].str.decode('utf-8').values.tolist())
    elif (j[0] > 9.7) and (j[0]<=10.5):
        massi_be1 = np.float(j[0])
        size_be1 = np.log10(j[1])
        ID_below1.append(df0.loc[((mass == massi_be1) &(size ==size_be1)) , 'IDs'].str.decode('utf-8').values.tolist())


#print('above',ID_above)
#print('below', ID_below)
from functools import reduce
import operator

r = ID_above
above_flat = reduce(operator.concat, r)
t = ID_below
below_flat = reduce(operator.concat, t)

print(len(above_flat))
print(len(below_flat))

l = ID_above4
above_flat4 = reduce(operator.concat, l)
k = ID_below4
below_flat4 = reduce(operator.concat, k)

m = ID_above3
above_flat3 = reduce(operator.concat, m)
n = ID_below3
below_flat3 = reduce(operator.concat, n)

q = ID_above2
above_flat2 = reduce(operator.concat, q)
s = ID_below2
below_flat2 = reduce(operator.concat, s)

c = ID_above1
above_flat1 = reduce(operator.concat, c)
d = ID_below1
below_flat1 = reduce(operator.concat, d)


print(len(above_flat4))
print(len(below_flat4))
#if type(i) == i else i
"""
mass_mask0 = (masses > 11.00)
print(sum(mass_mask0))

mask0 = (cross_product < 0)
print(sum(mask0))

cat_mask0 = pd.DataFrame(cat[mass_mask0])
circ_radii0 = cat_mask0['R_c_size']
print('lencircradii', len(circ_radii0))
ind0 = cat[mass_mask0].set_index(cat_mask0['IDs'].str.decode("utf-8").str.rstrip())

index0 = np.log10(circ_radii0[mask0]).index.to_list()

mask20 = (cross_product > 0)
index20 = np.log10(circ_radii0[mask20][mass_mask0]).index.to_list()

#print(len(cat_mask))
for i0 in index0:
    stack_obs_above.append(cat_mask0['IDs'].str.decode("utf-8").str.rstrip()[i0])
print(stack_obs_above)

for id0 in index20:
    stack_obs_below.append(cat_mask0['IDs'].str.decode("utf-8").str.rstrip()[id0])
print(stack_obs_below)

len_below = len(stack_obs_below)
len_above = len(stack_obs_above)

print(len_below, len_above)
#for ID in stack_obs_above:
#    z = ID_.loc[ID, 'zspec']

#below = 32 objects being stacked
#above = 28 objects being stacked

input()

"""

stacking_above = stacks(above_flat)
new_wavs = np.arange(2400, 4200, 1.25)

def plot_stackssingle(stack, name, color):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack*10**18, color=color, lw=1.1 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=17)
    plt.xlim(2350, 4250)
    #plt.ylim(0 ,2.0)
    plt.title('median stacks '+ str(name) +' line \n 9.7 < (log(M*) <= 12.0)', size =18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacking_'+str(name)+'_relation_aboveall.pdf')
    plt.close()
    #plt.show()

plot_stackssingle(stacking_above, 'above', 'r')


stacking_below = stacks(below_flat)

plot_stackssingle(stacking_below,'below', 'k')


def plot_stacks(stack1, stack2):
    plt.figure(figsize=(20,8))
    plt.plot(new_wavs, stack1*10**18, color="r", lw=1.3, ls ='-', label = f'above relation (N = {len(ID_above)})')
    plt.plot(new_wavs, stack2*10**18, color="k", lw=1.4, label = f'below relation (N = {len(ID_below)})')
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2350, 4240)
    #plt.ylim(0. ,1.75)
    plt.legend(fontsize=14)
    plt.title('Median stacks above and below normalised Shen et al. 2003 ETG relation', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('stacks_abovebelow_shennorm_aboveall.pdf')
    plt.close()
    #plt.show()

plot_stacks(stacking_above, stacking_below)



#def plot_histo():



"""
masses1 = cat['log10(M*)']

mass_mask1 = (masses1 > 10.25) & (masses1 <= 10.5)
mass_mask2 = (masses1 > 10.5) & (masses1 <= 10.75)
mass_mask3 = (masses1 > 10.75) & (masses1 <= 11.00)
mass_mask4 = (masses1 > 11.0)
dataframe = pd.DataFrame(cat)
mask_above1 = cross[mass_mask1] < 0
mask_above2 = cross[mass_mask2] < 0
mask_above3 = cross[mass_mask3] < 0
mask_above4 = cross[mass_mask4] < 0
cat_mask1 = pd.DataFrame(cat[mass_mask1])
cat_mask2 = pd.DataFrame(cat[mass_mask2])
cat_mask3 = pd.DataFrame(cat[mass_mask3])
cat_mask4 = pd.DataFrame(cat[mass_mask4])
circ_radii1 = cat_mask1['R_c_size']
circ_radii2 = cat_mask2['R_c_size']
circ_radii3 = cat_mask3['R_c_size']
circ_radii4 = cat_mask4['R_c_size']
#print(len(circ_radii1))
ind1_w = cat[mass_mask1].set_index(cat_mask1['IDs'].str.decode("utf-8").str.rstrip())
ind2_w = cat[mass_mask2].set_index(cat_mask2['IDs'].str.decode("utf-8").str.rstrip())
ind3_w = cat[mass_mask3].set_index(cat_mask3['IDs'].str.decode("utf-8").str.rstrip())
ind4_w = cat[mass_mask4].set_index(cat_mask4['IDs'].str.decode("utf-8").str.rstrip())

index1 = np.log10(circ_radii1[mask_above1]).index.to_list()
index2 = np.log10(circ_radii2[mask_above2]).index.to_list()
index3 = np.log10(circ_radii3[mask_above3]).index.to_list()
index4 = np.log10(circ_radii4[mask_above4]).index.to_list()

mask_below1 = cross[mass_mask1] > 0
mask_below2 = cross[mass_mask2] > 0
mask_below3 = cross[mass_mask3] > 0
mask_below4 = cross[mass_mask4] > 0

ind1 = np.log10(circ_radii1[mask_below1]).index.to_list()
ind2 = np.log10(circ_radii2[mask_below2]).index.to_list()
ind3 = np.log10(circ_radii3[mask_below3]).index.to_list()
ind4 = np.log10(circ_radii4[mask_below4]).index.to_list()


stack_above1 = []
stack_above2 = []
stack_above3 = []
stack_above4 = []

for i1 in index1:
    stack_above1.append(cat_mask1['IDs'].str.decode("utf-8").str.rstrip()[i1])
for i2 in index2:
    stack_above2.append(cat_mask2['IDs'].str.decode("utf-8").str.rstrip()[i2])
for i3 in index3:
    stack_above3.append(cat_mask3['IDs'].str.decode("utf-8").str.rstrip()[i3])
for i4 in index4:
    stack_above4.append(cat_mask4['IDs'].str.decode("utf-8").str.rstrip()[i4])

stack_below1 = []
stack_below2 = []
stack_below3 = []
stack_below4 = []

for id1 in ind1:
    stack_below1.append(cat_mask1['IDs'].str.decode("utf-8").str.rstrip()[id1])
for id2 in ind2:
    stack_below2.append(cat_mask2['IDs'].str.decode("utf-8").str.rstrip()[id2])
for id3 in ind3:
    stack_below3.append(cat_mask3['IDs'].str.decode("utf-8").str.rstrip()[id3])
for id4 in ind4:
    stack_below4.append(cat_mask4['IDs'].str.decode("utf-8").str.rstrip()[id4])

print(stack_below4)
print(stack_above4)


stack_mass1_below = stacks(stack_below1)
stack_mass1_above = stacks(stack_above1)
stack_mass2_below = stacks(stack_below2)
stack_mass2_above = stacks(stack_above2)
stack_mass3_below = stacks(stack_below3)
stack_mass3_above = stacks(stack_above3)
stack_mass4_below = stacks(stack_below4)
stack_mass4_above = stacks(stack_above4)

def plot_all_mass_bins(stack_mass1above, stack_mass2above, stack_mass3above,stack_mass4above, stack_mass1below, stack_mass2below, stack_mass3below, stack_mass4below ):
    fig = plt.figure(figsize=(100,100))
    fig.add_subplot(111, frameon=True)
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=8)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=8)
    #gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols = 2,nrows=2,sharex=True, sharey=True)
    fig.suptitle('Median stacks above and below normalised Shen et al. 2003 ETG relation', size = 8)# excluding possible AGN (CDFS + UDS)')
    ax1.plot(new_wavs, stack_mass1below*10**18, color="k", lw=0.9, ls ='-', label = f'below line (N = {len(stack_below1)})')
    ax1.plot(new_wavs, stack_mass1above*10**18, color="r", lw=0.9, ls ='-', label = f'above line (N = {len(stack_above1)})')
    ax2.plot(new_wavs, stack_mass2below*10**18, color="k", lw=0.9, ls ='-', label = f'below line (N = {len(stack_below2)})')
    ax2.plot(new_wavs, stack_mass2above*10**18, color="r", lw=0.9, ls ='-', label = f'above line (N = {len(stack_above2)})')
    ax3.plot(new_wavs, stack_mass3below*10**18, color="k", lw=0.9, ls ='-', label = f'below line (N = {len(stack_below3)})')
    ax3.plot(new_wavs, stack_mass3above*10**18, color="r", lw=0.9, ls ='-', label = f'above line (N = {len(stack_above3)})')
    ax4.plot(new_wavs, stack_mass4below*10**18, color="k", lw=0.9, ls ='-', label = f'below line (N = {len(stack_below4)})')
    ax4.plot(new_wavs, stack_mass4above*10**18, color="r", lw=0.9, ls ='-', label = f'above line (N = {len(stack_above4)})')
    axs = ax1, ax2, ax3, ax4
    #for ax in axs.flat:
    #    ax.label_outer()

    for ax in axs:
        ax.legend(fontsize=5, loc='upper left')
        #ax.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=8)
        #ax.set_ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=8)
        #ax.set_xlim(2300, 4250)
        ax.tick_params(labelcolor='k', labelsize='small')
        #ax.set_ylim(0. ,2.5)
    #ax4.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=8)
    #ax3.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=8)
    #ax3.set_ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=8)
    #ax1.set_ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=8)
    ax1.text(.99,.95,'(10.25 < log(M*) <= 10.5)',
        horizontalalignment='right',
        transform=ax1.transAxes, size=6)
    ax2.text(.99,.95,'(10.5 < log(M*) <= 10.75)',
        horizontalalignment='right',
        transform=ax2.transAxes, size=6)
    ax3.text(.99,.95,'(10.75 < log(M*) < 11.0)',
        horizontalalignment='right',
        transform=ax3.transAxes, size=6)
    ax4.text(.99,.95,'(log(M*) >= 11.0)',
        horizontalalignment='right',
        transform=ax4.transAxes, size=6 )
    #plt.legend(fontsize=8)
    #set_shared_ylabel(axs, 'Wavelength ($\mathrm{\AA}$)','Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$')

    plt.savefig('subplot_stacks_massbins.pdf')

    plt.close()

plot_all_mass_bins(stack_mass1_above, stack_mass2_above, stack_mass3_above, stack_mass4_above, stack_mass1_below, stack_mass2_below, stack_mass3_below, stack_mass4_below)

#for ID in list_IDs:
    #new_re.append(ID_.loc[ID, 're'])
    #new_UV.append(IDs.loc[ID, 'UV_colour_50'])
    #new_VJ.append(IDs.loc[ID, 'VJ_colour_50'])
"""
