import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from stacks import stacks

my_sample = Table.read('FirstProjectCatalogs/xmatch_spec_derived237objs.fits').to_pandas()
ross_sample = Table.read('FirstProjectCatalogs/vandels_dr4_passive_flag34.fits').to_pandas()



my_ID = my_sample['ID_1']
ross_ID = ross_sample['ID']


my_RA, my_DEC = my_sample['RA_1'].values, my_sample['DEC_1'].values

ross_RA, ross_DEC = ross_sample['RA'].values, ross_sample['DEC'].values

ID_not = []
#list(set(list1).difference(list2))
missing = list(set(ross_ID).difference(my_ID))

#print(missing)

#SIZES CHECKING

concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df = pd.DataFrame(concat_3dhst)
ID_list = np.array(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
ID_ = df.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))


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

fig, ax1 = plt.subplots(figsize=[12,7])
list_IDs = []
new_UV = []
new_VJ = []
new_re = []

for i in ID_list:
    for j in ID_pipes:
        j = j.decode("utf-8")
        if i == j:
            list_IDs.append(i)
            #new_UV.append(UV[i])
            #new_VJ.append(VJ[i])
            #new_re.append(re[i])
print(len(list_IDs))

for ID in list_IDs:
    new_re.append(ID_.loc[ID, 're'])
    new_UV.append(IDs.loc[ID, 'UV_colour_50'])
    new_VJ.append(IDs.loc[ID, 'VJ_colour_50'])

print(len(new_re), len(new_UV))




"""
x = np.linspace(0., 3.9, 1000)
Y = 0.88*x + 0.69
ax1.plot(x, Y, linestyle = '-', lw = 2., color = 'k')
if [Y <= 1.3]:
    xmax = (1.3 - 0.69)/0.88
    xmin = 0
    ax1.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 2.)#1.3

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.69
    Ymax = 0.88*len(x) +0.69

    ax1.vlines(x=1.6, ymin=Ymin, ymax=3.5, ls = '-', lw = 2.)#1.3



im = ax1.scatter(new_VJ, new_UV, s=30, c=new_re, cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'size (re)', size=10)
ax1.set_xlim(0.4, 3.)
ax1.set_ylim(0, 2.5)
ax1.set_title('UVJ diagram for VANDELS objects 1 < z < 2.5', size = 25)
ax1.set_xlabel(r'V-J', size=10)
ax1.set_ylabel(r'U-V',size=10)
ax1.set_xticks(np.arange(0.5, 3., 0.5))
ax1.tick_params(axis = 'both', labelsize =10, size = 10)

plt.show()
"""
