import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import matplotlib
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from x_ray_check import stacks
from astropy.stats import sigma_clip
#from stacks import stacks
plt.rc('text', usetex=True)
vandels_cat_new = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()

UV = vandels_cat_new["UV_colour_50"]
VJ = vandels_cat_new["VJ_colour_50"]
SSFR = vandels_cat_new["ssfr_50"]
objID = vandels_cat_new["#ID"]
stel_mas = vandels_cat_new["stellar_mass_50"]
age = vandels_cat_new["exponential:age_50"]
age1 = vandels_cat_new["mass_weighted_age_50"]

UV_neg_err =  UV - vandels_cat_new["UV_colour_16"]
UV_pos_err = vandels_cat_new["UV_colour_84"] -UV
#print(UV_pos_err[0],UV[0], UV_neg_err[0])

VJ_neg_err = VJ - vandels_cat_new["VJ_colour_16"]
VJ_pos_err = vandels_cat_new["VJ_colour_84"] - VJ

#print(VJ_pos_err[0],VJ[0],VJ_neg_err[0])

#fig, ax = plt.subplots()
xerr = [VJ_neg_err, VJ_pos_err]
c = age1
csfr = SSFR

"""
fig, ax = plt.subplots()

x = np.linspace(0., 3.9, 1000)

Y = 0.88*x + 0.49
Ycut = 0.88*x + 0.69
Y1 = -0.88*x + 3.7
Y2 = -0.88*x + 1.9
Y3 = 0.88*x + 1.1

Ys = -0.88*x + 1.9 + 0.45
Yss = Ys +0.45
Ysss = Yss +0.45
x1  = 3.21/1.76
x2 = 1.41/1.76
x3 = 0.8/1.76
x4 = 2.6/1.76

c1 = 4.19/2
c2 = 1.59/2
c3 = 4.8/2
c4 = 3/2
c_new = 2.195
x_new = 3.01/1.76

#plt.plot(x, Y, linestyle = '-', lw = 1.3, color = 'k')
#plt.plot(x, Ycut, linestyle = '--', lw = 1.3, color = 'k')
plt.plot(x[(Ys<Y3) & (Ys > Ycut)], Ys[(Ys<Y3) & (Ys > Ycut)], linestyle = '-', lw = 1.3, color = 'r')
plt.plot(x[(Yss<Y3) & (Yss > Ycut)], Yss[(Yss<Y3) & (Yss > Ycut)], linestyle = '-', lw = 1.3, color = 'r')
plt.plot(x[(Ysss<Y3) & (Ysss > Ycut)], Ysss[(Ysss<Y3) & (Ysss > Ycut)], linestyle = '-', lw = 1.3, color = 'r')

bin1 = []
bin2 = []
bin3 =[]
bin4 =[]
###### bins using the 0.88x + 0.49 line ########
#m1, n1, m2, n2 = x2, c2, 1.25/1.76, 1.725
#m3,n3, m4, n4= 1.86/1.76, 1.42, 1.7/1.76, 1.95
#m5, n5, m6, n6 = 2.31/1.76, 1.645, 2.15/1.76, 2.175
#m7, n7, m8, n8 = 2.76/1.76, 1.87, 2.6/1.76, 2.4
#bottom left top right corners of each rectangle

### using the stricter cut of 0.88x + 0.69 ################

m1, n1, m2, n2 = 1.21/1.76, 1.295, 1.25/1.76, 1.725
m3,n3, m4, n4= 1.66/1.76, 1.52, 1.7/1.76, 1.95
m5, n5, m6, n6 = 2.11/1.76, 1.745, 2.15/1.76, 2.175
m7, n7, m8, n8 = 2.56/1.76, 1.97, 2.6/1.76, 2.4

for i in range(len(UV)):
    point = Point(VJ[i], UV[i])
    #print(point1) polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
    polygon = Polygon([(m1, n1), (x4, c4), (m2, n2), (m3, n3)])
    polygon2 = Polygon([(m3, n3), (m2, n2), (m4, n4), (m5, n5)])
    polygon3 = Polygon([(m5, n5), (m4, n4), (m6, n6), (m7, n7)])
    polygon4 = Polygon([(m7, n7), (m6, n6), (m8, n8), (x1, c1)])
    if polygon.contains(point):
        bin1.append(ID[i])
    elif polygon2.contains(point):
        bin2.append(ID[i])
    elif polygon3.contains(point):
        bin3.append(ID[i])
    elif polygon4.contains(point):
        bin4.append(ID[i])
    #print(ID[i])
print(len(bin1), len( bin2), len(bin3),len(bin4))


#print( bin2, bin3, bin4)
lim = (Y < c1) & (x > x2)
lim2 = (x > x3) & (Y1 > c_new ) & (Y1 < c3)
#lim3 =  (Y3 < c3 ) & (Y3 > c4) & (x > x4)
lim4 =   (Y2 < c4) & (Y2 > 1.3) & (x < x2)

lim5 =  (Y3 > c4) & (Y3<c3)

if [Y <= 1.3]:
    xmax = (1.3 - 0.49)/0.88
    xmin = 0
    plt.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 1.3)

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.49
    Ymax = 0.88*len(x) +0.49

    plt.vlines(x=1.6, ymin=Ymin, ymax=3.5, ls = '-', lw = 1.3)

plt.plot(x,Y, linestyle = '-', lw = 1.3, color = 'k' )

plt.plot(x[(x <1.6) & (Y > 1.3)], Y[(x <1.6) & (Y > 1.3)], linestyle = '-', lw = 1.3, color = 'k')
#plt.plot(x[(x <1.6) & (Ycut > 1.3)], Ycut[(x<1.6) & (Ycut > 1.3)], ls = '--', lw = 1.5, color='k')
plt.plot(x[(0.682 < x) & (x<x_new) & (Ycut < c_new)], Ycut[(0.682 < x) & (x < x_new) & (Ycut < c_new)], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim2], Y1[lim2], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim4], Y2[lim4], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim5], Y3[lim5], ls = '-', lw = 1.3, color='r')


"""

both_xray = Table.read('FirstProjectCatalogs/concat_possible_xray_matches_massi.fits').to_pandas()
df1 = pd.DataFrame(both_xray)#, index = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8")))
ID_list1 = df1.set_index(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))

agn = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))
print(len(agn))
passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
df2 = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )
ID_list2 = df2.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))

all_obs = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))

missing = list(set(all_obs).difference(agn))
print(len(missing))
new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)
#print(new_df)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))

print('ID_LIST length = ', len(ID_list))

input()

def plot_stacks(new_wavs, med_stack, bin_number):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, med_stack*10**18, color="black", lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    #plt.xlim(2300, 4250)
    plt.ylim(0 ,2.0)
    plt.title('Median Stacked Spectra of galaxies 1 < z < 2.5', size=17)
    plt.savefig('stack_plot_bin'+str(bin_number)+ 'stricter.pdf')
    plt.close()


fig, ax1 = plt.subplots(figsize=[20,15])

# equivalent but more general
ax2 = fig.add_subplot(411, position=[0.36, 0.48, 0.38, 0.11])
ax3 = fig.add_subplot(412, position=[0.36, 0.37, 0.38, 0.11],sharex=ax2, sharey=ax2)
ax4 = fig.add_subplot(413, position=[0.36, 0.26, 0.38, 0.11],sharex=ax2, sharey=ax2)
ax5 = fig.add_subplot(414, position=[0.36, 0.15, 0.38, 0.11],sharex=ax2, sharey=ax2)

x = np.linspace(0., 3.9, 1000)

Y = 0.88*x + 0.69
Y1 = -0.88*x + 3.7
Y2 = -0.88*x + 1.9
Y3 = 0.88*x + 1.1

Ys = -0.88*x + 1.9 + 0.45
Yss = Ys +0.45
Ysss = Yss +0.45
x1  = 3.21/1.76
x2 = 1.41/1.76
x3 = 0.8/1.76
x4 = 2.6/1.76

c1 = 4.19/2
c2 = 1.59/2
c3 = 4.8/2
c4 = 3/2

ax1.plot(x, Y, linestyle = '-', lw = 2., color = 'k')
ax1.plot(x[(Ys<Y3) & (Ys > Y)], Ys[(Ys<Y3) & (Ys > Y)], linestyle = '-', lw = 2., color = 'r')
ax1.plot(x[(Yss<Y3) & (Yss > Y)], Yss[(Yss<Y3) & (Yss > Y)], linestyle = '-', lw = 2., color = 'r')
ax1.plot(x[(Ysss<Y3) & (Ysss > Y)], Ysss[(Ysss<Y3) & (Ysss > Y)], linestyle = '-', lw = 2., color = 'r')

bin1 = []
bin2 = []
bin3 =[]
bin4 =[]

### using the stricter cut of 0.88x + 0.69 ################

m1, n1, m2, n2 = 1.21/1.76, 1.295, 1.25/1.76, 1.725
m3,n3, m4, n4= 1.66/1.76, 1.52, 1.7/1.76, 1.95
m5, n5, m6, n6 = 2.11/1.76, 1.745, 2.15/1.76, 2.175
m7, n7, m8, n8 = 2.56/1.76, 1.97, 2.6/1.76, 2.4

#bottom left top right corners of each rectangle
#m1, n1, m2, n2 = x2, c2, 1.25/1.76, 1.725
#m3,n3, m4, n4= 1.86/1.76, 1.42, 1.7/1.76, 1.95
#m5, n5, m6, n6 = 2.31/1.76, 1.645, 2.15/1.76, 2.175
#m7, n7, m8, n8 = 2.76/1.76, 1.87, 2.6/1.76, 2.4
q,r = 3.01/1.76, 2.195

for i in range(len(UV)):
    point = Point(VJ[i], UV[i])
    #print(point1) polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
    polygon = Polygon([(m1, n1), (x3, c4), (m2, n2), (m3, n3)])
    polygon2 = Polygon([(m3, n3), (m2, n2), (m4, n4), (m5, n5)])
    polygon3 = Polygon([(m5, n5), (m4, n4), (m6, n6), (m7, n7)])
    polygon4 = Polygon([(m7, n7), (m6, n6), (m8, n8), (q,r)])
    if polygon.contains(point):
        bin1.append(objID[i])
    elif polygon2.contains(point):
        bin2.append(objID[i])
    elif polygon3.contains(point):
        bin3.append(objID[i])
    elif polygon4.contains(point):
        bin4.append(objID[i])
    #print(ID[i])
print(len(bin1), len(bin2), len(bin3),len(bin4))
#print(bin1)

new_ID1 = []
new_ID2 = []
new_ID3 = []
new_ID4 = []

for s in bin1:
    new_ID1.append(s.decode("utf-8").rstrip())
for j in bin2:
    new_ID2.append(j.decode("utf-8").rstrip())
for k in bin3:
    new_ID3.append(k.decode("utf-8").rstrip())
for l in bin4:
    new_ID4.append(l.decode("utf-8").rstrip())

#print(new_ID1, new_ID2, new_ID3, new_ID4)

#all_IDs_strict_UVJ = new_ID1 + new_ID2 +new_ID3 + new_ID4


c1 = 2.195
x1 = 3.01/1.76
#print( bin2, bin3, bin4)
lim = (Y < r) & (x > x2)
lim2 = (x > x3) & (Y1 > r ) & (Y1 < c3)
#lim3 =  (Y3 < c3 ) & (Y3 > c4) & (x > x4)
lim4 =   (Y2 < c4) & (Y2 > 0.682) & (x < m1)

lim5 =  (Y3 > c4) & (Y3<c3)

if [Y <= 1.3]:
    xmax = (1.3 - 0.69)/0.88
    xmin = 0
    ax1.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 2.)#1.3

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.69
    Ymax = 0.88*len(x) +0.69

    ax1.vlines(x=1.6, ymin=Ymin, ymax=3.5, ls = '-', lw = 2.)#1.3

#plt.plot(x[lims1], Y[lims1], ls = '-', lw = 0.9, color='r')
ax1.plot(x[lim], Y[lim], ls = '-', lw = 2., color='r',)
ax1.plot(x[lim2], Y1[lim2], ls = '-', lw = 2., color='r')
ax1.plot(x[lim4], Y2[lim4], ls = '-', lw = 2., color='r')
ax1.plot(x[lim5], Y3[lim5], ls = '-', lw = 2., color='r')


im = ax1.scatter(VJ, UV, s=150, c=age1, cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5 )
cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar.set_label(r'Age (Gyr)', size=25)
ax1.set_xlim(0.4, 3.)
ax1.set_ylim(0, 2.5)
ax1.set_title('UVJ diagram for VANDELS objects 1 $<$ z $<$ 2.5', size = 25)
ax1.set_xlabel(r'V - J', size=25)
ax1.set_ylabel(r'U - V',size=25)
ax1.set_xticks(np.arange(0.5, 3., 0.5))
ax1.tick_params(axis = 'both', labelsize =20, size = 20)

wavs_stack = np.arange(2400, 4200, 1.25)


#median ages going up the passive box in bins
df1 = pd.DataFrame(vandels_cat_new)
#df2 = pd.DataFrame(passive_cut)

concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df3 = pd.DataFrame(concat_3dhst)
size_list = np.array(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))

ID_ = df3.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))

ID_pipes = vandels_cat_new['#ID'].values
IDs = df1.set_index(s.decode('utf-8') for s in vandels_cat_new['#ID'])
#IDs = df1.set_index(s.decode('utf-8') for s in vandels_cat_new['#ID'])
all_ages = []
all_masses = []
all_sizes = []
list_IDs = []
for i in size_list:
    for j in ID_pipes:
        j = j.decode("utf-8")
        if i == j:
            list_IDs.append(i)

for ID in list_IDs:
    if ID not in agn:
        all_ages.append(IDs.loc[ID, "mass_weighted_age_50"])
        #print(all_ages)
        all_masses.append(ID_list.loc[ID, "log10(M*)"])
        #print(all_masses)
        all_sizes.append(ID_.loc[ID, 're'])

print(len(all_sizes))
median_age = np.nanmedian(all_ages)
print('median age:', median_age , 'Gyrs')
median_mass = np.nanmedian(all_masses)
print('median mass:', median_mass)
median_size = np.nanmedian(all_sizes)
print('median size:', median_size , 'kpc')

for xray in agn:
    if xray in new_ID1:
        new_ID1.remove(xray)
    if xray in new_ID2:
        new_ID2.remove(xray)
    if xray in new_ID3:
        new_ID3.remove(xray)
    if xray in new_ID4:
        new_ID4.remove(xray)

all_IDs_strict_UVJ = new_ID1 + new_ID2 +new_ID3 + new_ID4

#med_spec_units, colour_index, D4000_index, Mg_UV_index, H_delta_EW, ages, masses, IDplz = stacks(all_IDs_strict_UVJ)

data_indices = {'C(29-33)': colour_index, 'Mg_UV': Mg_UV_index, 'H_delta_EW': H_delta_EW, 'D4000': D4000_index, 'ages': ages, 'masses':masses, 'ID': IDplz}
df_over = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses', 'ID'])
df_less = pd.DataFrame(data_indices, columns = ['C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses', 'ID'])
df_less = df_less.groupby((df_less['D4000'] <= 1.35)).get_group(True)

df_over = df_over.groupby((df_over['D4000'] > 1.35)).get_group(True)
#print(df_over)

ID_less_than_135 = df_less['ID'].values
ID_more_than_135 = df_over['ID'].values

#print('IDS for less and more D4000 = 1.35', ID_less_than_135, ID_more_than_135)
#ages = df_over['ages']
#D4000_above = df_over['D4000']
#H_delta_EW_above = df_over['H_delta_EW']
#ages_above = df_over['ages']




#print(len(colour_index),len(D4000_index))
"""
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
figs1, axs1 = plt.subplots(figsize=[12,8.5])
ims1 = axs1.scatter(colour_index, Mg_UV_index, s=160, c=ages, cmap=plt.cm.magma, marker='o', linewidth=0.5 , alpha = 0.9)
#ax1.errorbar(colour_index, Mg_UV_index ,marker='s', xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
cbars = figs1.colorbar(ims1, ax=axs1)
cbars.set_label(r'Age (Gyr)', size=12)
#ax1.scatter(colour_index, Mg_UV_index ,marker='o',color = 'k',zorder = 1)
axs1.set_ylim(0.78, 2.3)
axs1.set_xlim(0.2, 1.2)
axs1.set_xlabel("C(29-33)", size=13)
axs1.set_ylabel("MgUV", size=13)
#plt.title('', size =14)# excluding possible AGN (CDFS + UDS)')
#plt.savefig('stricterUVJ_mgUVvC29_33_cbar.pdf')
plt.close()
"""

IDs_all_stacks = ID_list.index.values
IDs_not_strict = []
for ID in IDs_all_stacks:
    if ID not in all_IDs_strict_UVJ:
        IDs_not_strict.append(ID)


#print(len(IDs_not_strict))

#input()

med_spec_units_not, colour_index_not, D4000_index_not, Mg_UV_index_not, H_delta_EW_not, ages_not, masses_not, Idplzwhy = stacks(IDs_not_strict)
data_indices_not = {'C(29-33)': colour_index_not, 'Mg_UV': Mg_UV_index_not, 'H_delta_EW': H_delta_EW_not, 'D4000': D4000_index_not, 'ages': ages_not, 'masses':masses_not, "ID":Idplzwhy}
df_new1 = pd.DataFrame(data_indices_not, columns = [ 'C(29-33)', 'Mg_UV', 'H_delta_EW', 'D4000', 'ages', 'masses', 'ID'])
d4000_not = df_new1['D4000']
hdelta_not = df_new1['H_delta_EW']
ages_not = df_new1['ages']
masses_not = df_new1['masses']


"""
figs2, axs2 = plt.subplots(figsize=[12,8.5])
ims2 = axs2.scatter(D4000_above, H_delta_EW_above, s=200, c = ages_above, cmap=plt.cm.magma, marker = 'o',linewidth=1.0, alpha = 0.9, label = 'Objects above 0.88(V - J) + 0.69 cut', zorder = 0)

#ims2 = axs2.scatter(d4000_not, hdelta_not , s=200, c=ages_not, cmap=plt.cm.magma, marker='*',linewidth=0.5, alpha = 0.9, label = 'Objects below 0.88(V - J) + 0.69 cut' , zorder =1, edgecolors='k')
#ax1.errorbar(colour_index, Mg_UV_index ,marker='s', xerr = 0.1*np.ones(len(colour_index))*colour_index, yerr = 0.1*np.ones(len(Mg_UV_index))*Mg_UV_index, ecolor='k', color ='k', linestyle =" ", elinewidth = 0.4, capsize = 4, mew = 0.3, zorder = 1)
cbars2 = figs2.colorbar(ims2, ax=axs2)
cbars2.set_label(r'Age (Gyr)', size=12)
#df_new = df_new.groupby(df_new['D4000']<1.34).get_group(True)
#df_new = df_new.groupby(df_new['H_delta_EW']>=4.5).get_group(True)
#D4000_below = df_new['D4000']
#H_delta_EW_below = df_new['H_delta_EW']
#ages_below= df_new['ages']
#axs2.scatter(D4000_above, H_delta_EW_above, s=200, c = ages_above, cmap=plt.cm.magma, marker = 'o',linewidth=1.0, alpha = 0.9, label = 'Objects above 0.88(V - J) + 0.69 cut', zorder = 0)
#ax2.scatter(D4000_index, H_delta_EW, marker='o',color = 'k',zorder = 1)
axs2.set_xlabel("D$_n$4000", size=13)
axs2.set_ylabel("EW(H$\delta$)", size=13)
axs2.legend()
#plt.savefig('onlystrictstars_hdeltad4000test.pdf')
plt.close()
"""


print('plot done')
input()

age_bin1 = []
for ID in new_ID1:
    age_bin1.append(IDs.loc[ID, "mass_weighted_age_50"])
age_bin2 = []
for ID2 in new_ID2:
    age_bin2.append(IDs.loc[ID2, "mass_weighted_age_50"])
age_bin3 = []
for ID3 in new_ID3:
    age_bin3.append(IDs.loc[ID3, "mass_weighted_age_50"])
age_bin4 = []
for ID4 in new_ID4:
    age_bin4.append(IDs.loc[ID4, "mass_weighted_age_50"])

#print(np.median(age_bin1), np.median(age_bin2), np.median(age_bin3), np.median(age_bin4))

#print(len(new_ID1), len(new_ID2), len(new_ID3), len(new_ID4))



med_stack1, colour_index1, D4000_index1, Mg_UV_index1, H_delta_EW1, ages1, masses1, IDplz1 = stacks(new_ID1)#, ID_list)
#plot_stacks(new_waves, med_stack1, 9)

med_stack2, colour_index2, D4000_index2, Mg_UV_index2, H_delta_EW2, ages2, masses2, IDplz1 = stacks(new_ID2)#, ID_list)
#plot_stacks(new_waves, med_stack2, 10)

med_stack3, colour_index3, D4000_index3, Mg_UV_index3, H_delta_EW3, ages3, masses3, IDplz1 = stacks(new_ID3)#, ID_list)
#plot_stacks(new_waves, med_stack3, 11)

med_stack4, colour_index4, D4000_index4, Mg_UV_index4, H_delta_EW4, ages4, masses4, IDplz1 = stacks(new_ID4)#, ID_list)

#print(new_ID4)
#figs, axs = plt.subplots(figsize=[12,8.5])
#axs.plot(med_stack1*10**18)


ax2.plot(wavs_stack, med_stack1*10**18, 'k',lw=0.8, label = f'N = {len(new_ID1)} ')
ax3.plot(wavs_stack, med_stack2*10**18, 'k',lw=0.8, label = f'N = {len(new_ID2)} ')
ax4.plot(wavs_stack, med_stack3*10**18, 'k',lw=0.8, label = f'N = {len(new_ID3)} ')
ax5.plot(wavs_stack, med_stack4*10**18, 'k',lw=0.8, label = f'N = {len(new_ID4)} ')

ax5.set_xlabel('Wavelength ($\mathrm{\AA}$)', size=15)
ax2.tick_params(axis = 'y', labelsize = 10, size=6)
ax3.tick_params(axis = 'y', labelsize = 10, size=6)
ax4.tick_params(axis = 'y', labelsize = 10, size=6)
ax5.tick_params(axis = 'both', labelsize = 10, size=5)
ax2.set_xlim(2380,4230)
ax3.set_xlim(2380,4230)
ax4.set_xlim(2380,4230)
ax5.set_xlim(2380,4230)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)
ax2.legend(loc = 'upper left')
ax3.legend(loc = 'upper left')
ax4.legend(loc = 'upper left')
ax5.legend(loc = 'upper left')
#ax2.set_ylabel('Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$', size=10)
ax2.set_title('Median Stacked Spectra of galaxies 1 $<$ z $<$ 2.5', size=15)

#ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#ax2 = fig.add_axes([0.72, 0.72, 0.16, 0.16])

#plt.savefig('UVJ_with_stacked_spectra_stricter_cut_legendtest.pdf')
plt.close()

print('next')
input()
cat2 = Table.read("Re_cat.fits").to_pandas()
df2 = pd.DataFrame(cat2)
df2 = df2.groupby(df2['log10(M*/Msun)']>10.4).get_group(True)
IDs = all_IDs_strict_UVJ

R_e = df2["Re_kpc"]
R_e_errs = df2["Re_kpc_errs"]
mass = df2["log10(M*/Msun)"]
ssfr_50 = df2["SSFR"]
size  = np.array(np.log10(R_e))
#IDs_me = df2.set_index(cat2["IDs"])

concat_3dhst = Table.read('FirstProjectCatalogs/concat_3dhst_passive_match.fits').to_pandas()
df = pd.DataFrame(concat_3dhst)
ID_list1 = np.array(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
ID_ = df.set_index(concat_3dhst['FIELD'].str.decode("utf-8").str.rstrip() + concat_3dhst['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + concat_3dhst['CAT'].str.decode("utf-8"))
agn = np.array(both_xray['FIELD'].str.decode("utf-8").str.rstrip()+ both_xray['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + both_xray['CAT'].str.decode("utf-8"))
ID_pipes = vandels_cat_new['#ID'].values
IDs_pipes = df1.set_index(s.decode('utf-8') for s in vandels_cat_new['#ID'])

new_IDs = []
new_redshifts = []
q_ratio = []
ssfr2 = []
ages = []
masses = []
re = []
re_errs =[]

list_IDs = []
for i in ID_list1:
    for j in IDs:
        if i == j:
            list_IDs.append(i)
print(len(list_IDs))
hdelta = []
d4000 = []
c29_33 = []
mgUV = []
RA = []
DEC = []

for ID in list_IDs:
    if ID not in agn:
        new_redshifts.append(ID_.loc[ID, "zspec"])
        ssfr2.append(IDs_pipes.loc[ID, "ssfr_50"])
        ages.append(IDs_pipes.loc[ID, "mass_weighted_age_50"])
        masses.append(ID_list.loc[ID, "log10(M*)"])
        #print(all_masses)
        RA.append(ID_list.loc[ID, 'RA_1'])
        DEC.append(ID_list.loc[ID, 'DEC_1'])
        re.append(ID_.loc[ID, 're'])
        re_errs.append(ID_.loc[ID, 'dre'])
        q_ratio.append(ID_.loc[ID, 'q'])
        new_IDs.append(ID)
#objects_list = list(set(np.array(df2["IDs"].values)).intersection(IDs))
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import scipy

cosmo  = FlatLambdaCDM(H0=70, Om0=0.3)
Mpc_to_kpc = 1000
masses = np.array(masses)
arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(np.array(new_redshifts))
print('max=', max(re))

Re_kpc = (np.array(re)*u.arcsec)/arcsec_per_kpc
Re_kpc_errs = (np.array(re_errs)*u.arcsec)/arcsec_per_kpc
from astropy.table import Table
from astropy.io import fits

Rc = (np.sqrt(np.array(q_ratio))*Re_kpc) /u.kpc
Rc_errs = (np.sqrt(np.array(q_ratio))*Re_kpc_errs) /u.kpc

col1 = fits.Column(name='IDs', format='30A', array=new_IDs)
col8 = fits.Column(name='RA', format = 'E', array = RA)
col9 = fits.Column(name='DEC', format = 'E', array = DEC)
col2 = fits.Column(name='redshifts', format='E', array=new_redshifts)
col3 = fits.Column(name='age', format='E',  array=ages)
col4 = fits.Column(name='log10(M*/Msun)', format='E', array=masses)
col5 = fits.Column(name='Re_kpc', format='E', array=Re_kpc/u.kpc)
col6 = fits.Column(name='Re_kpc_errs', format='E', array=Re_kpc_errs/u.kpc)
col7 = fits.Column(name ='SSFR', format = 'E', array = ssfr2)

hdu = fits.BinTableHDU.from_columns([col1, col8, col9, col2, col3, col4, col5, col6, col7])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
#file =  "Re_cat_strict_UVJ_cut.fits"
#hdu.writeto(file)
# only 68 objects in the strict passive box where all objects are above o.88x +0.69, x<1.6, y>1.3

cat3 = Table.read("Re_cat_strict_UVJ_cut.fits").to_pandas()
df3 = pd.DataFrame(cat3)
df3 = df3.groupby(df3['log10(M*/Msun)']>10.4).get_group(True)
strict_masses = df3["log10(M*/Msun)"]
strict_sizes= df3["Re_kpc"]
R_e_errs = df3["Re_kpc_errs"]
redshifts = df3['redshifts']
age=df3['age']
log_A_model = np.arange(-0.3, 0.3, 0.01)
best_chi2 = np.inf
masses = np.array(strict_masses.values)
def vdw_relation(logA, alpha, x_values):
    logR_eff = logA + alpha*(np.log10((10**x_values)/(5*10**10)))
    return logR_eff

for cvals2 in range(len(log_A_model)):
    c_vals2 =log_A_model[cvals2]

    vdw_model = vdw_relation(c_vals2, 0.76, strict_masses)
    diffs2 = vdw_model - np.log10(strict_sizes)
    #print(diffs)
    #print(y_model, '\n', np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq2 = np.sum((diffs2**2)/((4.34*(R_e_errs/strict_sizes))))

    if chisq2 < best_chi2:
        best_chi_vdw = chisq2
        best_c_vdw = c_vals2

print(f'best_c_vdw: {best_c_vdw} \n best_chi_vdw: {best_chi_vdw}')

"""
c__shen = np.arange(-8.0, -3.0, 0.01)
model_shen_mask = (masses > 10.4)
best_chi  = np.inf

errs = 0.434*(Rc_errs/Rc)

for cvals in range(len(c_model)):
    c_vals = c_model[cvals]

    y_model = 0.56*masses + c_vals
    diffs = y_model - np.log10(Rc)
    #print(diffs)
    #print(y_model, np.log10(R_c))#(0.434*(Rc_errs/R_c)
    chisq = np.sum((diffs**2)/((10*errs)**2))

    if chisq < best_chi:
        best_c = c_vals
        best_chi = chisq
#print(vdw_model)
"""

redshifts = df3['redshifts']
alpha = 0.76
log_A = 0.22
x = np.linspace(9.5, 11.5, len(strict_sizes))
log_Reff = log_A + alpha*np.log10((10**x)/(5*10**10))

IDs = cat3['IDs'].values
vdw_norm_model = vdw_relation(best_c_vdw, 0.76, x)

print(IDs)



input()
#y_model = 0.56*x + (best_c)
index = np.log10(strict_sizes).index.to_list()
#print(index)
IDs = cat3['IDs']
#IDs = d4000lessIDs
size = np.array(np.log10(strict_sizes))

vdw_norm_model = np.array(vdw_norm_model)
mask = (size>vdw_norm_model)
mask1 = (size<vdw_norm_model)
index_masked = np.log10(strict_sizes)[mask].index.to_list()
index_masked2 = np.log10(strict_sizes)[mask1].index.to_list()

IDs_above = IDs[index_masked].str.decode("utf-8").str.rstrip().values
IDs_below= IDs[index_masked2].str.decode("utf-8").str.rstrip().values
all_IDs = np.concatenate((IDs_above, IDs_below), axis = None)

new_wavs = np.arange(2400, 4200, 1.25)
#med_stacks, colour_index, D4000_index, Mg_UV_index, H_delta_EW, ages = stacks(all_IDs)
#shen_model = 0.56*x + (best_c)

fig2, ax2 = plt.subplots(figsize=[12,8.5])
im2 = ax2.scatter(strict_masses, np.log10(strict_sizes), s=200, c=age, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.8)
cbar2 = fig2.colorbar(im2, ax=ax2)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbar2.set_label(r'Age (Gyr)', size=12)
#ax1.plot(x, y_model, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
ax2.plot(x, vdw_norm_model, 'k',ls = '--', alpha = 0.6, lw = 1.2, label=f'A. van der Wel., 2014 z = 1.25 ETG relation \n normalised by f = {round((best_c_vdw-0.22),2)}')
#ax1.plot(x, log_Reff, 'r' ,lw = 1.2, label= 'v. der Wel z = 1.25 ETG relation (not normalised)')
#ax1.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
ax2.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax2.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(prop={'size': 10})
plt.title(r'3D-HST $\mathrm{log_{10}{(R_{e}/kpc)}} \ $ v $ \ \mathrm{log_{10}{(M*/M_{\odot})}}$ (N = 55)', size = 13)
plt.xlim(10.3, 11.4)
plt.ylim(-0.4,1.1)
#plt.savefig('RevM_cbar_ages_strictUVJ.pdf')
plt.close()


print('plots you need today massi')
print('next')
input()

low_lim3, upp_lim3 = 10.5, 10.75
low_lim, upp_lim = 10.75, 11.0
low_lim2, upp_lim2 = 11.0, 11.3
low_lim5, upp_lim5 = 10.4, 11.0
# plotting UVJ strict selected stacks up the UVJ passive box, in mass bins

from stacking_code import stacks
cat3 = Table.read("Re_cat_strict_UVJ_cut.fits").to_pandas()


def stack_lims(lower_lim, higher_lim):
    #mass_mask = (masses>lower_lim) & (masses <= higher_lim)
    #df4['IDs'].str.decode("utf-8").str.rstrip().values == ID_less_than_135
    df4 = pd.DataFrame(cat3)
    df4 = df4.groupby((df4['log10(M*/Msun)']>lower_lim)&(df4['log10(M*/Msun)']<=higher_lim)).get_group(True)
    #df4 = df4.set_index(s.rstrip() for s in d4000lessIDs)
    #R_e = df0["Re_kpc"]
    stricter_masses = df4["log10(M*/Msun)"]
    stricter_sizes= df4["Re_kpc"]
    R_e_errs_4 = df4["Re_kpc_errs"]
    redshifts_4 = df4['redshifts']
    size_4 = np.array(np.log10(stricter_sizes))
    index_4 = np.log10(stricter_sizes).index.to_list()
    x_array_4 = np.linspace(lower_lim, higher_lim, len(size_4))
    vdw_norm_model_4 = best_c_vdw + alpha*np.log10((10**x_array_4)/(5*10**10))
    mask_4 = (size_4>vdw_norm_model_4)
    mask1_4 = (size_4<vdw_norm_model_4)
    index_masked_mass_4 = np.log10(stricter_sizes)[mask_4].index.to_list()
    index_masked2_mass_4 = np.log10(stricter_sizes)[mask1_4].index.to_list()

    IDs_above_4 = IDs[index_masked_mass_4].str.decode("utf-8").str.rstrip().values
    IDs_below_4 = IDs[index_masked2_mass_4].str.decode("utf-8").str.rstrip().values
    d4000lessIDs_above = []

    for i in ID_less_than_135:
        for j in IDs_above_4:
            #j = j.decode("utf-8")
            if i == j:
                d4000lessIDs_above.append(i)

    d4000lessIDs_below= []

    for k in ID_less_than_135:
        for l in IDs_below_4:
            #l = l.decode("utf-8")
            if k == l:
                d4000lessIDs_below.append(k)

    #d4000moreIDs = []
    """
    for k in ID_more_than_135:
        for l in IDs:
            l = l.decode("utf-8")
            if k == l:
                d4000moreIDs.append(k)
    """
    """for id in ID_less_than_135:
        if id in IDs:
            IDs = id"""

    print('please work less above:', d4000lessIDs_above , '\n please work less below:', d4000lessIDs_below)

    stacking_above_4 = stacks(d4000lessIDs_above)
    stacking_below_4 = stacks(d4000lessIDs_below)
    len_above, len_below = len(d4000lessIDs_above), len(d4000lessIDs_below)
    #len_above, len_below = len(IDs_above_4), len(IDs_below_4)
    #stacking_above_4 = stacks(IDs_above_4)
    #stacking_below_4 = stacks(IDs_below_4)
    stacking_both_4 = stacking_above_4, stacking_below_4
    len_IDs = len_above, len_below

    return stacking_both_4, len_IDs


#stacking_both_105_1075, len_IDs = stack_lims(low_lim3, upp_lim3)
#stacks_above_105_1075, stacks_below_105_1075 = stacking_both_105_1075[0], stacking_both_105_1075[1]
#stacking_both_1075_11, len_IDs_2 = stack_lims(low_lim, upp_lim)
#stacks_above_1075_11, stacks_below_1075_11 = stacking_both_1075_11[0], stacking_both_1075_11[1]
#stacking_both_11_113, len_IDs_3= stack_lims(low_lim2, upp_lim2)
#stacks_above_11_113, stacks_below_11_113 = stacking_both_11_113[0], stacking_both_11_113[1]
stacking_both_105_11, len_IDs_5= stack_lims(low_lim5, upp_lim5)
stacks_above_105_11, stacks_below_105_11 = stacking_both_105_11[0], stacking_both_105_11[1]


def plot_stacks(stack1, stack2, len_above, len_below, lim1, lim2):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(figsize=(20,8))
    #plt.plot(new_wavs, stack2*10**18, color="k", lw=1.5, label = f'below relation (N = {len(IDs_below_1075_11)})')
    plt.plot(new_wavs, stack1, color="r", lw=1.5, ls ='-', label = f'above relation (N = {len_above})')
    plt.plot(new_wavs, stack2, color="k", lw=1.5, label = f'below relation (N = {len_below})')
    plt.plot(new_wavs, np.zeros(len(new_wavs)), color = 'grey', ls = '--', lw = 1.3)
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=17)
    plt.ylabel("Flux", size =17)#$(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^{-1}})}$", size=17)
    plt.xlim(2400, 4200)
    plt.ylim(-0.5 ,3.4)
    plt.legend(fontsize=14, loc = 'upper left')
    plt.title(f'Median stacks above and below normalised van der Wel ETG relation \n ({lim1} $<$ (log(M*) $\leq$ {lim2})', size = 18)# excluding possible AGN (CDFS + UDS)')
    plt.savefig('abovebelow_vdw_'+str(lim1)+'_'+str(lim2)+'_stack_strict_UVJ_d4000less135.pdf')
    plt.close()
    #plt.show()

#plot_stacks(stacks_above_105_1075, stacks_below_105_1075, len_IDs[0], len_IDs[1], low_lim3, upp_lim3)
#plot_stacks(stacks_above_1075_11, stacks_below_1075_11, len_IDs_2[0], len_IDs_2[1], low_lim, upp_lim)
#plot_stacks(stacks_above_11_113, stacks_below_11_113, len_IDs_3[0], len_IDs_3[1], low_lim2, upp_lim2)
plot_stacks(stacks_above_105_11, stacks_below_105_11, len_IDs_5[0], len_IDs_5[1], low_lim5, upp_lim5)










"""
figss, (axss, axss3) = plt.subplots(1, 2, figsize=(14,6))
imss = axss.scatter(strict_masses, np.log10(strict_sizes), s=130, c=D4000_index, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbarss = fig.colorbar(imss, ax=axss)
cbarss.set_label(r'$D_{n}4000$', size=12)
axss.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
axss.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
axss.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)
axss.set_xlim(10.3, 11.3)
axss.set_ylim(-0.4, 1.1)
imss3 = axss3.scatter(strict_masses, np.log10(strict_sizes), s=130, c=H_delta_EW, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )
cbarss3 = fig.colorbar(imss3, ax=axss3)
cbarss3.set_label(r'$EW(H\delta)$', size=12)
axss3.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)# label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
axss3.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
axss3.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
axss3.set_xlim(10.3, 11.3)
axss3.set_ylim(-0.4, 1.1)
plt.savefig("Re_v_M_cbar_both_strict_UVJ_cut.pdf")
plt.close()
"""
"""
figss2, axss2 = plt.subplots(figsize=[12,8.5])
imss2 = axss2.scatter(strict_masses, np.log10(strict_sizes), s=130, c=H_delta_EW, cmap=plt.cm.magma, marker='o', edgecolors='black',linewidth=0.5 )
cbarss2 = figss2.colorbar(imss2, ax=axss2)
#cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
cbarss2.set_label(r'$EW(H\delta)$', size=12)
#ax1.plot(x, y_model, 'k', lw = 0.7, label= 'Shen et al. local ETG relation')
axss2.plot(x, vdw_norm_model, 'k', lw = 1.2, label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
#ax1.plot(x, log_Reff, 'r' ,lw = 1.2, label= 'v. der Wel z = 1.25 ETG relation (not normalised)')
#ax1.scatter(masses, np.log10(Rc), marker='o', s=20, c='r', edgecolors='k')
axss2.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
axss2.set_ylabel(r'$\mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(prop={'size': 10})
plt.title('3D-HST log10(Re) versus log10(M*/Msun)', size = 13)
plt.xlim(10.3, 11.4)
#plt.show()
plt.savefig('Re_v_M*_cbar_Hdelta_stricterUVJ.pdf')
plt.close()
"""
"""
fig5, (ax4, ax5) = plt.subplots(1, 2, figsize=(14,6))
im4 = ax4.scatter((np.log10(strict_sizes)- np.median(np.log10(strict_sizes))), H_delta_EW, s=130, c=D4000_index, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbar4 = fig5.colorbar(im4, ax=ax4)
cbar4.set_label(r'$D_{n}4000$', size=12)
#ax.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax4.set_xlabel(r'$\Delta \mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax4.set_ylabel(r'$EW(H\delta)$', size = 12)
#ax4.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)
#ax4.set_xlim(10.3, 11.3)
#ax4.set_ylim(-0.4, 1.1)
im5 = ax5.scatter((np.log10(strict_sizes)- np.median(np.log10(strict_sizes))), D4000_index, s=130, c=H_delta_EW, cmap=plt.cm.magma, marker='o',linewidth=0.5, alpha = 0.9 )#edgecolors='black'
cbar5 = fig5.colorbar(im5, ax=ax5)
cbar5.set_label(r'$EW(H\delta)$', size=12)
#ax5.plot(x, vdw_norm_model, 'k', lw = 1., ls = '--', alpha = 0.5)# label=f'v. der Wel z = 1.25 ETG relation normalised by f = {round((best_c_vdw-0.22),2)}')
#ax5.set_xlabel(r'$\mathrm{log_{10}{(M*/M_{\odot})}}$', size = 12)
ax5.set_xlabel(r'$\Delta \mathrm{log_{10}{(R_{e}/kpc)}}$', size = 12)
ax5.set_ylabel(r'$D_{n}4000$', size = 12)
#ax5.set_xlim(10.3, 11.3)
#ax5.set_ylim(-0.4, 1.1)
plt.savefig("deltaRe_v_d4000andhdelta_strictUVJ.pdf")
plt.close()

"""

#stack the objects like you did in stacks_M_Re.py and plot them and also do the sigma_50 plots for the stricter UVJ regime

#then for the indices plots - do a plot with EW v D4000 for all strict uvj colours one colour and then everything else in another colour
