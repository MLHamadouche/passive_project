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

passive_cut = Table.read('FirstProjectCatalogs/x_match_final_passive_sample_edit.fits').to_pandas()
df2 = pd.DataFrame(passive_cut)#, index = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8")) )
ID_list2 = df2.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))

all_obs = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + passive_cut['CAT'].str.decode("utf-8"))

missing = list(set(all_obs).difference(agn))

new_df=pd.concat([df1,df2]).drop_duplicates(subset = 'ID_1', keep=False)
#print(new_df)

ID_list = new_df.set_index(new_df['FIELD'].str.decode("utf-8").str.rstrip() + new_df['ID_1'].astype(str).str.pad(6, side='left', fillchar='0') + new_df['CAT'].str.decode("utf-8"))


def plot_stacks(new_wavs, med_stack, bin_number):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, med_stack*10**18, color="black", lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    #plt.xlim(2300, 4250)
    plt.ylim(0 ,2.0)
    plt.title('Median Stacked Spectra of galaxies 1 < z < 2.5', size=17)
    plt.savefig('stack_plot_bin'+str(bin_number)+ '.pdf')
    plt.close()


fig, ax1 = plt.subplots(figsize=[20,15])

# equivalent but more general
ax2 = fig.add_subplot(411, position=[0.36, 0.48, 0.38, 0.11])
ax3 = fig.add_subplot(412, position=[0.36, 0.37, 0.38, 0.11],sharex=ax2, sharey=ax2)
ax4 = fig.add_subplot(413, position=[0.36, 0.26, 0.38, 0.11],sharex=ax2, sharey=ax2)
ax5 = fig.add_subplot(414, position=[0.36, 0.15, 0.38, 0.11],sharex=ax2, sharey=ax2)

x = np.linspace(0., 3.9, 1000)

Y = 0.88*x + 0.49
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

#m1, n1, m2, n2 = 1.21/1.76, 1.295, 1.25/1.76, 1.725
#m3,n3, m4, n4= 1.66/1.76, 1.52, 1.7/1.76, 1.95
#m5, n5, m6, n6 = 2.11/1.76, 1.745, 2.15/1.76, 2.175
#m7, n7, m8, n8 = 2.56/1.76, 1.97, 2.6/1.76, 2.4

#bottom left top right corners of each rectangle
m1, n1, m2, n2 = x2, c2, 1.25/1.76, 1.725
m3,n3, m4, n4= 1.86/1.76, 1.42, 1.7/1.76, 1.95
m5, n5, m6, n6 = 2.31/1.76, 1.645, 2.15/1.76, 2.175
m7, n7, m8, n8 = 2.76/1.76, 1.87, 2.6/1.76, 2.4

for i in range(len(UV)):
    point = Point(VJ[i], UV[i])
    #print(point1) polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
    polygon = Polygon([(m1, n1), (x4, c4), (m2, n2), (m3, n3)])
    polygon2 = Polygon([(m3, n3), (m2, n2), (m4, n4), (m5, n5)])
    polygon3 = Polygon([(m5, n5), (m4, n4), (m6, n6), (m7, n7)])
    polygon4 = Polygon([(m7, n7), (m6, n6), (m8, n8), (x1, c1)])
    if polygon.contains(point):
        bin1.append(objID[i])
    elif polygon2.contains(point):
        bin2.append(objID[i])
    elif polygon3.contains(point):
        bin3.append(objID[i])
    elif polygon4.contains(point):
        bin4.append(objID[i])
    #print(ID[i])
print(len(bin1), len( bin2), len(bin3),len(bin4))
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

#print( bin2, bin3, bin4)
lim = (Y < c1) & (x > x2)
lim2 = (x > x3) & (Y1 > c1 ) & (Y1 < c3)
#lim3 =  (Y3 < c3 ) & (Y3 > c4) & (x > x4)
lim4 =   (Y2 < c4) & (Y2 > c2) & (x < x2)

lim5 =  (Y3 > c4) & (Y3<c3)

if [Y <= 1.3]:
    xmax = (1.3 - 0.49)/0.88
    xmin = 0
    ax1.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 2.)#1.3

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.49
    Ymax = 0.88*len(x) +0.49

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
ax1.set_title('UVJ diagram for VANDELS objects 1 < z < 2.5', size = 25)
ax1.set_xlabel(r'V-J', size=25)
ax1.set_ylabel(r'U-V',size=25)
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

age_bin1 = []
for ID in new_ID1:
    age_bin1.append(IDs.loc[ID, "exponential:age_50"])
age_bin2 = []
for ID2 in new_ID2:
    age_bin2.append(IDs.loc[ID2, "exponential:age_50"])
age_bin3 = []
for ID3 in new_ID3:
    age_bin3.append(IDs.loc[ID3, "exponential:age_50"])
age_bin4 = []
for ID4 in new_ID4:
    age_bin4.append(IDs.loc[ID4, "exponential:age_50"])

print(np.median(age_bin1), np.median(age_bin2), np.median(age_bin3), np.median(age_bin4))

print(len(new_ID1), len(new_ID2), len(new_ID3), len(new_ID4))

input()

med_stack1 = stacks(new_ID1, ID_list)
#plot_stacks(new_waves, median_stack, 9)

med_stack2 = stacks(new_ID2, ID_list)
#plot_stacks(new_waves, med_stack2, 10)

med_stack3 = stacks(new_ID3, ID_list)
#plot_stacks(new_waves, med_stack3, 11)

med_stack4 = stacks(new_ID4, ID_list)

#print(new_ID4)
#plot_stacks(new_waves, med_stack4, 12)


ax2.plot(wavs_stack, med_stack1*10**18, 'g',lw=0.8 )
ax3.plot(wavs_stack, med_stack2*10**18, 'g',lw=0.8 )
ax4.plot(wavs_stack, med_stack3*10**18, 'g',lw=0.8)
ax5.plot(wavs_stack, med_stack4*10**18, 'g',lw=0.8)

#ax2.tick_params(axis='x',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#      bottom='off',      # ticks along the bottom edge are off
#        top='off',         # ticks along the top edge are off
#    labelbottom='off'  # labels along the bottom edge are off)
#   )
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
#ax2.set_ylabel('Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$', size=10)
ax2.set_title('Median Stacked Spectra of galaxies 1 < z < 2.5', size=15)

#ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#ax2 = fig.add_axes([0.72, 0.72, 0.16, 0.16])

plt.savefig('UVJ_with_stacked_spectra.pdf')
plt.close()
