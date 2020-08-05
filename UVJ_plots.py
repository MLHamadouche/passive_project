import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import matplotlib
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from stacks import stacks

vandels_cat = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()

UV = vandels_cat["UV_colour_50"]
VJ = vandels_cat["VJ_colour_50"]
SSFR = vandels_cat["ssfr_50"]
ID = vandels_cat["#ID"]
stel_mas = vandels_cat["stellar_mass_50"]
age = vandels_cat["exponential:age_50"]
age1 = vandels_cat["mass_weighted_age_50"]

#color='darkseagreen',
"""
fig, ax = plt.subplots()

im=ax.scatter(VJ, UV, s = 15,c= c,cmap=plt.cm.BuPu_r, alpha=0.8, marker='o', edgecolors='black',  linewidth=0.5)
plt.xlabel(r'V-J')
plt.ylabel(r'U-V')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'ssfr')
#plt.title(r'UVJ diagram for passive objects in VANDELS catalogue 1 < z < 1.5')
plt.savefig('UVJ_ssfrcbar.pdf')
plt.close()
"""

UV_neg_err =  UV - vandels_cat["UV_colour_16"]
UV_pos_err = vandels_cat["UV_colour_84"] -UV
#print(UV_pos_err[0],UV[0], UV_neg_err[0])

VJ_neg_err = VJ - vandels_cat["VJ_colour_16"]
VJ_pos_err = vandels_cat["VJ_colour_84"] - VJ

#print(VJ_pos_err[0],VJ[0],VJ_neg_err[0])

#fig, ax = plt.subplots()
xerr = [VJ_neg_err, VJ_pos_err]
c = age1
csfr = SSFR
fig, ax = plt.subplots()

x = np.linspace(0., 3.9, 237)

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

plt.plot(x, Y, linestyle = '-', lw = 1.3, color = 'k')
plt.plot(x[(Ys<Y3) & (Ys > Y)], Ys[(Ys<Y3) & (Ys > Y)], linestyle = '-', lw = 1.3, color = 'r')
plt.plot(x[(Yss<Y3) & (Yss > Y)], Yss[(Yss<Y3) & (Yss > Y)], linestyle = '-', lw = 1.3, color = 'r')
plt.plot(x[(Ysss<Y3) & (Ysss > Y)], Ysss[(Ysss<Y3) & (Ysss > Y)], linestyle = '-', lw = 1.3, color = 'r')

bin1 = []
bin2 = []
bin3 =[]
bin4 =[]
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
lim2 = (x > x3) & (Y1 > c1 ) & (Y1 < c3)
#lim3 =  (Y3 < c3 ) & (Y3 > c4) & (x > x4)
lim4 =   (Y2 < c4) & (Y2 > c2) & (x < x2)

lim5 =  (Y3 > c4) & (Y3<c3)

if [Y <= 1.3]:
    xmax = (1.3 - 0.49)/0.88
    xmin = 0
    plt.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 1.3)

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.49
    Ymax = 0.88*len(x) +0.49

    plt.vlines(x=1.6, ymin=Ymin, ymax=3.5, ls = '-', lw = 1.3)


#plt.plot(x[lims1], Y[lims1], ls = '-', lw = 0.9, color='r')
plt.plot(x[lim], Y[lim], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim2], Y1[lim2], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim4], Y2[lim4], ls = '-', lw = 1.3, color='r')
plt.plot(x[lim5], Y3[lim5], ls = '-', lw = 1.3, color='r')


im=ax.scatter(VJ, UV, s = 12, c = csfr, cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5)
plt.xlabel(r'V-J')
plt.ylabel(r'U-V')
plt.xlim(0.25,2.5)
plt.xticks(np.arange(0., 2.5, 0.5))
plt.ylim(0.,2.5)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'$\mathrm{log_{10}(sSFR/yr)}$')
plt.title(r'UVJ diagram for VANDELS objects 1 < z < 1.5')
plt.savefig('UVJ_agecbar_box_zspec_test.pdf')
plt.close()
#plt.show()

passive_cut = Table.read('FirstProjectCatalogs/xmatch_spec_derived237objs.fits').to_pandas()
df = pd.DataFrame(passive_cut)
ID_list = df.set_index(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8"))

redshifts = passive_cut['z_spec']
#objects = np.array(passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')+ passive_cut['CAT'].str.decode("utf-8"))

new_ID = []
new_ID2 = []
new_ID3 = []
new_ID4 = []
for i in bin1:
    new_ID.append(i.decode("utf-8"))
for j in bin2:
    new_ID2.append(j.decode("utf-8"))
for k in bin3:
    new_ID3.append(k.decode("utf-8"))
for l in bin2:
    new_ID4.append(l.decode("utf-8"))

objects = new_ID
objects2 = new_ID2
objects3 = new_ID3
objects4 = new_ID4



new_wavs = np.arange(2400, 4200, 1.5) #new wavelength grid
#z_mask = (redshifts[objects] <= 1.5) & (redshifts[objects]>= 1.0)

med_stack = stacks(new_wavs, objects)

def plot_stacks(new_wavs, med_stack, bin_number):
    plt.figure(figsize=(15,7))
    plt.plot(new_wavs, med_stack*10**18, color="black", lw=1.5 )
    plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=15)
    plt.ylabel("Flux $(10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \\AA{^-1})}$", size=15)
    plt.xlim(2300, 4250)
    plt.title('Median Stacked Spectra of galaxies 1 < z < 1.5', size=17)
    plt.savefig('stack_plot_bin'+str(bin_number)+ '.pdf')

plot_stacks(new_wavs, med_stack, 1)
med_stack2 = stacks(new_wavs, objects2)
plot_stacks(new_wavs, med_stack2, 2)
med_stack3 = stacks(new_wavs, objects3)
plot_stacks(new_wavs, med_stack3, 3)

med_stack4 = stacks(new_wavs, objects4)
plot_stacks(new_wavs, med_stack4, 4)



#plt.plot(x, y2, ls = '-.', lw = 0.9, color='k', )
#plt.errorbar(x=VJ, y=UV, yerr= [UV_neg_err, UV_pos_err], xerr=xerr, color='darkseagreen',linestyle=None,
#fmt='o',ms=4, linewidth=1.,solid_capstyle='butt' ,mew=0.5,ecolor='lightcoral', alpha=0.5, mec='black' )
#plt.errorbar(x=VJ, y=UV, yerr= [UV_neg_err, UV_pos_err], xerr=xerr, color='darkseagreen',linestyle=None,
#fmt='o',ms=4, linewidth=1.,solid_capstyle='butt' ,mew=0.5,ecolor='lightcoral', alpha=0.5, mec='black' )
#plt.scatter(VJ, UV, s=12, alpha=0.5, color='darkseagreen', marker='o', edgecolors='black',  linewidth=0.5)
