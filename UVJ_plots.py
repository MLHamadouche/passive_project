import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import matplotlib
import extinction

vandels_cat = Table.read("pipes/cats/vandels_cat_zspec.fits").to_pandas()

UV = vandels_cat["UV_colour_50"]
VJ = vandels_cat["VJ_colour_50"]
SSFR = np.log(vandels_cat["sfr_50"])
stel_mas = vandels_cat["stellar_mass_50"]
age = vandels_cat["exponential:age_50"]
age1 = vandels_cat["mass_weighted_age_50"]
#color='darkseagreen',
"""
fig, ax = plt.subplots()
c = SSFR
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

fig, ax = plt.subplots()

x = np.linspace(0., 3.5, 1000)
Y = 0.88*x + 0.49
y2 = 0.88*x + 0.69
lims = (Y > 1.3) & (x<1.6)
plt.plot(x[lims], Y[lims], linestyle = '-', lw = 0.8, color = 'k')

if [Y <= 1.3]:
    xmax = (1.3 - 0.49)/0.88
    xmin = 0
    plt.hlines(y = 1.3, xmin = xmin, xmax = xmax, ls = '-', lw = 0.8)

if [x>=1.6]:
    Ymin = 0.88*1.6 +0.49
    Ymax = 0.88*len(x) +0.49

    plt.vlines(x=1.6, ymin=Ymin, ymax=3.5, ls = '-', lw = 0.8)
plt.plot(x, Y, ls = '-', lw = 0.9, color='k')
#plt.plot(x, y2, ls = '-.', lw = 0.9, color='k', )
#plt.errorbar(x=VJ, y=UV, yerr= [UV_neg_err, UV_pos_err], xerr=xerr, color='darkseagreen',linestyle=None,
#fmt='o',ms=4, linewidth=1.,solid_capstyle='butt' ,mew=0.5,ecolor='lightcoral', alpha=0.5, mec='black' )
#plt.errorbar(x=VJ, y=UV, yerr= [UV_neg_err, UV_pos_err], xerr=xerr, color='darkseagreen',linestyle=None,
#fmt='o',ms=4, linewidth=1.,solid_capstyle='butt' ,mew=0.5,ecolor='lightcoral', alpha=0.5, mec='black' )
#plt.scatter(VJ, UV, s=12, alpha=0.5, color='darkseagreen', marker='o', edgecolors='black',  linewidth=0.5)
im=ax.scatter(VJ, UV, s = 25,c= c,cmap=plt.cm.magma, marker='o', edgecolors='black',  linewidth=0.5)
plt.xlabel(r'V-J')
plt.ylabel(r'U-V')
plt.xlim(0,2.5)
plt.ylim(0.,3.0)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'log(M/M$_{\odot}$)')
cbar.set_label(r'Age (Gyrs)')
plt.title(r'UVJ diagram for VANDELS objects 1 < z < 1.5')
plt.savefig('UVJ_agecbar_box_zspec.pdf')
plt.close()
#plt.show()
