import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import matplotlib
import extinction

vandels_cat = Table.read("pipes/cats/vandels_cat.fits").to_pandas()

UV = vandels_cat["UV_colour_50"]
VJ = vandels_cat["VJ_colour_50"]
SSFR = np.log(vandels_cat["sfr_50"])
#color='darkseagreen',

fig, ax = plt.subplots()
c = SSFR
im=ax.scatter(VJ, UV, s = 15,c= c,cmap=plt.cm.BuPu_r, alpha=0.8, marker='o', edgecolors='black',  linewidth=0.5)
plt.xlabel(r'V-J')
plt.ylabel(r'U-V')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'ssfr')
plt.title(r'UVJ diagram for passive objects in VANDELS catalogue 1 < z < 1.5')
plt.savefig('UVJ_ssfrcbar.pdf')
plt.close()

UV_neg_err =  UV - vandels_cat["UV_colour_16"]
UV_pos_err = vandels_cat["UV_colour_84"] -UV
print(UV_pos_err[0],UV[0], UV_neg_err[0])

VJ_neg_err = VJ - vandels_cat["VJ_colour_16"]
VJ_pos_err = vandels_cat["VJ_colour_84"] - VJ

print(VJ_pos_err[0],VJ[0],VJ_neg_err[0])

#fig, ax = plt.subplots()
xerr = [VJ_neg_err, VJ_pos_err]
c = SSFR

plt.errorbar(x=VJ, y=UV, yerr= [UV_neg_err, UV_pos_err], xerr=xerr, color='darkseagreen',linestyle=None,
fmt='o',ms=4, linewidth=1.,solid_capstyle='butt' ,mew=0.5,ecolor='lightcoral', alpha=0.5, mec='black' )
#plt.scatter(VJ, UV, s=12, alpha=0.5, color='darkseagreen', marker='o', edgecolors='black',  linewidth=0.5)
#im=ax.scatter(VJ, UV, s = 15,c= c,cmap=plt.cm.BuPu_r, alpha=0.8, marker='o', edgecolors='black',  linewidth=0.5)
plt.xlabel(r'V-J')
plt.ylabel(r'U-V')
#cbar = fig.colorbar(im, ax=ax)
#cbar.set_label(r'ssfr')
plt.title(r'UVJ diagram for VANDELS objects 1 < z < 1.5')
plt.savefig('UVJ_werrs.pdf')
plt.close()
