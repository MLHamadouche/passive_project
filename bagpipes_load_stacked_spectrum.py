import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
from stacking_code import stacks
import LoadData as ld
import bagpipes as pipes

list_of_IDs = ['UDS-HST013753SELECT' ,'UDS-HST019329SELECT' ,'UDS-HST021218SELECT','UDS-HST024977SELECT' ,'CDFS-HST013637SELECT', 'CDFS-HST016360SELECT', 'CDFS-GROUND140661SELECT' ,'CDFS-HST026535SELECT', 'CDFS-HST001851SELECT']

def load_stacked_spectrum(id):

    stacked_spectrum, stack_errors = stacks(list_of_IDs)
    wavs = np.arange(2400, 4200, 1.25)

    flux_errs = stack_errors #need to ask about how to get median error spectrum out of the stacks

    spectrum = np.c_[wavs, stacked_spectrum, flux_errs]

    return spectrum

#list_of_IDs = ['UDS-HST013753SELECT' ,'UDS-HST019329SELECT' ,'UDS-HST021218SELECT','UDS-HST024977SELECT' ,'CDFS-HST013637SELECT', 'CDFS-HST016360SELECT', 'CDFS-GROUND140661SELECT' ,'CDFS-HST026535SELECT', 'CDFS-HST001851SELECT']
"""
def load_stacked_spectrum(list_of_IDs):
    list_of_IDs = ['UDS-HST013753SELECT' ,'UDS-HST019329SELECT' ,'UDS-HST021218SELECT','UDS-HST024977SELECT' ,'CDFS-HST013637SELECT', 'CDFS-HST016360SELECT', 'CDFS-GROUND140661SELECT' ,'CDFS-HST026535SELECT', 'CDFS-HST001851SELECT']
    stacked_spectrum = stacks(list_of_IDs)
    wavs = np.arange(2400, 4200, 1.25)

    flux_errs = np.ones(len(stacked_spectrum))*0.1 #need to ask about how to get median error spectrum out of the stacks

    spectrum = np.c_[wavs, stacked_spectrum, flux_errs]

    return spectrum
"""
dblplaw = {}
dblplaw["tau"] = (0., 15.)
dblplaw["alpha"] = (0.01, 1000.)
dblplaw["beta"] = (0.01, 1000.)
dblplaw["alpha_prior"] = "log_10"
dblplaw["beta_prior"] = "log_10"
dblplaw["massformed"] = (1., 15.)
dblplaw["metallicity"] = (0.1, 2.)
dblplaw["metallicity_prior"] = "log_10"

nebular = {}
nebular["logU"] = -3.

dust = {}
dust["type"] = "CF00"
dust["eta"] = 2.
dust["Av"] = (0., 2.0)
dust["n"] = (0.3, 2.5)
dust["n_prior"] = "Gaussian"
dust["n_prior_mu"] = 0.7
dust["n_prior_sigma"] = 0.3

fit_instructions = {}
fit_instructions["redshift"] = (0.75, 1.25)
fit_instructions["t_bc"] = 0.01
fit_instructions["redshift_prior"] = "Gaussian"
fit_instructions["redshift_prior_mu"] = 0.9
fit_instructions["redshift_prior_sigma"] = 0.05
fit_instructions["dblplaw"] = dblplaw
fit_instructions["nebular"] = nebular
fit_instructions["dust"] = dust
fit_instructions["veldisp"] = (1., 1000.)   #km/s
fit_instructions["veldisp_prior"] = "log_10"
calib = {}
calib["type"] = "polynomial_bayesian"

calib["0"] = (0.5, 1.5) # Zero order is centred on 1, at which point there is no change to the spectrum.
calib["0_prior"] = "Gaussian"
calib["0_prior_mu"] = 1.0
calib["0_prior_sigma"] = 0.25

calib["1"] = (-0.5, 0.5) # Subsequent orders are centred on zero.
calib["1_prior"] = "Gaussian"
calib["1_prior_mu"] = 0.
calib["1_prior_sigma"] = 0.25

calib["2"] = (-0.5, 0.5)
calib["2_prior"] = "Gaussian"
calib["2_prior_mu"] = 0.
calib["2_prior_sigma"] = 0.25

fit_instructions["calib"] = calib
mlpoly = {}
mlpoly["type"] = "polynomial_max_like"
mlpoly["order"] = 2
noise = {}
noise["type"] = "white_scaled"
noise["scaling"] = (1., 10.)
noise["scaling_prior"] = "log_10"
fit_instructions["noise"] = noise


id = '0001'

galaxy = pipes.galaxy(id, load_data =load_stacked_spectrum, photometry_exists=False, spectrum_exists = True)

fit = pipes.fit(galaxy, fit_instructions, run="spectroscopy")

fit.fit(verbose=False)

fig = fit.plot_spectrum_posterior(save=False, show=True)
fig = fit.plot_calibration(save=False, show=True)
fig = fit.plot_sfh_posterior(save=False, show=True)
fig = fit.plot_corner(save=False, show=True)
