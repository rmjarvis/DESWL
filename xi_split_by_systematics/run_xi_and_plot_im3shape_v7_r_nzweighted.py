import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pylab
import os
import sys
import pyfits
import glob
import treecorr


## Full im3shape_v7_r
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_03_z_13.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_03_z_13.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r.out'
treecorr.corr2(config)

####### AIRMASS ##############
## im3shape_v7_r Airmass Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_airmass_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_airmass_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_airmass_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_airmass_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r Airmass Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_airmass_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_airmass_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_airmass_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_airmass_r_lower_nzweighted.out'
treecorr.corr2(config)

####### EXPOSURE TIME ##############
## im3shape_v7_r exptime Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_exptime_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_exptime_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_exptime_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_exptime_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r exptime Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_exptime_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_exptime_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_exptime_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_exptime_r_lower_nzweighted.out'
treecorr.corr2(config)

####### FWHM ##############
## im3shape_v7_r fwhm Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_fwhm_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_fwhm_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_fwhm_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_fwhm_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r fwhm Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_fwhm_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_fwhm_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_fwhm_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_fwhm_r_lower_nzweighted.out'
treecorr.corr2(config)

####### Maglimit ##############
## im3shape_v7_r maglimit Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_maglimit_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_maglimit_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_maglimit_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_maglimit_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r maglimit Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_maglimit_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_maglimit_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_maglimit_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_maglimit_r_lower_nzweighted.out'
treecorr.corr2(config)

####### Sky Brightness ##############
## im3shape_v7_r skybrite Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skybrite_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skybrite_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skybrite_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skybrite_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r skybrite Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skybrite_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skybrite_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skybrite_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skybrite_r_lower_nzweighted.out'
treecorr.corr2(config)

####### Sky Sigma ##############
## im3shape_v7_r skysigma Upper
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skysigma_r_upper_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skysigma_r_upper_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skysigma_r_upper_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skysigma_r_upper_nzweighted.out'
treecorr.corr2(config)

## im3shape_v7_r skysigma Lower
config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_shearhsear_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skysigma_r_lower_nzweighted.dat'
config['gg_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skysigma_r_lower_nzweighted.out'
treecorr.corr2(config)

config = treecorr.read_config('/Users/drgk/DES/SV_tests/split_xi/sample.params_kappakappa_im3shape')
config['file_name'] = '/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skysigma_r_lower_nzweighted.dat'
config['kk_file_name'] = '/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skysigma_r_lower_nzweighted.out'
treecorr.corr2(config)


x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r.out')
r_gg = x[:,1]
xi_plus_gg = x[:,2]
xi_minus_gg = x[:,3]
sigma_gg = x[:,6]
weight_gg = x[:,7]
npairs_gg = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r.out')
r_kk = x[:,1]
xi_kk = x[:,2]
sigma_kk = x[:,3]
weight_kk = x[:,4]
npairs_kk = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_airmass_r_upper_nzweighted.out')
r_gg_airmass_upper = x[:,1]; xi_plus_gg_airmass_upper = x[:,2]
xi_minus_gg_airmass_upper = x[:,3]; sigma_gg_airmass_upper = x[:,6]
weight_gg_airmass_upper = x[:,7]; npairs_gg_airmass_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_airmass_r_upper_nzweighted.out')
r_kk_airmass_upper = x[:,1]; xi_kk_airmass_upper = x[:,2]
sigma_kk_airmass_upper = x[:,3]; weight_kk_airmass_upper = x[:,4]; npairs_kk_airmass_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_airmass_r_lower_nzweighted.out')
r_gg_airmass_lower = x[:,1]; xi_plus_gg_airmass_lower = x[:,2]
xi_minus_gg_airmass_lower = x[:,3]; sigma_gg_airmass_lower = x[:,6]
weight_gg_airmass_lower = x[:,7]; npairs_gg_airmass_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_airmass_r_lower_nzweighted.out')
r_kk_airmass_lower = x[:,1]; xi_kk_airmass_lower = x[:,2]
sigma_kk_airmass_lower = x[:,3]; weight_kk_airmass_lower = x[:,4]; npairs_kk_airmass_lower = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_exptime_r_upper_nzweighted.out')
r_gg_exptime_upper = x[:,1]; xi_plus_gg_exptime_upper = x[:,2]
xi_minus_gg_exptime_upper = x[:,3]; sigma_gg_exptime_upper = x[:,6]
weight_gg_exptime_upper = x[:,7]; npairs_gg_exptime_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_exptime_r_upper_nzweighted.out')
r_kk_exptime_upper = x[:,1]; xi_kk_exptime_upper = x[:,2]
sigma_kk_exptime_upper = x[:,3]; weight_kk_exptime_upper = x[:,4]; npairs_kk_exptime_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_exptime_r_lower_nzweighted.out')
r_gg_exptime_lower = x[:,1]; xi_plus_gg_exptime_lower = x[:,2]
xi_minus_gg_exptime_lower = x[:,3]; sigma_gg_exptime_lower = x[:,6]
weight_gg_exptime_lower = x[:,7]; npairs_gg_exptime_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_exptime_r_lower_nzweighted.out')
r_kk_exptime_lower = x[:,1]; xi_kk_exptime_lower = x[:,2]
sigma_kk_exptime_lower = x[:,3]; weight_kk_exptime_lower = x[:,4]; npairs_kk_exptime_lower = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_fwhm_r_upper_nzweighted.out')
r_gg_fwhm_upper = x[:,1]; xi_plus_gg_fwhm_upper = x[:,2]
xi_minus_gg_fwhm_upper = x[:,3]; sigma_gg_fwhm_upper = x[:,6]
weight_gg_fwhm_upper = x[:,7]; npairs_gg_fwhm_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_fwhm_r_upper_nzweighted.out')
r_kk_fwhm_upper = x[:,1]; xi_kk_fwhm_upper = x[:,2]
sigma_kk_fwhm_upper = x[:,3]; weight_kk_fwhm_upper = x[:,4]; npairs_kk_fwhm_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_fwhm_r_lower_nzweighted.out')
r_gg_fwhm_lower = x[:,1]; xi_plus_gg_fwhm_lower = x[:,2]
xi_minus_gg_fwhm_lower = x[:,3]; sigma_gg_fwhm_lower = x[:,6]
weight_gg_fwhm_lower = x[:,7]; npairs_gg_fwhm_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_fwhm_r_lower_nzweighted.out')
r_kk_fwhm_lower = x[:,1]; xi_kk_fwhm_lower = x[:,2]
sigma_kk_fwhm_lower = x[:,3]; weight_kk_fwhm_lower = x[:,4]; npairs_kk_fwhm_lower = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_maglimit_r_upper_nzweighted.out')
r_gg_maglimit_upper = x[:,1]; xi_plus_gg_maglimit_upper = x[:,2]
xi_minus_gg_maglimit_upper = x[:,3]; sigma_gg_maglimit_upper = x[:,6]
weight_gg_maglimit_upper = x[:,7]; npairs_gg_maglimit_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_maglimit_r_upper_nzweighted.out')
r_kk_maglimit_upper = x[:,1]; xi_kk_maglimit_upper = x[:,2]
sigma_kk_maglimit_upper = x[:,3]; weight_kk_maglimit_upper = x[:,4]; npairs_kk_maglimit_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_maglimit_r_lower_nzweighted.out')
r_gg_maglimit_lower = x[:,1]; xi_plus_gg_maglimit_lower = x[:,2]
xi_minus_gg_maglimit_lower = x[:,3]; sigma_gg_maglimit_lower = x[:,6]
weight_gg_maglimit_lower = x[:,7]; npairs_gg_maglimit_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_maglimit_r_lower_nzweighted.out')
r_kk_maglimit_lower = x[:,1]; xi_kk_maglimit_lower = x[:,2]
sigma_kk_maglimit_lower = x[:,3]; weight_kk_maglimit_lower = x[:,4]; npairs_kk_maglimit_lower = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skybrite_r_upper_nzweighted.out')
r_gg_skybrite_upper = x[:,1]; xi_plus_gg_skybrite_upper = x[:,2]
xi_minus_gg_skybrite_upper = x[:,3]; sigma_gg_skybrite_upper = x[:,6]
weight_gg_skybrite_upper = x[:,7]; npairs_gg_skybrite_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skybrite_r_upper_nzweighted.out')
r_kk_skybrite_upper = x[:,1]; xi_kk_skybrite_upper = x[:,2]
sigma_kk_skybrite_upper = x[:,3]; weight_kk_skybrite_upper = x[:,4]; npairs_kk_skybrite_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skybrite_r_lower_nzweighted.out')
r_gg_skybrite_lower = x[:,1]; xi_plus_gg_skybrite_lower = x[:,2]
xi_minus_gg_skybrite_lower = x[:,3]; sigma_gg_skybrite_lower = x[:,6]
weight_gg_skybrite_lower = x[:,7]; npairs_gg_skybrite_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skybrite_r_lower_nzweighted.out')
r_kk_skybrite_lower = x[:,1]; xi_kk_skybrite_lower = x[:,2]
sigma_kk_skybrite_lower = x[:,3]; weight_kk_skybrite_lower = x[:,4]; npairs_kk_skybrite_lower = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skysigma_r_upper_nzweighted.out')
r_gg_skysigma_upper = x[:,1]; xi_plus_gg_skysigma_upper = x[:,2]
xi_minus_gg_skysigma_upper = x[:,3]; sigma_gg_skysigma_upper = x[:,6]
weight_gg_skysigma_upper = x[:,7]; npairs_gg_skysigma_upper = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skysigma_r_upper_nzweighted.out')
r_kk_skysigma_upper = x[:,1]; xi_kk_skysigma_upper = x[:,2]
sigma_kk_skysigma_upper = x[:,3]; weight_kk_skysigma_upper = x[:,4]; npairs_kk_skysigma_upper = x[:,5]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/gg_galshear_im3shape_v7_r_skysigma_r_lower_nzweighted.out')
r_gg_skysigma_lower = x[:,1]; xi_plus_gg_skysigma_lower = x[:,2]
xi_minus_gg_skysigma_lower = x[:,3]; sigma_gg_skysigma_lower = x[:,6]
weight_gg_skysigma_lower = x[:,7]; npairs_gg_skysigma_lower = x[:,8]

x = np.loadtxt('/Users/drgk/DES/SV_tests/split_xi/kk_galshear_im3shape_v7_r_skysigma_r_lower_nzweighted.out')
r_kk_skysigma_lower = x[:,1]; xi_kk_skysigma_lower = x[:,2]
sigma_kk_skysigma_lower = x[:,3]; weight_kk_skysigma_lower = x[:,4]; npairs_kk_skysigma_lower = x[:,5]

##### CORRECT FOR M 
xi_plus_gg_corrected = xi_plus_gg/xi_kk
xi_minus_gg_corrected = xi_minus_gg/xi_kk
xi_plus_gg_airmass_upper_corrected = xi_plus_gg_airmass_upper/xi_kk_airmass_upper
xi_plus_gg_airmass_lower_corrected = xi_plus_gg_airmass_lower/xi_kk_airmass_lower
xi_plus_gg_fwhm_upper_corrected = xi_plus_gg_fwhm_upper/xi_kk_fwhm_upper
xi_plus_gg_fwhm_lower_corrected = xi_plus_gg_fwhm_lower/xi_kk_fwhm_lower
xi_plus_gg_maglimit_upper_corrected = xi_plus_gg_maglimit_upper/xi_kk_maglimit_upper
xi_plus_gg_maglimit_lower_corrected = xi_plus_gg_maglimit_lower/xi_kk_maglimit_lower
xi_plus_gg_exptime_upper_corrected = xi_plus_gg_exptime_upper/xi_kk_exptime_upper
xi_plus_gg_exptime_lower_corrected = xi_plus_gg_exptime_lower/xi_kk_exptime_lower
xi_plus_gg_skybrite_upper_corrected = xi_plus_gg_skybrite_upper/xi_kk_skybrite_upper
xi_plus_gg_skybrite_lower_corrected = xi_plus_gg_skybrite_lower/xi_kk_skybrite_lower
xi_plus_gg_skysigma_upper_corrected = xi_plus_gg_skysigma_upper/xi_kk_skysigma_upper
xi_plus_gg_skysigma_lower_corrected = xi_plus_gg_skysigma_lower/xi_kk_skysigma_lower

################# PLOT FULL CATALOGUE ###################

plt.figure()
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg,xi_minus_gg_corrected,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$ Full')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
#plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
plt.title('Im3shape v7 r, infoflag==0, errorflag==0, $0.3 < z < 1.3$, NBC Corrected',fontsize=18)
plt.show()

################# PLOT SYSTEMATICS ###################
# Airmass
plt.figure()
plt.subplot(3,2,1)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_airmass_upper_corrected,sigma_gg_airmass_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Airmass Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_airmass_lower_corrected,sigma_gg_airmass_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Airmass Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
plt.title('Im3shape v7 r, NBC Corrected, $n(z)$ re-weighted',fontsize=18)
plt.show()


# Exposure Time
#plt.figure()
plt.subplot(3,2,2)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_exptime_upper_corrected,sigma_gg_exptime_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_exptime_lower_corrected,sigma_gg_exptime_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Exposure Time',fontsize=18)
plt.show()


# FWHM
#plt.figure()
plt.subplot(3,2,3)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_fwhm_upper_corrected,sigma_gg_fwhm_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ FWHM Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_fwhm_lower_corrected,sigma_gg_fwhm_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ FWHM Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, FWHM',fontsize=18)
plt.show()


# Magnitude Limit
#plt.figure()
plt.subplot(3,2,4)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_maglimit_upper_corrected,sigma_gg_maglimit_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_maglimit_lower_corrected,sigma_gg_maglimit_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Magnitude Limit',fontsize=18)
plt.show()


# Sky Brightness
#plt.figure()
plt.subplot(3,2,5)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_skybrite_upper_corrected,sigma_gg_skybrite_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skybrite_lower_corrected,sigma_gg_skybrite_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Brightness',fontsize=18)
plt.show()


# Sky Sigma
#plt.figure()
plt.subplot(3,2,6)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_skysigma_upper_corrected,sigma_gg_skysigma_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skysigma_lower_corrected,sigma_gg_skysigma_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Sigma',fontsize=18)
plt.show()


################# PLOT SYSTEMATICS - LOG SCALE ###################
# Airmass
plt.figure()
plt.subplot(3,2,1)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_airmass_upper_corrected,sigma_gg_airmass_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Airmass Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_airmass_lower_corrected,sigma_gg_airmass_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Airmass Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
plt.title('Im3shape v7 r, NBC Corrected, $n(z)$ re-weighted',fontsize=18)
plt.show()


# Exposure Time
#plt.figure()
plt.subplot(3,2,2)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_exptime_upper_corrected,sigma_gg_exptime_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_exptime_lower_corrected,sigma_gg_exptime_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Exposure Time',fontsize=18)
plt.show()


# FWHM
#plt.figure()
plt.subplot(3,2,3)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_fwhm_upper_corrected,sigma_gg_fwhm_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ FWHM Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_fwhm_lower_corrected,sigma_gg_fwhm_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ FWHM Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.yscale('log',nonposy='clip')
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, FWHM',fontsize=18)
plt.show()


# Magnitude Limit
#plt.figure()
plt.subplot(3,2,4)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_maglimit_upper_corrected,sigma_gg_maglimit_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_maglimit_lower_corrected,sigma_gg_maglimit_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Magnitude Limit',fontsize=18)
plt.show()


# Sky Brightness
#plt.figure()
plt.subplot(3,2,5)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_skybrite_upper_corrected,sigma_gg_skybrite_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skybrite_lower_corrected,sigma_gg_skybrite_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Brightness',fontsize=18)
plt.show()


# Sky Sigma
#plt.figure()
plt.subplot(3,2,6)
plt.errorbar(r_gg,xi_plus_gg_corrected,sigma_gg,color='blue',marker='o',linestyle='',label='$\\xi_+$ Full')
plt.errorbar(r_gg*1.02,xi_plus_gg_skysigma_upper_corrected,sigma_gg_skysigma_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skysigma_lower_corrected,sigma_gg_skysigma_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.ylim([0,0.00009])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Sigma',fontsize=18)
plt.show()



sys.exit()

#################### PLOT RESIDUALS #####################

# Airmass
plt.figure()
plt.subplot(3,2,1)
plt.errorbar(r_gg*1.02,xi_plus_gg_airmass_upper-xi_plus_gg,sigma_gg_airmass_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Airmass Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_airmass_lower-xi_plus_gg,sigma_gg_airmass_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Airmass Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
plt.title('Im3shape v7 r',fontsize=18)
plt.show()

# Exposure Time
#plt.figure()
plt.subplot(3,2,2)
plt.errorbar(r_gg*1.02,xi_plus_gg_exptime_upper-xi_plus_gg,sigma_gg_exptime_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_exptime_lower-xi_plus_gg,sigma_gg_exptime_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Exposure Time',fontsize=18)
plt.show()


# FWHM
#plt.figure()
plt.subplot(3,2,3)
plt.errorbar(r_gg*1.02,xi_plus_gg_fwhm_upper-xi_plus_gg,sigma_gg_fwhm_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ FWHM Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_fwhm_lower-xi_plus_gg,sigma_gg_fwhm_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ FWHM Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, FWHM',fontsize=18)
plt.show()


# Magnitude Limit
#plt.figure()
plt.subplot(3,2,4)
plt.errorbar(r_gg*1.02,xi_plus_gg_maglimit_upper-xi_plus_gg,sigma_gg_maglimit_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_maglimit_lower-xi_plus_gg,sigma_gg_maglimit_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Magnitude Limit',fontsize=18)
plt.show()


# Sky Brightness
#plt.figure()
plt.subplot(3,2,5)
plt.errorbar(r_gg*1.02,xi_plus_gg_skybrite_upper-xi_plus_gg,sigma_gg_skybrite_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skybrite_lower-xi_plus_gg,sigma_gg_skybrite_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Brightness',fontsize=18)
plt.show()


# Sky Sigma
#plt.figure()
plt.subplot(3,2,6)
plt.errorbar(r_gg*1.02,xi_plus_gg_skysigma_upper-xi_plus_gg,sigma_gg_skysigma_upper,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Upper')
plt.errorbar(r_gg*1.04,xi_plus_gg_skysigma_lower-xi_plus_gg,sigma_gg_skysigma_lower,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$\\xi_{+} - \\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Sigma',fontsize=18)
plt.show()



#################### PLOT RESIDUALS - LOG SCALE #####################

# Airmass
plt.figure()
plt.subplot(3,2,1)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_airmass_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_airmass_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ Airmass Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_airmass_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_airmass_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ Airmass Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
plt.title('Im3shape v7 r',fontsize=18)
plt.show()

# Exposure Time
#plt.figure()
plt.subplot(3,2,2)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_exptime_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_exptime_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_exptime_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_exptime_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ Exposure Time Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Exposure Time',fontsize=18)
plt.show()


# FWHM
#plt.figure()
plt.subplot(3,2,3)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_fwhm_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_fwhm_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ FWHM Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_fwhm_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_fwhm_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ FWHM Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, FWHM',fontsize=18)
plt.show()


# Magnitude Limit
#plt.figure()
plt.subplot(3,2,4)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_maglimit_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_maglimit_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_maglimit_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_maglimit_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ Magnitude Limit Lower')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Magnitude Limit',fontsize=18)
plt.show()


# Sky Brightness
#plt.figure()
plt.subplot(3,2,5)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_skybrite_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_skybrite_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_skybrite_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_skybrite_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Brightness Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Brightness',fontsize=18)
plt.show()


# Sky Sigma
#plt.figure()
plt.subplot(3,2,6)
plt.errorbar(r_gg*1.02,abs(xi_plus_gg_skysigma_upper-xi_plus_gg)/xi_plus_gg,sigma_gg_skysigma_upper/xi_plus_gg,color='green',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Upper')
plt.errorbar(r_gg*1.04,abs(xi_plus_gg_skysigma_lower-xi_plus_gg)/xi_plus_gg,sigma_gg_skysigma_lower/xi_plus_gg,color='orange',marker='o',linestyle='',label='$\\xi_+$ Sky Sigma Lower')
#plt.errorbar(r_gg,xi_minus_gg,sigma_gg,color='red',marker='o',linestyle='',label='$\\xi_-$')
plt.plot(np.linspace(0.0001,300),np.zeros(np.shape(np.linspace(0.0001,300))),color='black',linestyle='--')
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.xlim([1,210])
plt.legend(loc=0,frameon=False)
plt.xlabel('$\\theta$ [arcmin]',fontsize=18)
plt.ylabel('$|\\xi_{+} - \\xi_{+,full}|/\\xi_{+,full}$',fontsize=18)
#plt.title('Im3shape v7 r, Sky Sigma',fontsize=18)
plt.show()





