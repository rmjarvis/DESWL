import os
import yaml, argparse, sys, logging , pyfits
import numpy as np
import shutil, warnings; warnings.simplefilter("once")

logger = logging.getLogger("homogenise_nz")
logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s ","%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.propagate = False


def get_weights(split_list_cat,redshift_col='ZB',normed_hist=True,label='systematics',plots=True):
	"""
	@param split_list_cat list columns for redshift of objects, each one corresponding to a bin split by systematics.
	Example split_list_cat = [cat_airmass_low_z,cat_airmass_high_z].
	For example cat_airmass = np.array(pyfits.getdata('DESXXXX-XXXX_cat.fits'))
	@param redshift_col each of these catalogs has to have a column corresponding to redshift point estimate.
	@param normed_hist if to use reweighting with normalised histograms (not much difference)
    @label title to pass for plots (example: 'airmass')
	@return return a list with weight vector for each systematic bin
	"""
    if plots: import pylab

    bins_photoz = np.linspace(config['photoz_min'],config['photoz_max'],config['photoz_nbins'])

    n_bins = len(split_list_cat)
    list_hist_z = []

    normed_hist=True

    logger.info('getting n(z) for %d systematic bins' % n_bins)
    for isb,vbin in enumerate(split_list_cat):

        nbin = len(vbin)       

        hist_z, bins_z = np.histogram(vbin,bins=bins_photoz,normed=normed_hist)
        list_hist_z.append(hist_z)
        if plots:  = pl.plot(get_bins_centers(bins_z),hist_z,label=r'bin %d '%(isb))

        logger.info('bin %d n_gal %d' % (isb,nbin))

    mean_z = np.mean(np.array(list_hist_z),axis=0)
    if plots: pl.plot(get_bins_centers(bins_z),mean_z,label='mean',lw=5)
    if plots: pl.legend()
    if plots: pl.xlabel('z')
    if plots: pl.ylabel('number of galaxies')
    if plots: pl.title('n(z) for %s bins' % label)
    if plots: filename_fig = 'nz_for_bins.%s.eps' % label
    if plots: pl.savefig(filename_fig)
    if plots: logger.info('wrote %s' % (filename_fig))

    if plots: pl.figure()
    list_weights = []
    logger.info('getting weights for %d systematic bins' % n_bins)
    for isb,vbin in enumerate(split_list_cat):

        nbin = len(res_bin)       

        hist_z, bins_z = np.histogram(vbin,bins=bins_photoz,normed=normed_hist)

        weights = mean_z/hist_z

        list_weights.append(weights)

        if plots: pl.plot(get_bins_centers(bins_z),weights,label=r'bin %d'%(isb))
        
        logger.info('bin %d n_gal %d' % (isb,nbin))

    if plots: pl.xlabel('z')
    if plots: pl.ylabel('weight value')
    if plots: pl.title('w(z) values for %s bins' % label)
    if plots: pl.legend()
    if plots: filename_fig = 'wz_for_bins.%s.eps' % label
    if plots: logger.info('wrote %s' % (filename_fig))


    if plots: pl.figure()
    if plots: pl.plot(get_bins_centers(bins_z),mean_z,label='mean',lw=5)
    digitize_bins = bins_photoz
    list_weights_use = []
    logger.info('getting reweighted histogrgams for %d systematic bins' % n_bins)
    for isb,vbin in enumerate(split_list_cat):

        nbin = len(res_bin)       

        inds = np.digitize(vbin, digitize_bins)

        weights_in_bin = list_weights[isb].copy()
        weights_in_bin = np.append(weights_in_bin,0)
        weights_in_bin = np.insert(weights_in_bin,0,0)
        weights_use = weights_in_bin[inds]
        list_weights_use.append(weights_use)

        hist_z, bins_z = np.histogram(vbin,bins=bins_photoz,weights=weights_use,normed=normed_hist)

        dx=0.01
        if plots: pl.plot(get_bins_centers(bins_z)+isb*dx,hist_z,'o',label=r'bin %d'%isb)       
        logger.info('bin %d n_gal %d' % (isb,nbin))

    if plots: pl.legend()
    if plots: pl.xlabel('photo-z')
    if plots: pl.ylabel('number of galaxies')
    if plots: pl.title('n(z) using weights')
    if plots: filename_fig = 'reweighted.%s.eps' % label
    if plots: logger.info('wrote %s' % (filename_fig))    

    return list_weights_use






def get_bins_centers(bins_edges,constant_spacing=True):

    # ensure np array
    bins_edges=np.array(bins_edges)
    bins_centers=np.zeros(len(bins_edges)-1)

    # get distance
    if constant_spacing:
        dx = bins_edges[1] - bins_edges[0]
        # print len(bins_edges)
        bins_centers = bins_edges[:-1] + dx/2.
    else:

        for be in range(len(bins_edges)-1):
            bins_centers[be] = np.mean([bins_edges[be],bins_edges[be+1]])      

    # print len(bins_centers)
    return bins_centers

   

