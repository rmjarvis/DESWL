import os
import yaml, argparse, sys, logging , pyfits
import numpy as np
import shutil, warnings; warnings.simplefilter("once")

logger = logging.getLogger("homogenise_nz")
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s ","%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.propagate = False

# global index of coadd objects ids for store
all_ids=None

def return_pdf_array(obj_ids,identifier,filename_photoz_h5,weight=None):
    """
    functions returns the stacked pfd of array of  coadd_objects_id's
    
    input:
    obj_ids : numpy type (int) array with coadd_objects_id
    identifier:
    string, either SKYNET,TPZ,ANNZ2,ZEBRA
    
    returns panda data frame with columsn 'z' and 'pdf'
    
    """
    import pandas as pd

    if identifier == 'SKYNET':
        z_max = 2.0
        nbins = 200
        key = 'SkyNet_v_1_2'
    elif identifier == 'SKYNET_MAG_AUTO':
        z_max = 2.0
        nbins = 200
        key = 'SkyNet_v_12_mag_auto'
    elif identifier == 'TPZ':
        z_max = 1.5
        nbins = 200
        key = 'TPZ'
    elif identifier == 'ANNZ':
        z_max = 1.5
        nbins = 150
        key = 'ANNZ'
    elif identifier == 'ZEBRA':
        z_max = 2.0
        nbins = 201
        key = 'ZEBRA'
    else:
        raise KeyError('Identifier must be in [SKYNET, SKYNET_MAG_AUTO, TPZ, ANNZ, ZEBRA]')
    
    pdf_array = _return_array(nbins, filename_photoz_h5, key, obj_ids,weight)
        
    if identifier == 'TPZ':
        dz = 0.007475
        z_edges = np.linspace(0.005-dz/2, 1.5-dz/2, 201)
    else:
        z_edges = np.linspace(0.0, z_max, nbins + 1)
    z_values = (z_edges[1:] + z_edges[:-1]) / 2.0
    
    return pdf_array, z_values
    
def _return_array(nbins, pdf_file, key, obj_ids,weight=None):
    """docstring for _return_stack"""

    import numpy as np

    # first use pandas to get the list of columns
    import pandas 
    store=pandas.HDFStore(pdf_file)
    row=store.select(key,start=0,stop=1)
    cols=[cn for cn in row.columns]
    cols.remove('z_mean')
    cols.remove('z_peak')
    n_cols = len(cols)
    logger.info('found %d columns',n_cols)
    sorting = np.argsort(cols)
    invsort = np.argsort(sorting)
    store.close()

    # we will use the h5py for slicing 
    # another wierd thing is that if I use h5py first, without opening and closing with pandas, it throws exception
    # apparently its a bug in current version

    import h5py
    store=h5py.File(pdf_file)

    # first create index 
    global all_ids
    if all_ids==None: 
        logger.info('loading indexing')
        all_ids = np.array(store[key]['table'][:]['index'])
    else:
        logger.info('using existing indexing')
    
    # sort coadd ids 
    all_index_sorting = np.argsort(all_ids)
    all_ids_sorted = all_ids[all_index_sorting]

    # sort the objects array
    obj_index_sorting = np.argsort(obj_ids)
    obj_index_invsorting = np.argsort(obj_index_sorting)
    obj_ids_sorted = obj_ids[obj_index_sorting]

    # match - get indices of coadd ids
    logger.debug('matching')
    obj_index_sorted = np.searchsorted(all_ids_sorted,obj_ids_sorted)
    obj_index = all_index_sorting[obj_index_sorted]
    n_gals_total = len(obj_index)

    # choose partitioning
    n_per_part = 100
    n_parts=n_gals_total/n_per_part+1

    # get container for parts
    list_obj_pdf = []

    # select data
    for ni in range(n_parts):
        istart = ni*n_per_part
        iend = (ni+1)*n_per_part
        obj_index_part = obj_index[istart:iend]
        obj_pdf = np.array(store[key]['table'][list(obj_index_part)]['values_block_0'])
        obj_pdf = obj_pdf[:,invsort]
        list_obj_pdf.append(obj_pdf)
        if ni % int(n_parts/10) == 0 : logger.info('read %7d/%7d objects' , ni*n_per_part, n_gals_total )

    # in place join parts
    obj_pdf = np.concatenate(list_obj_pdf,axis=0)

    # test
    # obj_pdf = all_ids[list(obj_index)]

    # in place copy to save memory
    obj_pdf = obj_pdf[obj_index_invsorting]

    # normalise and reweight 
    obj_pdf = obj_pdf * (weight/obj_pdf.sum(axis=1))[:,np.newaxis]

    logger.info('got array with %d galaxies, total memory size %2.3f GB' % (obj_pdf.shape[0],obj_pdf.nbytes/1.0e9))

    # close the sore
    store.close()

    # no point in upgrading to 32, waste of memory
    return obj_pdf.astype('float32')

def _normalize_pdf(pdf_all,dz):
    """
    returns normalized pdf
    """    
    area = np.trapz(pdf_all, dx=dz)
    return pdf_all/area


def get_weights(split_list_cat,target_nz_index=-1,photoz_min=0.2,photoz_max=1.2,photoz_nbins=50,normed_hist=True,label='systematics',plots=False):
    """
    @param split_list_cat list columns for redshifts and statistical weights of objects, each one corresponding to a bin split by systematics.
    Example split_list_cat = [ (cat_airmass_low_z,cat_airmass_low_w) , (cat_airmass_high_z,cat_airmass_high_w) ].
    For example cat_airmass = np.array(pyfits.getdata('DESXXXX-XXXX_cat.fits'))
    @param target_nz_index n(z) of this bin will be used as a target nz for reweighting of the other bins. If target_nz_index=-1 (default), then mean of all bins will be used.
    @param redshift_col each of these catalogs has to have a column corresponding to redshift point estimate.
    @param normed_hist if to use reweighting with normalised histograms (not much difference)
    @label title to pass for plots (example: 'airmass')
    @return return a list with weight vector for each systematic bin
    """
    if plots: import pylab as pl

    bins_photoz = np.linspace(photoz_min,photoz_max,photoz_nbins)

    n_bins = len(split_list_cat)
    list_hist_z = []


    logger.info('getting n(z) for %d systematic bins' % n_bins)
    for isb,vbin in enumerate(split_list_cat):

        col_redshift = vbin[0]
        col_statistical_weight = vbin[1]
        nbin = len(col_redshift)       

        hist_z, bins_z = np.histogram(col_redshift,bins=bins_photoz,normed=normed_hist,weights=col_statistical_weight)
        list_hist_z.append(hist_z)
        if plots:  pl.plot(get_bins_centers(bins_z),hist_z,label=r'bin %d '%(isb))

        logger.info('bin %d n_gal %d' % (isb,nbin))

    if target_nz_index==-1:
        mean_z = np.mean(np.array(list_hist_z),axis=0)
    else:
        try:
            mean_z = np.array(list_hist_z[target_nz_index])
        except:
            raise Exception('target_nz_index may be larger than the number of bins lin split_list_cat')


    if plots: pl.plot(get_bins_centers(bins_z),mean_z,label='target',lw=5,zorder = 0)
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

        col_redshift = vbin[0]
        col_statistical_weight = vbin[1]
        nbin = len(col_redshift)       

        hist_z, bins_z = np.histogram(col_redshift,bins=bins_photoz,normed=normed_hist,weights=col_statistical_weight)

        nz_weights = mean_z/hist_z

        list_weights.append(nz_weights)

        if plots: pl.plot(get_bins_centers(bins_z),nz_weights,label=r'bin %d'%(isb))
        
        logger.info('bin %d n_gal %d' % (isb,nbin))

    if plots: pl.xlabel('z')
    if plots: pl.ylabel('weight value')
    if plots: pl.title('w(z) values for %s bins' % label)
    if plots: pl.legend()
    if plots: filename_fig = 'wz_for_bins.%s.eps' % label
    if plots: pl.savefig(filename_fig)
    if plots: logger.info('wrote %s' % (filename_fig))


    if plots: pl.figure()
    if plots: pl.plot(get_bins_centers(bins_z),mean_z,label='target',lw=5,zorder = 0)
    digitize_bins = bins_photoz
    list_weights_use = []
    logger.info('getting reweighted histogrgams for %d systematic bins' % n_bins)
    for isb,vbin in enumerate(split_list_cat):

        col_redshift = vbin[0]
        col_statistical_weight = vbin[1]
        nbin = len(col_redshift)       

        inds = np.digitize(col_redshift, digitize_bins)

        weights_in_bin = list_weights[isb].copy()
        weights_in_bin = np.append(weights_in_bin,0)
        weights_in_bin = np.insert(weights_in_bin,0,0)
        weights_use = weights_in_bin[inds]
        list_weights_use.append(weights_use)

        hist_z, bins_z = np.histogram(col_redshift,bins=bins_photoz,weights=weights_use*col_statistical_weight,normed=normed_hist)

        dx=0.01
        if plots: pl.plot(get_bins_centers(bins_z)+isb*dx,hist_z,'o',label=r'bin %d'%isb)       
        logger.info('bin %d n_gal %d' % (isb,nbin))

    if plots: pl.legend()
    if plots: pl.xlabel('photo-z')
    if plots: pl.ylabel('number of galaxies')
    if plots: pl.title('n(z) using weights')
    if plots: filename_fig = 'reweighted.%s.eps' % label
    if plots: pl.savefig(filename_fig)
    if plots: logger.info('wrote %s - run pl.show() if needed' % (filename_fig))    

    return list_weights_use


def get_weights_fullPZ(split_list_cat,z_values,target_nz_index=-1,normed_hist=True,label='systematics',plots=False):
    """
    @param split_list_cat list columns for redshifts and statistical weights of objects, each one corresponding to a bin split by systematics.
    Example split_list_cat = [ (cat_airmass_low_z,cat_airmass_low_w) , (cat_airmass_high_z,cat_airmass_high_w) ].
    For example cat_airmass = np.array(pyfits.getdata('DESXXXX-XXXX_cat.fits'))
    @param target_nz_index n(z) of this bin will be used as a target nz for reweighting of the other bins. If target_nz_index=-1 (default), then mean of all bins will be used.
    @param redshift_col each of these catalogs has to have a column corresponding to redshift point estimate.
    @param normed_hist if to use reweighting with normalised histograms (not much difference)
    @label title to pass for plots (example: 'airmass')
    @return return a list with weight vector for each systematic bin
    """
    if plots: import pylab as pl

    n_bins = len(split_list_cat)
    list_hist_z = []

    logger.info('getting target n(z) for %d systematic bins' % n_bins)

    if target_nz_index==-1:

        raise Exception('not implemented yet target_nz_index==-1')

    elif target_nz_index in range(len(split_list_cat)):

        
        target_pz = split_list_cat[target_nz_index].sum(axis=0)
        normalisation = np.sum(target_pz)
        target_pz /= normalisation

        logger.info('n_gals in the target bin %d, normalisation %2.2f',split_list_cat[target_nz_index].shape[0],normalisation)

    else:
        raise Exception('wrong target_nz_index %d' , target_nz_index)
    
    if plots:  pl.plot(z_values,target_pz,'c',label='target',lw=5)

    list_bin_pz = [None]*n_bins
    for isb,vbin in enumerate(split_list_cat):

        list_bin_pz[isb] = split_list_cat[isb].sum(axis=0)
        list_bin_pz[isb] /= np.sum(list_bin_pz[isb])

        if plots:  pl.plot(z_values,list_bin_pz[isb],label='bin %d' % isb)

    if plots: pl.legend()
    if plots: pl.xlabel('z')
    if plots: pl.ylabel('stacked p(z)')
    if plots: pl.title('p(z) for %s bins' % label)
    if plots: filename_fig = 'nz_for_bins.%s.eps' % label
    if plots: pl.savefig(filename_fig)
    if plots: logger.info('wrote %s' % (filename_fig))

    if plots: pl.figure()
    if plots: pl.plot(z_values,target_pz,'c',label='target',lw=5)

    list_weights_use = []

    n_z_bins = split_list_cat[0].shape[1]

    sigma = 1e-5
    import scipy.sparse.linalg
    
    for isb,vbin in enumerate(split_list_cat):

        R = split_list_cat[isb] / np.sum(split_list_cat[isb].flatten())
        # t = R.sum(axis=0) - target_pz # wtf?
        t = target_pz - R.sum(axis=0)# wtf?
        # t = target_pz # wtf?


        # weight_pz, residuals, rank, singular_values = np.linalg.lstsq(R.T,target_pz)
        # np.dot( np.linalg.inv( (np.dot(R,R.T)+np.eye(R.shape[0])*sigma) ) )

        weight_pz, istop, itn, normr, normar, norma, conda, normx = scipy.sparse.linalg.lsmr(R.T, t, damp=sigma, atol=1e-06, btol=1e-06, conlim=100000000.0, maxiter=None, show=False)

        weight_pz += 1

        logger.info('getting weights for bin %d n_gals=%5d median_w=% 2.2e norm_w =% 2.2e' , isb, R.shape[0], np.median(weight_pz), normx )

        # weights_kron = np.kron(np.ones([1,n_z_bins]),weight_pz[:,None])
        weighted_stack  = np.dot(R.T,weight_pz)
        weighted_stack /= np.sum(weighted_stack)

        if plots:  pl.plot(z_values,weighted_stack,label=r'bin %d '%(isb))

        # pl.figure()
        # pl.plot(z_values,t,lw=5)
        # pl.plot(z_values,weighted_stack)

        list_weights_use.append(weight_pz)

    if plots: pl.legend()
    if plots: pl.xlabel('z')
    if plots: pl.ylabel('stacked p(z)')
    if plots: pl.title('p(z) for %s bins' % label)
    if plots: filename_fig = 'reweighted_pz.%s.eps' % label
    if plots: pl.savefig(filename_fig)
    if plots: logger.info('wrote %s' % (filename_fig))

    if plots: pl.figure()
    pl.hist(list_weights_use[0],200)

    if plots: pl.figure()
    pl.hist(list_weights_use[1],200)

    if plots: pl.figure()
    pl.hist(list_weights_use[2],200)

    pl.show()

    import pdb; pdb.set_trace()
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

   

