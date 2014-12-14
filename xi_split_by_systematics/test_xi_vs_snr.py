import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ: mpl.use('tkagg')
import yaml, argparse, sys, logging , pyfits, fitsio
import numpy as np
import matplotlib.pyplot as pl; print 'using matplotlib backend' , pl.get_backend()
pl.rcParams['image.interpolation'] = 'nearest' ;
import shutil, tktools, warnings; warnings.simplefilter("once")

logger = logging.getLogger("test_xi_vs_snr")
logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s ","%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.propagate = False

prob_z = None

def add_col(rec, name, arr, dtype=None):

    import numpy
    arr = numpy.asarray(arr)
    if dtype is None:
        dtype = arr.dtype
    newdtype = numpy.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = numpy.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec

def load_data():

    import glob;
    columns_list = [ config['columns'][name] for name in config['columns'].keys() ]
    filelist_results = glob.glob(config['results_files'])
    list_data = []
    logger.info('found %d files' % len(filelist_results))
    for file_results in filelist_results[args.first:(args.first+args.num)]:

        res=fitsio.read(file_results)
        warnings.warn('%s'%str(res.dtype.names))
        res=res[columns_list]
        for col in columns_list:
            if col not in res.dtype.names:
                logger.error('column %s not found in file %s',col,file_results)

        res.dtype.names=config['columns'].keys()

        if config['method'] == 'im3shape':
            warnings.warn('using 1 as m1 and m2')
            res = add_col(res,'m1',res['m'])
            res = add_col(res,'m2',res['m'])
        elif config['method'] == 'ngmix':
            warnings.warn('using m=mean(m1,m2)')
            m=(res['m1']+res['m2'])/2. -1
            res['e1'] = -res['e1'].copy()
            res = add_col(res,'m',m)
            res = add_col(res,'c1',np.zeros(len(res)))
            res = add_col(res,'c2',np.zeros(len(res)))
            w=1./(2.*0.16**2.+res['cov11']+res['cov22']+2.*res['cov12'])
            res = add_col(res,'w', w)
        else:
            raise Exception('unknown shear catalogs %s',config['method'])

        res_cut = apply_cuts(res)
        logger.info('%s size=%d/%d' % (file_results,len(res_cut),len(res)) )
        list_data.append(res_cut)

    res = np.concatenate(list_data)
    logger.info('total data size=%d' % (len(res)))

    return res  
        
def apply_cuts(res):

    if config['method'] == 'im3shape':

        # cut_vals = config['info_cuts']
        select = (res['error_flag'] == 0) & (res['info_flag'] == 0) & (res['ra']>56) & (res['ra']<94) & (res['de']<-42) & (res['sg']==1)
        logger.debug('cuts done, catalog size %d' % (len(np.nonzero(select)[0]) ))
        res_cut = res[select]

    elif config['method'] == 'ngmix':    

        select = (res['flags'] == 0) & (res['arate'] > 0.4) & (res['arate'] < 0.6) & (res['snrt']>3) & (res['ra']>56) & (res['ra']<94) & (res['de']<-42) & (res['sg']==1)
        logger.debug('cuts done, catalog size %d' % (len(np.nonzero(select)[0]) ))
        res_cut = res[select]

    else:
        raise Exception('unknown shear catalogs %s',config['method'])


   
    return res_cut

def get_weights():

    import homogenise_nz

    res = load_data()

    list_snr_bins_cats = []

    for isb,snr_bin in enumerate(config['snr_bins']):

        select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1]) 
        res_bin = res[select]
        res_bin_z = (res_bin['photoz'], res_bin['w'])

        list_snr_bins_cats.append(res_bin_z)

    label = 'snr.%s' % config['method']
    list_weights = homogenise_nz.get_weights(list_snr_bins_cats,target_nz_index=0,label=label,photoz_min=0.3,photoz_max=1.3,plots=True)


    import cPickle as pickle
    filename_pickle = 'weights.%s.pp2' % (config['method'])
    pickle.dump(list_weights,open(filename_pickle,'w'))
    print 'saved ', filename_pickle

    import pdb; pdb.set_trace()
    pl.show()

def get_weights_fullPZ():

    import homogenise_nz
 
    res = load_data()

    list_snr_bins_cats = []

    for isb,snr_bin in enumerate(config['snr_bins']):

        select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1]) 
        res_bin = res[select]
      
        import DES_pdf_stacker
        pdf_array, z_values = DES_pdf_stacker.return_pdf_array(res_bin['coadd_id'],config['name_pzcode'],config['filename_photoz_h5'],weight=res_bin['w'])
        # pdf_array, z_values = homogenise_nz.return_pdf_array(res_bin['coadd_id'],config['name_pzcode'],config['filename_photoz_h5'],weight=res_bin['w'])

        list_snr_bins_cats.append(pdf_array)

    label = 'snr.%s.%s' % (config['method'],config['name_pzcode'])
    list_weights = homogenise_nz.get_weights_fullPZ(list_snr_bins_cats,z_values=z_values,target_nz_index=0,label=label,plots=True,sigma_regularisation = 1e-9)

    import cPickle as pickle
    filename_pickle = 'weights.%s.pp2' % (config['method'])
    pickle.dump(list_weights,open(filename_pickle,'w'))
    print 'saved ', filename_pickle

    import pdb; pdb.set_trace()
    pl.show()

def get_xi_vs_snr():

    res = load_data()

    n_snr_bins = len(config['snr_bins'])

    list_bin_xi = []

    filename_pickle = 'weights.%s.pp2' % (config['method'])
    import cPickle as pickle; 
    list_weights = pickle.load(open(filename_pickle))

    for isb,snr_bin in enumerate(config['snr_bins']):

            select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1] )

            res_bin = res[select]
            nbin = len(res_bin)       

            weights_use = list_weights[isb]

            if len(weights_use) != len(res_bin):
                raise Exception('different number of galaxies in weights %d and data %d - did you mess with settings somewhere in the middle of the process?' % (len(weights_use),len(res_bin)))


            filename_xi = 'list_bin_xi.%s.' % config['method']
            if args.no_weights:
                weights_use = np.ones_like(weights_use)
                filename_xi = filename_xi + 'no-homogen.'
            else:
                filename_xi = filename_xi + 'homogen.'
            if args.use_responsivity:
                filename_xi = filename_xi + 'resp.'
            else:
                filename_xi = filename_xi + 'no-resp.'            
            filename_xi = filename_xi + 'pp2'

            logr,xip,xim,stdxi,mean_responsivity = get_corr(res_bin,weights_use)    

            hist_e1, bins_e1 = np.histogram(res_bin['e1'],bins=np.linspace(-1,1,100),normed=True)
            hist_e2, bins_e2 = np.histogram(res_bin['e2'],bins=np.linspace(-1,1,100),normed=True)
            list_bin_xi.append((logr,xip,xim,stdxi,nbin,mean_responsivity))

    import cPickle as pickle
    pickle.dump(list_bin_xi,open(filename_xi,'w'),protocol=2)
    logger.info('pickled %s',filename_xi)

def plot_xi_vs_snr():

    def _get_bin_centers(bins_edges):
        dx = bins_edges[1] - bins_edges[0]
        bins_centers = bins_edges[:-1] + dx/2.
        return bins_centers

    import cPickle as pickle

    titlestr = config['method'] + ': '

    filename_xi = 'list_bin_xi.%s.' % config['method']
    if args.no_weights:
        filename_xi = filename_xi + 'no-homogen.'
        titlestr += 'no n(z) reweighting '
    else:
        filename_xi = filename_xi + 'homogen.'
        titlestr += 'using n(z) reweighting '

    if args.use_responsivity:
        filename_xi = filename_xi + 'resp.'
        titlestr += 'and using responsivity/nbc'
    else:
        filename_xi = filename_xi + 'no-resp.'
        titlestr += 'and no responsivity/nbc'
    filename_xi = filename_xi + 'pp2'

    list_bin_xi=pickle.load(open(filename_xi))

    dx=0.05
    for isb,snr_bin in enumerate(config['snr_bins']):

            logr,xip,xim,stdxi,nbin,mean_responsivity = list_bin_xi[isb]

            pl.figure(1)
            pl.errorbar(np.exp(logr)*(isb*dx+1),xip/1e-5,yerr=stdxi/1e-5,label=r're[$\xi_+$] SNR$\in[%2.0f,%2.0f]$'%(snr_bin[0],snr_bin[1]),fmt='.')
            pl.axhline(0,c='k')
            pl.xscale('log',nonposx='clip'); 
            pl.yscale('symlog', nonposy='clip',linthreshy=1e-5); 
            pl.ylim([1e-3,1e2]);    
            pl.xlabel(r'$\theta$ [arcmin]');     
            pl.ylabel(r're[$\xi_{+}$] / $10^{-5}$') ;    
            pl.legend(framealpha=0.0,frameon=False, loc='lower center')
            # pl.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
            pl.title(titlestr)

            pl.figure(2)
            pl.errorbar(np.exp(logr)*(isb*dx+1),xim/1e-5,yerr=stdxi/1e-5,label=r're[$\xi_-$] SNR$\in[%2.0f,%2.0f]$'%(snr_bin[0],snr_bin[1]),fmt='.')
            pl.axhline(0,c='k')
            pl.xscale('log',nonposx='clip'); 
            pl.yscale('symlog', nonposy='clip',linthreshy=1e-5); 
            pl.ylim([1e-3,1e2]);    
            pl.xlabel(r'$\theta$ [arcmin]');     
            pl.ylabel(r're[$\xi_{-}$] / $10^{-5}$') ;    
            pl.legend(framealpha=0.0,frameon=False, loc='lower center')
            pl.title(titlestr)

            pl.figure(3)
            pl.errorbar(np.exp(logr)*(isb*dx+1),mean_responsivity,label=r'responsivity SNR$\in[%2.0f,%2.0f]$'%(snr_bin[0],snr_bin[1]),fmt='.')
            pl.axhline(0,c='k')
            pl.xscale('log',nonposx='clip'); 
            pl.yscale('symlog', nonposy='clip',linthreshy=1e-5); 
            pl.ylim([1e-3,1e2]);    
            pl.xlabel(r'$\theta$ [arcmin]');     
            pl.ylabel(r're[$\xi_{-}$] / $10^{-5}$') ;    
            pl.legend(framealpha=0.0,frameon=False, loc='lower center')
            pl.title(titlestr)


            print 'SNR=%2.0f - %2.0f' % (snr_bin[0],snr_bin[1])
            for ir in range(len(logr)):
                print 'logr=% 2.2e \t xi+=% 2.4e \t xi-=% 2.4e \t +/- % 2.4e' % (np.exp(logr[ir]),xip[ir],xim[ir],stdxi[ir])

    
    # f.subplots_adjust(hspace=0,wspace=0)

    # pl.figure()
    # logr,xip1,xim,stdxi,nbin = list_bin_xi[0][1]
    # logr,xip2,xim,stdxi,nbin = list_bin_xi[0][3]
    # pl.plot(np.exp(logr),xip1/xip2,'-o')
    # pl.xscale('log',nonposx='clip'); 
    # pl.xlabel(r'$\theta$ [arcmin]');     
    # pl.ylabel(r're[$\xi_{+}^{SNR=10-20} / \xi_{+}^{SNR=40-}$]') ;
    # pl.title('SNR=10-20 vs SNR=40- ')    

    
    pl.show()


    import pdb; pdb.set_trace()

def get_corr(res,weights):

    import treecorr
    nbins=config['treecorr_n_bins']
    min_sep=config['treecorr_min_sep']
    max_sep=config['treecorr_max_sep']
    sep_units='arcmin'

    shape_cat = treecorr.Catalog(g1=res['e1']-res['c1'], g2=res['e2']-res['c2'], ra= res['ra'], dec= res['de'], ra_units='deg', dec_units='deg' , w = weights * res['w'])
    gg = treecorr.GGCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units,verbose=2)  
    gg.process(shape_cat)

    responsivity_type = 'scalar'

    if responsivity_type == 'scalar': 
        responsivity_cat = treecorr.Catalog(k=(1+res['m']), ra= res['ra'], dec= res['de'], ra_units='deg', dec_units='deg' , w = weights * res['w'])
        mm = treecorr.KKCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units,verbose=2)  
        mm.process(responsivity_cat)
        corr_xip = mm.xi
        corr_xim = mm.xi
    elif responsivity_type == 'vector':
        responsivity_cat = treecorr.Catalog(g1=(1+res['m1']), g2=(1+res['m2']), ra= res['ra'], dec= res['de'], ra_units='deg', dec_units='deg' , w = weights * res['w'])
        mm = treecorr.GGCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units,verbose=2)  
        mm.process(responsivity_cat)
        corr_xip = mm.xip
        corr_xim = mm.xim

    if args.use_responsivity==True:
        xip = gg.xip / corr_xip            # The real part of xi+
        xim = gg.xim / corr_xim         # The real part of xi-
        warnings.warn('using responsivity ')
        logger.info('responsivity %s',str(corr_xip))
    else:
        xip = gg.xip             # The real part of xi+
        xim = gg.xim             # The real part of xi-

    logr = gg.logr          # The nominal center of each bin
    meanlogr = gg.meanlogr  # The mean <log(r)> within the bins
    varxi = gg.varxi        # The variance of each xi+ or xi- value taking into account shape noise only
    stdxi = np.sqrt(varxi)
    mean_responsivity = mm.xi

    return logr,xip,xim,stdxi,mean_responsivity


def main():

    valid_actions = ['plot_xi_vs_snr','get_xi_vs_snr','get_weights','get_weights_fullPZ']

    description = 'test_xi_vs_snr'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', type=str, default='config.yaml' , action='store', help='filename of file containing config')
    parser.add_argument('-n', '--num',type=int,  default=2, action='store', help='number of results files to use')
    parser.add_argument('-f', '--first',type=int, default=0, action='store', help='number of first result file to open')
    parser.add_argument('-a','--actions', nargs='+', action='store', help='which actions to run, available: %s' % str(valid_actions) )
    parser.add_argument('-rd','--results_dir', action='store', help='where results files are' , default='results/' )
    parser.add_argument('--no_weights', action='store_true', help='if to use redshift weights' , default=False )
    parser.add_argument('--use_responsivity', action='store_true', help='if to use responsivity' , default=False )

    global args

    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logger.setLevel(logging_level)

    global config
    config = yaml.load(open(args.filename_config))
    if args.actions==None:
        logger.error('no action specified, choose from %s and use -a {action}' % valid_actions)
        return
    for action in valid_actions:
        if action in args.actions:
            logger.info('executing %s' % action)
            exec action+'()'
    for ac in args.actions:
        if ac not in valid_actions:
            print 'invalid action %s' % ac

if __name__=='__main__':
    main()
