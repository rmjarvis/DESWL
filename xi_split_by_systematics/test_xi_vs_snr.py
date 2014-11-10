import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ: mpl.use('tkagg')
import yaml, argparse, sys, logging , pyfits, fitsio
import numpy as np
import matplotlib.pyplot as pl; print 'using matplotlib backend' , pl.get_backend()
pl.rcParams['image.interpolation'] = 'nearest' ;
import shutil, plotstools, warnings; warnings.simplefilter("once")

logger = logging.getLogger("test_xi_vs_snr")
logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s ","%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.propagate = False

prob_z = None

def load_data():

    import glob;
    columns_list = [ config['columns'][name] for name in config['columns'].keys() ]
    filelist_results = glob.glob('%s/DES*.gz'%config['results_dir'])
    list_data = []
    logger.info('found %d files' % len(filelist_results))
    for file_results in filelist_results[args.first:(args.first+args.num)]:

        res=fitsio.read(file_results)
        warnings.warn('%s'%str(res.dtype.names))
        res=res[columns_list]
        res.dtype.names=config['columns'].keys()

        if config['method'] == 'im3shape':
            warnings.warn('using 1 as m1 and m2')
            res['m1'] = 1
            res['m2'] = 1

        res_cut = apply_cuts(res)
        logger.info('%s size=%d/%d' % (file_results,len(res_cut),len(res)) )
        list_data.append(res_cut)

    res = np.concatenate(list_data)
    logger.info('total data size=%d' % (len(res)))

    return res  
        
def apply_cuts(res):

    if config['method'] == 'im3shape':

        cut_vals = config['info_cuts']
        select = (res['error_flag'] == 0) 
        logger.debug('select on error_flag, catalog size %d' % (len(np.nonzero(select)[0]) ))

        for cv in cut_vals:

            select_cut =  (res['info_flag'] & cv) == 0
            select = select & select_cut
            logger.debug('select on %d, catalog size %d' % (cv,len(np.nonzero(select)[0]) ))

        res_cut = res[select]

    else:

        select = (res['error_flag'] == 0) 
        logger.debug('select on error_flag, catalog size %d' % (len(np.nonzero(select)[0]) ))
        res_cut = res[select]
    
    return res_cut

def get_weights():

    import homogenise_nz

    res = load_data()

    list_snr_bins_cats = []

    for isb,snr_bin in enumerate(config['snr_bins']):

        select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1]) 
        res_bin = res[select]
        res_bin_z = res_bin['photoz']

        list_snr_bins_cats.append(res_bin_z)

    list_weights = homogenise_nz.get_weights(list_snr_bins_cats,label='snr')

    import cPickle as pickle
    filename_pickle = 'weights.pp2'
    pickle.dump(list_weights,open(filename_pickle,'w'))
    print 'saved ', filename_pickle


    import pdb; pdb.set_trace()



def get_xi_vs_snr():

    res = load_data()

    n_photoz_bins = len(config['photoz_bins'])
    n_snr_bins = len(config['snr_bins'])

    list_bin_xi = []

    filename_pickle = 'weights.pp2'
    import cPickle as pickle; 
    list_weights = pickle.load(open(filename_pickle))

    for isb,snr_bin in enumerate(config['snr_bins']):

            select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1] )

            res_bin = res[select]
            nbin = len(res_bin)       

            weights_use = list_weights[isb]

            if len(weights_use) != len(res_bin):
            	raise Exception('different number of galaxies in weights and data - did you mess with the setting somewhere in the middle of the process?' % (len(list_weights),len(res_bin)))

            if args.no_weights:
                dummy_weights = np.ones_like(weights_use)
                logr,xip,xim,stdxi = get_corr(res_bin,dummy_weights)    
                filename_xi = 'list_bin_xi.noweights.pp2'
            else:
                logr,xip,xim,stdxi = get_corr(res_bin,weights_use)    
                filename_xi = 'list_bin_xi.pp2'

            hist_e1, bins_e1 = np.histogram(res_bin['e1'],bins=np.linspace(-1,1,100),normed=True)
            hist_e2, bins_e2 = np.histogram(res_bin['e2'],bins=np.linspace(-1,1,100),normed=True)
            list_bin_xi.append((logr,xip,xim,stdxi,nbin))

    import cPickle as pickle
    pickle.dump(list_bin_xi,open(filename_xi,'w'),protocol=2)
    logger.info('pickled %s',filename_xi)

def plot_xi_vs_snr():

    def _get_bin_centers(bins_edges):
        dx = bins_edges[1] - bins_edges[0]
        bins_centers = bins_edges[:-1] + dx/2.
        return bins_centers

    import cPickle as pickle

    if args.no_weights:
        filename_xi = 'list_bin_xi.noweights.pp2'
    else:
        filename_xi = 'list_bin_xi.pp2'

    list_bin_xi=pickle.load(open(filename_xi))

    n_photoz_bins = len(config['photoz_bins'])

    dx=0.05
    for isb,snr_bin in enumerate(config['snr_bins']):

            logr,xip,xim,stdxi,nbin = list_bin_xi[isb]

            pl.figure(1)
            pl.errorbar(np.exp(logr)*(isb*dx+1),xip/1e-5,yerr=stdxi/1e-5,label=r're[$\xi_+$] SNR$\in[%2.0f,%2.0f]$'%(snr_bin[0],snr_bin[1]),fmt='.')
            pl.axhline(0,c='k')
            pl.xscale('log',nonposx='clip'); 
            pl.yscale('symlog', nonposy='clip',linthreshy=1e-5); 
            pl.ylim([1e-3,1e2]);    
            pl.xlabel(r'$\theta$ [arcmin]');     
            pl.ylabel(r're[$\xi_{+}$] / $10^{-5}$') ;    
            pl.legend(loc='lower center',ncol=2,mode='expand') 
            # pl.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

            pl.figure(2)
            pl.errorbar(np.exp(logr)*(isb*dx+1),xim/1e-5,yerr=stdxi/1e-5,label=r're[$\xi_-$] SNR$\in[%2.0f,%2.0f]$'%(snr_bin[0],snr_bin[1]),fmt='.')
            pl.axhline(0,c='k')
            pl.xscale('log',nonposx='clip'); 
            pl.yscale('symlog', nonposy='clip',linthreshy=1e-5); 
            pl.ylim([1e-3,1e2]);    
            pl.xlabel(r'$\theta$ [arcmin]');     
            pl.ylabel(r're[$\xi_{-}$] / $10^{-5}$') ;    
            pl.legend(loc='lower center',ncol=2,mode='expand') 

            print 'SNR=%2.0f - %2.0f' % (snr_bin[0],snr_bin[1])
            for ir in range(len(logr)):
                print 'logr=% 2.2e \t xi+=% 2.4e \t xi-=% 2.4e \t +/- % 2.4e' % (np.exp(logr[ir]),xip[ir],xim[ir],stdxi[ir])

    
    # f.subplots_adjust(hspace=0,wspace=0)

    pl.figure()
    logr,xip1,xim,stdxi,nbin = list_bin_xi[0][1]
    logr,xip2,xim,stdxi,nbin = list_bin_xi[0][3]
    pl.plot(np.exp(logr),xip1/xip2,'-o')
    pl.xscale('log',nonposx='clip'); 
    pl.xlabel(r'$\theta$ [arcmin]');     
    pl.ylabel(r're[$\xi_{+}^{SNR=10-20} / \xi_{+}^{SNR=40-}$]') ;
    pl.title('SNR=10-20 vs SNR=40- ')    

    
    pl.show()


    import pdb; pdb.set_trace()

def get_corr(res,weights):

    import treecorr
    nbins=config['treecorr_n_bins']
    min_sep=config['treecorr_min_sep']
    max_sep=config['treecorr_max_sep']
    sep_units='arcmin'

    shape_cat = treecorr.Catalog(g1=res['e1'], g2=res['e2'], ra= res['ra'], dec= res['de'], ra_units='deg', dec_units='deg' , w = weights)
    gg = treecorr.GGCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units,verbose=2)  
    gg.process(shape_cat)

    responsivity_cat = treecorr.Catalog(g1=res['m1'], g2=res['m2'], ra= res['ra'], dec= res['de'], ra_units='deg', dec_units='deg' , w = weights)
    mm = treecorr.GGCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units,verbose=2)  
    mm.process(responsivity_cat)

    if args.use_responsivity==True:
        xip = gg.xip / mm.xip*2.            # The real part of xi+
        xim = gg.xim / mm.xim*2.            # The real part of xi-
        warnings.warn('using responsivity')
    else:
        xip = gg.xip             # The real part of xi+
        xim = gg.xim             # The real part of xi-

    logr = gg.logr          # The nominal center of each bin
    meanlogr = gg.meanlogr  # The mean <log(r)> within the bins
    varxi = gg.varxi        # The variance of each xi+ or xi- value taking into account shape noise only
    stdxi = np.sqrt(varxi)

    return logr,xip,xim,stdxi


def main():

    valid_actions = ['plot_xi_vs_snr','get_xi_vs_snr','get_weights']

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
