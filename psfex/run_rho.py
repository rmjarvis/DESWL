#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits*',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')

    # Options
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')

    args = parser.parse_args()
    return args


def measure_rho(ra,dec,e1,e2,s,m_e1,m_e2,m_s,max_sep):
    """Compute the rho statistics
    """
    import numpy
    import treecorr
    import galsim

    #print 'e1 = ',e1
    #print 'm_e1 = ',m_e1
    #print 'de1 = ',e1-m_e1
    #print 'mean e1 = ',numpy.mean(e1)
    #print 'mean de1 = ',numpy.mean(e1-m_e1)
    #print 'e2 = ',e2
    #print 'm_e2 = ',m_e2
    #print 'de2 = ',e2-m_e2
    #print 'mean e2 = ',numpy.mean(e2)
    #print 'mean de2 = ',numpy.mean(e2-m_e2)
    #print 's = ',s
    #print 'm_s = ',m_s
    #print 'ds = ',s-m_s
    #print 'mean s = ',numpy.mean(s)
    #print 'mean ds = ',numpy.mean(s-m_s)

    # From Barney's paper: http://arxiv.org/pdf/0904.3056v2.pdf
    # rho1 = < (e-em)* (e-em) >     Barney originally called this D1.
    # rho2 = Re < e* (e-em) >       Barney's D2 is actually 2x this.
    # rho3 = < (s^2-sm^2)/s^2 (s^2-sm^2)/s^2 >  Not in Barney's paper

    #print 'ra = ',ra
    #print 'dec = ',dec
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                            g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             g1=(e1-m_e1), g2=(e2-m_e2))
    dscat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             k=(s**2-m_s**2)/s**2)

    rho1 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho1.process(decat)
    #print 'rho1 = ',rho1.xip
    #print 'rho1.xip_im = ',rho1.xip_im
    #print 'rho1.xim = ',rho1.xim
    #print 'rho1.xim_im = ',rho1.xim_im
    #print 'rho1.sigma = ',numpy.sqrt(rho1.varxi)

    rho2 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho2.process(ecat, decat)
    #print 'rho2 = ',rho2.xip
    #print 'rho2.xip_im = ',rho2.xip_im
    #print 'rho2.xim = ',rho2.xim
    #print 'rho2.xim_im = ',rho2.xim_im
    #print 'rho2.sigma = ',numpy.sqrt(rho2.varxi)

    rho3 = treecorr.KKCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho3.process(dscat)
    #print 'rho3 = ',rho3.xi
    #print 'rho3.sigma = ',numpy.sqrt(rho3.varxi)

    return rho1,rho2,rho3


def main():
    import os
    import glob
    import galsim
    import json
    import numpy
    import astropy.io.fits as pyfits

    args = parse_args()

    # Make the work directory if it does not exist yet.
    work = os.path.expanduser(args.work)
    print 'work dir = ',work
    try:
        if not os.path.isdir(work):
            os.makedirs(work)
    except OSError as e:
        print "Ignore OSError from makedirs(work):"
        print e
        pass

    datadir = '/astro/u/astrodat/data/DES'

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        cat_file = os.path.join(exp_dir, exp + "_psf.fits")
        with pyfits.open(cat_file) as pyf:
            data = pyf[1].data

        ccdnums = numpy.unique(data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        stats = []

        for ccdnum in ccdnums:
            print '\nProcessing ', ccdnum

            mask = (data['ccdnum'] == ccdnum) & (data['flag'] == 0)
            if mask.sum() == 0:
                print '   All objects with this ccdnum are flagged.'
                continue

            ra = data['ra'][mask]
            dec = data['dec'][mask]
            e1 = data['e1'][mask]
            e2 = data['e2'][mask]
            size = data['size'][mask]

            psfex_e1 = data['psfex_e1'][mask]
            psfex_e2 = data['psfex_e2'][mask]
            psfex_size = data['psfex_size'][mask]

            rho1, rho2, rho3 = measure_rho(ra,dec,e1,e2,size,
                                           psfex_e1,psfex_e2,psfex_size,
                                           max_sep=20)

            k10arcmin = int(round(numpy.log(10 / 0.5)/0.1))
            #if numpy.abs(rho2.xip[k10arcmin]) > 1.e-4:
            if False:
                print '   rho2 at 10 arcmin = ',rho2.xip[k10arcmin]
                print '   flag this chip as bad in blacklist (TODO!)'
            else:
                stats.append([ 
                    int(ccdnum),
                    rho1.meanlogr.tolist(),
                    rho1.xip.tolist(),
                    rho1.xim.tolist(),
                    rho2.xip.tolist(),
                    rho2.xim.tolist(),
                    rho3.xi.tolist(),
                    ])
            print 'len stats = ',len(stats)

            if args.single_ccd:
                break

        mask = (data['flag'] == 0)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        ra = data['ra'][mask]
        dec = data['dec'][mask]
        e1 = data['e1'][mask]
        e2 = data['e2'][mask]
        size = data['size'][mask]

        psfex_e1 = data['psfex_e1'][mask]
        psfex_e2 = data['psfex_e2'][mask]
        psfex_size = data['psfex_size'][mask]

        rho1, rho2, rho3 = measure_rho(ra,dec,e1,e2,size,
                                       psfex_e1,psfex_e2,psfex_size,
                                       max_sep=100)

        desdm_e1 = data['desdm_e1'][mask]
        desdm_e2 = data['desdm_e2'][mask]
        desdm_size = data['desdm_size'][mask]

        drho1, drho2, drho3 = measure_rho(ra,dec,e1,e2,size,
                                          desdm_e1,desdm_e2,desdm_size,
                                          max_sep=100)

        #print 'rho1.xip = ',rho1.xip
        #print 'rho1.xip_im = ',rho1.xip_im
        #print 'rho1.xim = ',rho1.xim
        #print 'rho1.xim_im = ',rho1.xim_im
        #print 'rho2.xip = ',rho2.xip
        #print 'rho2.xip_im = ',rho2.xip_im
        #print 'rho2.xim = ',rho2.xim
        #print 'rho2.xim_im = ',rho2.xim_im
        #print 'rho3.xi = ',rho3.xi
        # Write out the interesting stats for this ccd into a file, which we can 
        # then pull all together into a single FITS catalog later.
        stat_file = os.path.join(exp_dir, exp + ".json")
        stats.append([
            int(expnum),
            rho1.meanlogr.tolist(),
            rho1.xip.tolist(),
            rho1.xip_im.tolist(),
            rho1.xim.tolist(),
            rho1.xim_im.tolist(),
            rho1.varxi.tolist(),
            rho2.xip.tolist(),
            rho2.xip_im.tolist(),
            rho2.xim.tolist(),
            rho2.xim_im.tolist(),
            rho2.varxi.tolist(),
            rho3.xi.tolist(),
            rho3.varxi.tolist(),
            drho1.meanlogr.tolist(),
            drho1.xip.tolist(),
            drho1.xip_im.tolist(),
            drho1.xim.tolist(),
            drho1.xim_im.tolist(),
            drho1.varxi.tolist(),
            drho2.xip.tolist(),
            drho2.xip_im.tolist(),
            drho2.xim.tolist(),
            drho2.xim_im.tolist(),
            drho2.varxi.tolist(),
            drho3.xi.tolist(),
            drho3.varxi.tolist(),
            ])
        print 'len stats = ',len(stats)
        with open(stat_file,'w') as f:
            json.dump(stats, f)

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
