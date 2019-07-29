#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

 
import os
import glob
import galsim
import json
import numpy
import astropy.io.fits as pyfits

def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')
    parser.add_argument('--output_dir', default=None,
                        help='location of output directory (default: {work}/{exp}/)')

    # Exposure inputs
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


def measure_rho(ra,dec,e1,e2,s,m_e1,m_e2,m_s,max_sep, tag=None):
    """Compute the rho statistics
    """
    import treecorr

    de1 = e1-m_e1
    de2 = e2-m_e2
    dt = (s**2-m_s**2)/s**2

    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             k=dt, g1=dt*e1, g2=dt*e2)

    ecat.name = 'ecat'
    decat.name = 'decat'
    dtcat.name = 'dtcat'
    if tag is not None:
        for cat in [ ecat, decat, dtcat ]:
            cat.name = tag + ":"  + cat.name

    min_sep = 0.5
    bin_size = 0.5
    bin_slop = 0.1

    results = []
    for (cat1, cat2) in [ (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:
        print 'Doing correlation of %s vs %s'%(cat1.name, cat2.name)

        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        results.append(rho)

    return results


def main():

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
        try:
            expnum = int(exp[6:])
        except:
            expnum = 0
        print 'expnum = ',expnum

        if args.output_dir is None:
            exp_dir = os.path.join(work,exp)
        else:
            exp_dir = args.output_dir
        print 'exp_dir = ',exp_dir

        cat_dir = os.path.join(work,'psf_cats')
        print os.path.join(cat_dir, '*%s*_psf.fits'%exp)
        cat_files = glob.glob(os.path.join(cat_dir, '*%s*_psf.fits'%exp))
        if len(cat_files) == 1:
            print 'cat_file = ',cat_files[0]
            with pyfits.open(cat_files[0]) as pyf:
                data = pyf[1].data
            all_data = [data]
        else:
            print 'cat_files = ',cat_files
            all_data = []
            for cat_file in cat_files:
                with pyfits.open(cat_file) as pyf:
                    all_data.append(pyf[1].data.copy())
                print 'len(all_data) = ',len(all_data)
                print 'last shape = ',all_data[-1].shape
            data = numpy.concatenate(all_data)
        print 'data.shape = ',data.shape

        ccdnums = numpy.unique(data['ccdnum'])
        print 'ccdnums = ',ccdnums

        stats = []

        for ccdnum in ccdnums:
            print '\nProcessing ', ccdnum

            mask = (data['ccdnum'] == ccdnum) & (data['flag'] == 0)
            print 'len(mask) = ',len(mask)
            print 'nonzero(mask) = ',numpy.sum(mask)
            if mask.sum() == 0:
                print '   All objects with this ccdnum are flagged.'
                continue

            ra = data['ra'][mask]
            dec = data['dec'][mask]
            e1 = data['e1'][mask]
            e2 = data['e2'][mask]
            size = data['size'][mask]

            psf_e1 = data['psf_e1'][mask]
            psf_e2 = data['psf_e2'][mask]
            psf_size = data['psf_size'][mask]

            rho1, rho2, rho3, rho4, rho5 = measure_rho(ra,dec,e1,e2,size,
                                                       psf_e1,psf_e2,psf_size,
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
                    rho3.xip.tolist(),
                    rho3.xim.tolist(),
                    rho4.xip.tolist(),
                    rho4.xim.tolist(),
                    rho5.xip.tolist(),
                    rho5.xim.tolist(),
                    ])
            print 'len stats = ',len(stats)

            if args.single_ccd:
                break

        print 'Do overall rho stats'
        mask = (data['flag'] == 0)
        print 'len(mask) = ',len(mask)
        print 'nonzero(mask) = ',numpy.sum(mask)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        ra = data['ra'][mask]
        dec = data['dec'][mask]
        e1 = data['e1'][mask]
        e2 = data['e2'][mask]
        size = data['size'][mask]

        psf_e1 = data['psf_e1'][mask]
        psf_e2 = data['psf_e2'][mask]
        psf_size = data['psf_size'][mask]

        rho1, rho2, rho3, rho4, rho5 = measure_rho(ra,dec,e1,e2,size,
                                                   psf_e1,psf_e2,psf_size,
                                                   max_sep=100)

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
            rho1.varxip.tolist(),
            rho1.varxim.tolist(),
            rho2.xip.tolist(),
            rho2.xip_im.tolist(),
            rho2.xim.tolist(),
            rho2.xim_im.tolist(),
            rho2.varxip.tolist(),
            rho2.varxim.tolist(),
            rho3.xip.tolist(),
            rho3.xip_im.tolist(),
            rho3.xim.tolist(),
            rho3.xim_im.tolist(),
            rho3.varxip.tolist(),
            rho3.varxim.tolist(),
            rho4.xip.tolist(),
            rho4.xip_im.tolist(),
            rho4.xim.tolist(),
            rho4.xim_im.tolist(),
            rho4.varxip.tolist(),
            rho4.varxim.tolist(),
            rho5.xip.tolist(),
            rho5.xip_im.tolist(),
            rho5.xim.tolist(),
            rho5.xim_im.tolist(),
            rho5.varxip.tolist(),
            rho5.varxim.tolist(),
            ])
        print 'len stats = ',len(stats)
        with open(stat_file,'w') as f:
            json.dump(stats, f)

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
