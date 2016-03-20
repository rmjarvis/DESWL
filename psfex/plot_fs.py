#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

import matplotlib
matplotlib.use('Agg') # needs to be done before import pyplot
import matplotlib.pyplot as plt
plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')

# Overall zeropoint adjustment from Eli.
zp = 5.3

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

def plot_fs(size, mag, flag, filename, modest=None):
    import numpy
    fig, ax = plt.subplots(1,1)

    # Mag adjustment from Eli.
    mag = mag + zp

    ax.set_xlim(10+zp, 24.)
    ax.set_ylim(0., 0.8)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel(r'$T = 2\sigma^2$ (arcsec${}^2$)')

    index = numpy.argsort(mag)
    n = numpy.arange(len(mag))

    for i in range(20):
        s = index[n%20==i]

        if modest is not None:
            mod = ax.scatter(mag[s][modest[s]&(flag[s]<4)], size[s][modest[s]&(flag[s]<4)],
                             facecolors='none', edgecolors='blue', s=50., lw=0.8, alpha=0.3)

        det = ax.scatter(mag[s][flag[s] == 0], size[s][flag[s] == 0],
                         color='black', marker='o', s=2.)
        cand = ax.scatter(mag[s][flag[s] == 1], size[s][flag[s] == 1],
                          color='magenta', marker='o', s=8.)
        psf = ax.scatter(mag[s][flag[s] == 2], size[s][flag[s] == 2],
                         color='green', marker='*', s=40.)

    ax.legend([det, cand, psf, mod],
              ['Detected Object', 'Candidate Star', 'PSF Star', 'Modest Star'],
              loc='upper left', frameon=True)

    plt.tight_layout()
    plt.savefig(filename)


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

    cat_dir = os.path.join(work,'psf_cats')

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))

        psf_file = os.path.join(cat_dir, exp + "_psf.fits")
        with pyfits.open(psf_file) as pyf:
            psf_data = pyf[1].data

        ccdnums = numpy.unique(psf_data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        stats = []

        for ccdnum in ccdnums:
            print '\nProcessing ', ccdnum

            mask = (psf_data['ccdnum'] == ccdnum)
            if mask.sum() == 0:
                print '   All objects with this ccdnum are flagged.'
                continue

            sex_file = os.path.join(input_dir, exp + "_%02d_cat.fits"%ccdnum)
            with pyfits.open(sex_file) as pyf:
                sex_data = pyf[2].data

            cat_file = os.path.join(exp_dir, exp + "_%02d_findstars.fits"%ccdnum)
            with pyfits.open(cat_file) as pyf:
                cat_data = pyf[1].data

            mag = cat_data['mag']
            size = 2. * cat_data['sigma0']**2
            # Flag definitions:
            # 0 = good detection
            # 1 = star candidate
            # 2 = final star
            # 4 = bad measurement
            flag = cat_data['star_flag']
            print 'ntot = ',len(cat_data)
            print 'nstar = ',numpy.sum(flag==1)
            flag[cat_data['size_flags']>0] = 4
            print 'nbad = ',numpy.sum(flag==4)

            x = cat_data['x']
            y = cat_data['y']

            psfex_file = os.path.join(exp_dir, exp + "_%02d_psfcat.used.fits"%ccdnum)
            with pyfits.open(psfex_file) as pyf:
                psfex_data = pyf[2].data

            xu = psfex_data['X_IMAGE']
            yu = psfex_data['Y_IMAGE']
            used = [ numpy.where((numpy.abs(xxu-x)<1) & (numpy.abs(yyu-y)<1))[0][0] for (xxu,yyu) in zip(xu,yu) ]
            flag[used] = 2
            print 'nused = ',numpy.sum(flag==2)

            xs = sex_data['X_IMAGE']
            ys = sex_data['Y_IMAGE']
            # Matching this is a bit tricky, since they are different detections, measurements
            # Just pick out the items with a single match within 1 pixel.
            m1 = [ numpy.where((numpy.abs(xxs-x)<1) & (numpy.abs(yys-y)<1))[0] for (xxs,yys) in zip(xs,ys) ]
            m2 = numpy.array([ len(a) for a in m1 ])
            cat_m = numpy.concatenate(numpy.array(m1)[m2==1])
            # Now cat_data[cat_m] gives rows that are in the original sextractor catalog.
            sex_m = numpy.where(m2==1)
            # And sex_data[sex_m] gives the corresponding rows in the sextractor catalog.
            modest = numpy.zeros_like(size, dtype=bool)
            bright_test = (sex_data[sex_m]['CLASS_STAR'] > 0.3) & (sex_data[sex_m]['MAG_AUTO'] < 18.0-zp)
            locus_test = (sex_data[sex_m]['SPREAD_MODEL'] + 3.*sex_data[sex_m]['SPREADERR_MODEL'] < 0.003)
            faint_psf_test = (sex_data[sex_m]['MAG_PSF'] > 30.-zp) & (sex_data[sex_m]['MAG_AUTO'] < 21.-zp)
            modest[cat_m] = (bright_test | locus_test) & ~faint_psf_test

            plot_file = os.path.join(exp_dir, exp + "_%02d.pdf"%ccdnum)
            plot_fs(size, mag, flag, plot_file, modest)

            if args.single_ccd:
                break

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
