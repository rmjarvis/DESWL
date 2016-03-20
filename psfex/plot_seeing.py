#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
import numpy
 
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


def add_to_list(filter, vlist, value):
    if 'griz' in vlist.keys():
        vlist['griz'].append(value)
    if filter not in vlist.keys():
        vlist[filter] = []
    if 'riz' in vlist.keys() and 'g' not in filter:
        vlist['riz'].append(value)
    vlist[filter].append(value)

def get_data(runs, exps, work, fwhm_list, med_size_list, mean_size_list, tag):

    import astropy.io.fits as pyfits
    import os

    expinfo_file = '/astro/u/mjarvis/work/exposure_info_' + tag + '.fits'
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    cat_dir = os.path.join(work,'psf_cats')

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        #print 'expnum = ',expnum

        if expnum not in expinfo['expnum']:
            print 'expnum is not in expinfo!'
            print 'expinfo[expnum] = ',expinfo['expnum']
            raise RuntimeError("Could not find information about this expnum")
        k = numpy.nonzero(expinfo['expnum'] == expnum)[0][0]
        #print 'k = ',k
        if expinfo['flag'][k] != 0:
            print 'exposure is flagged.'
            continue

        fwhm = expinfo['fwhm'][k]
        # This is in pixels.  Convert to arcsec:
        fwhm *= 0.263

        #if fwhm > 1.5:
            #print 'fwhm = %s > 1.5'%fwhm
            #continue

        cat_file = os.path.join(cat_dir, exp + "_psf.fits")
        try:
            with pyfits.open(cat_file) as pyf:
                data = pyf[1].data.copy()
        except:
            print 'Could not open cat_file %s.  Skipping this one.'%cat_file
            continue

        mask = (data['flag'] == 0)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            continue

        size = data['size'][mask] * 2.3548

        filter = expinfo['filter'][k]
        print 'filter = ',filter

        add_to_list(filter, fwhm_list, fwhm)
        add_to_list(filter, med_size_list, numpy.median(size))
        add_to_list(filter, mean_size_list, numpy.mean(size))

    print '\nFinished processing all exposures'


def plot_seeing(fwhm, tag=None):
    fig = plt.figure()
    plt.minorticks_on()
    
    ax = fig.add_subplot(111)

    print 'median seeing in g = ',numpy.median(fwhm['g'])
    print 'median seeing in r = ',numpy.median(fwhm['r'])
    print 'median seeing in i = ',numpy.median(fwhm['i'])
    print 'median seeing in z = ',numpy.median(fwhm['z'])

    riz = numpy.concatenate([fwhm['r'], fwhm['i'], fwhm['z']])
    print 'median seeing in riz = ',numpy.median(riz)

    nbins = 40
    range = (0.7, 1.7)
    #n, bins, p = ax.hist(riz, bins=nbins, range=range, histtype='step', fill=True,
                         #color='black', facecolor='cyan', label='riz')
    #n, bins, p = ax.hist(fwhm['g'], bins=bins, histtype='step', color='green', label='g')
    #n, bins, p = ax.hist(fwhm['r'], bins=bins, histtype='step', color='red', label='r')
    #n, bins, p = ax.hist(fwhm['i'], bins=bins, histtype='step', color='magenta', label='i')
    #n, bins, p = ax.hist(fwhm['z'], bins=bins, histtype='step', color='blue', label='z')
    width = (range[1]-range[0])/nbins
    n, bins, p = ax.hist([fwhm['z'],fwhm['i'],fwhm['r']], bins=nbins, range=range, 
                         histtype='barstacked', fill=True,
                         color=['black','purple','red'], width=width)
    ax.set_xlabel('Seeing FWHM (arcsec)')
    ax.set_ylabel('Number of exposures')
    ax.legend(reversed(p), ['r', 'i', 'z'], loc='upper right')
    ax.set_xlim(*range)
    plt.tight_layout()
    if tag is None:
        plt.savefig('seeing.pdf')
    else:
        plt.savefig('seeing_%s.pdf'%tag)
 

def main():
    import os

    args = parse_args()

    work = os.path.expanduser(args.work)
    print 'work dir = ',work

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    fwhm_list = {}
    med_size_list = {}
    mean_size_list = {}

    get_data(runs, exps, work, fwhm_list, med_size_list, mean_size_list, args.tag)

    print 'FWHM from FITS headers:'
    plot_seeing(fwhm_list)
    print 'Median seeing based on hsm measured size'
    plot_seeing(med_size_list, 'med')
    print 'Mean seeing based on hsm measured size'
    plot_seeing(mean_size_list, 'mean')


if __name__ == "__main__":
    main()
