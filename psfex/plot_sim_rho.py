#! /usr/bin/env python
# Program to plot rho statistics on PSFEx outputs.

import astropy.io.fits as pyfits
import numpy
import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os
import sys
import glob
import treecorr

plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot rho stats for sim')

    # Drectory arguments
    parser.add_argument('--psf_cats', default='./',
                        help='location of psf_cats directory')

    args = parser.parse_args()
    return args


def pretty_rho1(meanr, rho, sig, rho3=None, rho4=None):
    import matplotlib.patches as mp
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig[rho3>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig[rho3<0], color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig, color='green', marker='s')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig[rho4>0], color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig[rho4<0], color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig, color='red', marker='^')
    if rho3 is not None and rho4 is not None:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$', r'$\rho_4(\theta)$'],
                   loc='upper right', fontsize=24)
        plt.ylim( [1.e-9, 5.e-6] )
    else:
        plt.legend([rho1_line],
                   [r'$\rho_1(\theta)$'],
                   loc='upper right')
        plt.ylim( [1.e-9, 5.e-6] )
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\rho(\theta)$')
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, rho5=None):
    import matplotlib.patches as mp
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho2_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='green')
        plt.plot(meanr*1.03, -rho5, color='green', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig[rho5>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig[rho5<0], color='green', ls='', marker='s')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig, color='green', marker='s')
    if rho5 is not None:
        plt.legend([rho2_line, rho5_line],
                   [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                   loc='upper right', fontsize=24)
        plt.ylim( [1.e-7, 5.e-4] )
    else:
        plt.legend([rho2_line],
                   [r'$\rho_2(\theta)$'],
                   loc='upper right')
        plt.ylim( [1.e-7, 5.e-4] )
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\rho(\theta)$')
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def main():

    args = parse_args()

    # Make the work directory if it does not exist yet.
    cat_dir = os.path.expanduser(args.psf_cats)
    print 'psfcats dir = ',cat_dir

    print os.path.join(cat_dir, '*_psf.fits')
    cat_files = glob.glob(os.path.join(cat_dir, '*_psf.fits'))

    min_sep = 0.5
    max_sep = 20
    bin_size = 0.5
    bin_slop = 0.1

    rho1 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=2)
    rho2 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=2)
    rho3 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=2)
    rho4 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=2)
    rho5 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=2)

    rho1.read(os.path.join(cat_dir,'rho1.fits'))
    rho2.read(os.path.join(cat_dir,'rho2.fits'))
    rho3.read(os.path.join(cat_dir,'rho3.fits'))
    rho4.read(os.path.join(cat_dir,'rho4.fits'))
    rho5.read(os.path.join(cat_dir,'rho5.fits'))

    print 'rho1.xip = ',rho1.xip
    print 'rho1.xim = ',rho1.xim
    print 'npairs = ',rho1.npairs
    print 'rho2 = ',rho2.xip

    plt.clf()
    pretty_rho1(rho1.meanr,rho1.xip,numpy.sqrt(rho1.varxi),rho3.xip,rho4.xip)
    plt.savefig('rho1.pdf')

    plt.clf()
    pretty_rho2(rho2.meanr,rho2.xip,numpy.sqrt(rho2.varxi),rho5.xip)
    plt.savefig('rho2.pdf')


if __name__ == "__main__":
    main()
