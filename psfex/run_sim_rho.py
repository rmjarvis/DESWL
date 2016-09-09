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
import treecorr

def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate rho stats for sim')

    # Drectory arguments
    parser.add_argument('--psf_cats', default='./',
                        help='location of psf_cats directory')

    # Options
    parser.add_argument('--single', default=False, action='store_const', const=True,
                        help='Only do a single exposure (used for debugging)')

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    # Make the work directory if it does not exist yet.
    cat_dir = os.path.expanduser(args.psf_cats)
    print 'psfcats dir = ',cat_dir

    print os.path.join(cat_dir, '*_psf.fits')
    cat_files = sorted(glob.glob(os.path.join(cat_dir, '*-000001_psf.fits')))

    min_sep = 0.5
    max_sep = 20
    bin_size = 0.3
    bin_slop = 0.1

    rho1 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=0)
    rho2 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=0)
    rho3 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=0)
    rho4 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=0)
    rho5 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=bin_size, bin_slop=bin_slop, verbose=0)

    all_ra = []
    all_dec = []
    all_e1 = []
    all_e2 = []
    all_s = []

    for file in cat_files:
        print 'Processing ',file

        try:
            with pyfits.open(file) as pyf:
                data = pyf[1].data.copy()
        except:
            print 'Failed to read data from %s.  Skiping.'%file
            continue
        #print 'len(data) = ',len(data)
        #print 'shape = ',data.shape

        mask = data['flag'] == 0
        #print 'len(mask) = ',len(mask)
        #print 'nonzero(mask) = ',numpy.sum(mask)
        if mask.sum() == 0:
            print '   All objects in this file are flagged.'
            continue

        ra = data['ra'][mask]
        dec = data['dec'][mask]
        e1 = data['e1'][mask]
        e2 = data['e2'][mask]
        s = data['size'][mask]

        m_e1 = data['psf_e1'][mask]
        m_e2 = data['psf_e2'][mask]
        m_s = data['psf_size'][mask]

        all_ra.append(ra)
        all_dec.append(dec)
        all_e1.append(e1)
        all_e2.append(e2)
        all_s.append(s)

        de1 = e1-m_e1
        de2 = e2-m_e2
        dt = (s**2-m_s**2)/s**2

        #print 'e1 = ',e1
        #print 'e2 = ',e2
        #print 'de1 = ',de1
        #print 'de2 = ',de2
        #print 'dt = ',dt

        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                                 k=dt, g1=dt*e1, g2=dt*e2)

        ecat.name = 'ecat'
        decat.name = 'decat'
        dtcat.name = 'dtcat'

        for (cat1, cat2, rho) in [ (decat, decat, rho1),
                                   (ecat, decat, rho2),
                                   (dtcat, dtcat, rho3),
                                   (decat, dtcat, rho4),
                                   (ecat, dtcat, rho5) ]:
            #print 'Doing correlation of %s vs %s'%(cat1.name, cat2.name)

            if cat1 is cat2:
                rho.process_auto(cat1)
            else:
                rho.process_cross(cat1, cat2)

    ra = numpy.concatenate(all_ra)
    dec = numpy.concatenate(all_dec)
    e1 = numpy.concatenate(all_e1)
    e2 = numpy.concatenate(all_e2)
    s = numpy.concatenate(all_s)

    allcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
    varg = treecorr.calculateVarG(allcat)

    for rho in [ rho1, rho2, rho3, rho4, rho5 ]:
        rho.finalize(varg, varg)

    print '\nFinished processing all files'
    rho1.write(os.path.join(cat_dir,'rho1.fits'))
    rho2.write(os.path.join(cat_dir,'rho2.fits'))
    rho3.write(os.path.join(cat_dir,'rho3.fits'))
    rho4.write(os.path.join(cat_dir,'rho4.fits'))
    rho5.write(os.path.join(cat_dir,'rho5.fits'))


if __name__ == "__main__":
    main()
