#! /usr/bin/env python
# Program to computer rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

import os
import traceback
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Output directory
    parser.add_argument('--output', default='./',
                        help='location of outputs')
    # Exposure inputs
    parser.add_argument('--exp_match', default='',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')

    # options
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')

    args = parser.parse_args()
    return args


def parse_file_name(file_name):
    """Parse the PSFEx file name to get the root name and the chip number
    """
    base_file = os.path.split(file_name)[1]
    print 'base_file = ',base_file

    suff = '_psfcat.psf'
    print 'suff = ',suff
    if not base_file.endswith(suff):
        raise ValueError("Invalid file name "+file)
    l = len(suff)
    print 'l = ',l
    root = base_file[:-len(suff)]
    print 'root = ',root

    ccdnum = int(root.split('_')[-1])
    print 'ccdnum = ',ccdnum
    return root, ccdnum


def read_used(exp_dir, root):
    """Read in the .used.fits file that PSFEx generates with the list of stars that actually
    got used in making the PSFEx file.
    """
    file_name = os.path.join(exp_dir, root + '_psfcat.used.fits')
    import astropy.io.fits as pyfits
    import copy
    f = pyfits.open(file_name, memmap=False)
    data = copy.copy(f[2].data)
    f.close()
    # This has the following columns:
    # SOURCE_NUMBER: 1..n
    # EXTENSION_NUMBER: Seems to be all 1's.
    # CATALOG_NUMBER: Also all 1's.
    # VECTOR_CONTEXT: nx2 array.  Seems to be the same as x,y
    # X_IMAGE: x values  (0-2048)
    # Y_IMAGE: y values (0-4096)
    # DELTAX_IMAGE: All numbers with abs val < 1.  Probably centroid shifts.
    # DELTAY_IMAGE: All numbers with abs val < 1.  Probably centroid shifts.
    # NORM_PSF: Big numbers.  I'm guessing total flux?  Possibly weighted flux.
    # CHI2_PSF: Mostly a little more than 1.  Reduced chi^2 presumably.
    # RESI_PSF: Very small numbers. <1e-4 typically.  Presumably the rms residual.
    return data

def read_findstars(exp_dir, root):
    """Read in the findstars output file.
    """
    file_name = os.path.join(exp_dir, root + '_findstars.fits')
    import astropy.io.fits as pyfits
    import copy
    f = pyfits.open(file_name, memmap=False)
    data = copy.copy(f[1].data)
    f.close()
    # This has the following columns:
    # id: The original id from the SExtractor catalog
    # x: The x position
    # y: The y position
    # sky: The local sky value
    # noise: The estimated noise.  But these are all 0, so I think this isn't being calculated.
    # size_flags: Error flags that occurred when estimating the size
    # mag: The magnitude from SExtractor
    # sg: SExtractor's star/galaxy estimate.  Currently class_star.
    # sigma0: The shapelet sigma that results in a b_11 = 0 shapelet parameter.
    # star_flag: 1 if findstars thought this was a star, 0 otherwise.
    return data
 
def find_fs_index(used_data, fs_data):
    """Find the index in the fs_data records corresponding to each star in used_data.
    """
    import numpy
    index = numpy.zeros(len(used_data),dtype=int)
    used_x = used_data['X_IMAGE']
    used_y = used_data['Y_IMAGE']
    fs_x = fs_data['x']
    fs_y = fs_data['y']
    for i in range(len(used_data)):
        x = used_x[i]
        y = used_y[i]
        close = numpy.where( (fs_x > x-5.) &
                             (fs_x < x+5.) &
                             (fs_y > y-5.) &
                             (fs_y < y+5.) )[0]
        if len(close) == 0:
            print 'Could not find object near x,y = (%f,%f)'%(x,y)
            index[i] = -1
        elif len(close) == 1:
            index[i] = close[0]
        else:
            print 'Multiple objects found near x,y = (%f,%f)'%(x,y)
            amin = numpy.argmin((fs_x[close] - x)**2 + (fs_y[close] - y)**2)
            print 'Found minimum at ',close[amin],', where (x,y) = (%f,%f)'%(
                    fs_x[close[amin]], fs_y[close[amin]])
            index[i] = close[amin]
    return index


def main():
    import glob
    args = parse_args()

    # Make the output directory if it does not exist yet.
    try:
        os.mkdir(args.output)
    except OSError:
        pass

    datadir = '/astro/u/astrodat/data/DES'

    for run,exp in zip(args.runs,args.exps):

        print 'Start work on run, exp = ',run,exp

        # The "output" dir has the PSFEx files.
        exp_dir = os.path.join(args.output,exp)
        print 'exp_dir = ',exp_dir

        # Get the file names in that directory.
        files = glob.glob('%s/%s'%(exp_dir,'*.psf'))
        print 'files = ',files

        for file_name in files:
            print '\nProcessing ', file_name
            flag = 0

            try:
                root, ccdnum = parse_file_name(file_name)
            except:
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue
            print '   root, ccdnum = ',root,ccdnum

            fs_data = read_findstars(exp_dir, root)
            print '   n_tot = ',len(fs_data)

            print '   n_fs = ',fs_data['star_flag'].sum()
            mask = fs_data['star_flag'] == 1

            used_data = read_used(exp_dir, root)
            print '   n_used = ',len(used_data)
            if len(used_data) == 0:
                print '   No stars were used.'
                continue

            index = find_fs_index(used_data, fs_data)
            print '   index = ',index

            fs_xmin = fs_data['x'][mask].min()
            fs_xmax = fs_data['x'][mask].max()
            fs_ymin = fs_data['y'][mask].min()
            fs_ymax = fs_data['y'][mask].max()
            print '   bounds from findstars = ',fs_xmin,fs_xmax,fs_ymin,fs_ymax
            fs_area = (fs_xmax-fs_xmin)*(fs_ymax-fs_ymin)
            print '   area = ',fs_area

            used_xmin = used_data['X_IMAGE'].min()
            used_xmax = used_data['X_IMAGE'].max()
            used_ymin = used_data['Y_IMAGE'].min()
            used_ymax = used_data['Y_IMAGE'].max()
            print '   final bounds of used stars = ',used_xmin,used_xmax,used_ymin,used_ymax
            used_area = (used_xmax-used_xmin)*(used_ymax-used_ymin)
            print '   area = ',used_area
            print '   fraction used = ',float(used_area) / fs_area

            if args.single_ccd:
                break

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
