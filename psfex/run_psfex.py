#! /usr/bin/env python
# Program to run single file in psfex

import os
import traceback
 
# Define the parameters for the blacklist

# How many stars are too few or too many?
FEW_STARS = 50
MANY_STARS = 500
# How high is a high FWHM?  1.8 arcsec / 0.26 arcsec/pixel = 6.9 pixels
HIGH_FWHM = 6.9

# flag values
NO_STARS_FLAG = 1
TOO_FEW_STARS_FLAG = 2
TOO_MANY_STARS_FLAG = 4
TOO_HIGH_FWHM_FLAG = 8
ERROR_FLAG = 16

class NoStarsException(Exception):
    pass


def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Directory arguments
    parser.add_argument('--sex_dir', default='/astro/u/rarmst/soft/bin/',
                        help='location of sextrator executable')
    parser.add_argument('--psfex_dir', default='/astro/u/rarmst/soft/bin/',
                        help='location of psfex executable')
    parser.add_argument('--findstars_dir', default='/astro/u/mjarvis/bin',
                        help='location wl executables')
    parser.add_argument('--output', default='./',
                        help='location of outputs')
    parser.add_argument('--clear_output', default=False, action='store_const', const=True,
                        help='should the output directory be cleared before writing new files?')

    # Exposure inputs
    parser.add_argument('--exp_match', default='',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')

    # Configuration files
    parser.add_argument('--sex_config',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/default.sex',
                        help='sextractor config file')
    parser.add_argument('--psfex_config',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/new.psfex',
                        help='psfex config file')
    parser.add_argument('--findstars_config',
                        default='wl.config +wl_firstcut.config',
                        help='wl config file')
    parser.add_argument('--sex_params',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.param_psfex',
                        help='sextractor param file')
    parser.add_argument('--sex_filter',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.conv',
                        help='name of sextractor filter file')
    parser.add_argument('--sex_nnw',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.nnw',
                        help='name of sextractor star file')
    parser.add_argument('--tapebump_file',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/mask_ccdnum.txt',
                        help='name of tape bump file')
    parser.add_argument('--make_symlinks', default=1, type=int,
                        help='make symlinks from $DESDATA/EXTRA/red/$run/psfex-rerun/$exp')

    # options
    parser.add_argument('--rm_files', default=1, type=int,
                        help='remove unpacked files after finished')
    parser.add_argument('--run_psfex', default=1, type=int,
                        help='run psfex on files')
    parser.add_argument('--use_findstars', default=1, type=int,
                        help='use findstars results in psfex')
    parser.add_argument('--mag_cut', default=-1, type=float,
                        help='remove the top mags using mag_auto')
    parser.add_argument('--nbright_stars', default=10, type=int,
                        help='use median of this many brightest stars for min mag')
    parser.add_argument('--use_tapebumps', default=1, type=int,
                        help='avoid stars in or near tape bumps')
    parser.add_argument('--tapebump_extra', default=2, type=float,
                        help='How much extra room around tape bumps to exclude stars in units of FWHM')
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')


    args = parser.parse_args()
    return args


def read_tapebump_file(file_name):
    """Read and parse the tapebump file if we are going to need it.
    """
    import numpy
    raw_tbdata = numpy.genfromtxt(file_name, delimiter=',')
    # repackage this as  dict (key = ccdnum) of lists of tuples (ymin, xmin, ymax, xmax)
    tbdata = {}
    for d in raw_tbdata:
        ccdnum = int(d[0])
        if ccdnum not in tbdata:
            tbdata[ccdnum] = []
        tbdata[ccdnum].append( (int(d[1]), int(d[2]), int(d[3]), int(d[4])) )
    print 'read in tape bump file.  %d bumps for %d chips'%(len(raw_tbdata),len(tbdata))
    return tbdata

def exclude_tapebumps(tbd, data, extra):
    """
    Process the tapebump data for a particular chip.

    tbd = tbdata for a particular ccdnum,
    data = the input data for the stars
    extra = how much extra distance around the tape bumps to exclude stars in pixels
    """
    import numpy
    # I'm sure it doesn't matter, but add an extra 0.5 pixel to the slop because the tape bump
    # values are in terms of which pixels to include as part of the tape bump.  So the edges
    # are an extra 0.5 pixel outside of that.
    extra += 0.5
    x = data['X_IMAGE']
    y = data['Y_IMAGE']
    masks = [(y>tb[0]-extra) & (x>tb[1]-extra) & (y<tb[2]+extra) & (x<tb[3]+extra) for tb in tbd]
    mask = numpy.any(masks, axis=0)
    if sum(mask) > 0:
        print '   masking %d stars for being in or near a bump'%sum(mask)
        print '       tapebumps = ',tbd
        print '       excluded x,y = ',zip(x[mask],y[mask])
    return data[numpy.logical_not(mask)]


def log_blacklist(blacklist_file, run, exp, ccdnum, flag):
    try:
        with open(blacklist_file,'a') as f:
            f.write("%s %s %d %d\n"%(run,exp,ccdnum,flag))
    except OSError:
        print 'Error opening blacklist.  Wait and try again.'
        import time
        time.sleep(1)
        return log_blacklist(blacklist_file,run,exp,ccdnum,flag)


def parse_file_name(file_name):
    """Parse the file name to get the directory, the root name, and the chip number
    """
    # find out if the file is fpacked by the extension
    base_file = os.path.split(file_name)[1]

    # find the base filename
    if os.path.splitext(base_file)[1] == '.fz':
        base_file=os.path.splitext(base_file)[0]
    
    if os.path.splitext(base_file)[1] != '.fits':
        raise ValueError("Invalid file name "+file)
    root = os.path.splitext(base_file)[0]

    # Get the ccdnum.  We could do this by reading the file and looking for the 
    # CCDNUM keyword.  But it's simpler to just use the filename.  We should be
    # able to assume that the file names are not mangled.
    ccdnum = int(root.split('_')[-1])
    return root, ccdnum


def unpack_file(file_name, odir, logfile):
    """Create the unpacked file in the output directory if necessary.

    If the unpacked file already exists, then a link is made.
    Otherwise funpack is run, outputting the result into the output directory.
    """
    import os
    # find out if the file is fpacked by the extension
    base_file = os.path.split(file_name)[1]

    # find the base filename
    if os.path.splitext(base_file)[1] == '.fz':

        img_file = os.path.join(odir,os.path.splitext(base_file)[0])

        # If an unpacked file does not exist in the output directory 
        if not os.path.exists(img_file):
            print '   unpacking fz file'
            cmd = 'funpack -O {outf} {inf} >> {log} 2>&1'.format(
                outf=img_file, inf=file_name, log=logfile)
            os.system(cmd)
    
    else:
        # If the file is not fpacked, make a symlink into the output directory
        img_file = os.path.join(odir,base_file)
        if not os.path.exists(img_file):
            os.symlink(file_name,odir)

    return img_file

def read_image_header(img_file):
    """Read some information from the image header.

    Currently this is just the SATURATE and FWHM values.

    Returns sat, fwhm
    """
    import pyfits

    hdu = 0

    with pyfits.open(img_file) as pyf:
        sat = -1
        fwhm = 4.
        try:
            sat = pyf[hdu].header['SATURATE']
            fwhm = pyf[hdu].header['FWHM']
        except:
            raise RuntimeError("Cannot read header information from " + img_file)
    return sat, fwhm
 

def run_sextractor(odir, root, img_file, sat, fwhm, logfile,
                   sex_dir, sex_config, sex_params, sex_filter, sex_nnw):
    """Run sextractor, but only if the output file does not exist yet.
    """
    cat_file = os.path.join(odir,root+'_psfcat.fits')

    if not os.path.exists(cat_file):
        print '   running sextractor'

        cat_cmd = "{sex_dir}/sex {img_file}[0] -c {config} -CATALOG_NAME {cat_file} -CATALOG_TYPE FITS_LDAC -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {img_file}[2] -PARAMETERS_NAME {params} -FILTER_NAME {filter}  -STARNNW_NAME {nnw} -DETECT_MINAREA 3 -SEEING_FWHM {fwhm} -SATUR_LEVEL {sat} >> {log} 2>&1".format(
            sex_dir=sex_dir, img_file=img_file, config=sex_config,
            cat_file=cat_file, params=sex_params, filter=sex_filter,
            nnw=sex_nnw, fwhm=fwhm, log=logfile, sat=sat)
        os.system(cat_cmd)
    return cat_file

def run_findstars(odir, root, cat_file, logfile, fs_dir, fs_config):
    """Run findstars, and return a new updated catalog file to use.
    """
    import pyfits
    import numpy
    import copy
    star_file = odir+'/'+root+'_findstars.fits'    

    # run find stars
    print '   running findstars'
    findstars_cmd = '{fs_dir}/findstars {fs_config} root={root} cat_ext=_psfcat.fits stars_file={star_file} input_prefix={odir}/ >> {log} 2>&1'.format(
            fs_dir=fs_dir, fs_config=fs_config, root=root, star_file=star_file, 
            odir=odir, log=logfile)
    #print findstars_cmd
    os.system(findstars_cmd)

    # Make a mask based on which objects findstars decided are stars.
    with pyfits.open(star_file) as pyf:
        mask = pyf[1].data['star_flag']==1
    nstars = numpy.count_nonzero(mask)
    print '   found %d stars'%nstars
    if nstars == 0:
        # Can't really do any of the rest of this, so skip out to the end.
        raise NoStarsException()

    # Read the information from the initial catalog file, including the bogus first two hdus.
    with pyfits.open(cat_file) as pyf:
        # Need to make copy of these to not fail
        hdu1 = copy.copy(pyf[0])
        hdu2 = copy.copy(pyf[1])
        data = pyf[2].data[mask]

        # create new catalog file with only these entries
        hdu3 = pyfits.BinTableHDU(data)
        hdu3.name = 'LDAC_OBJECTS'
        list = pyfits.HDUList([hdu1, hdu2, hdu3])
        new_cat_file = cat_file.replace('psfcat','psfcat_findstars')
        list.writeto(new_cat_file,clobber=True)

    return new_cat_file, nstars

def remove_bad_stars(odir, root, ccdnum, cat_file, tbdata,
                     mag_cut, nbright_stars,
                     use_tapebumps, tapebump_extra, fwhm):
    """Remove stars that are considered bad for some reason.

    Currently these reasons include:
    - Magnitude indicates that the star is significantly contaminated by the brighter/fatter
      effect.
    - Star falls in or near the tape bumps.
    """
    import pyfits
    import numpy
    import copy

    # get the brightest 10 stars that have flags=0 and take the median just in case some
    # strange magnitudes were selected
    with pyfits.open(cat_file) as pyf:
        data = copy.copy(pyf[2].data)

    # Start with a basic FLAGS==0 mask:
    flags_mask = data['FLAGS']==0
    print '   nstars with FLAGS==0 = ',numpy.count_nonzero(flags_mask)
    data = data[flags_mask]

    # Start with the current name.  We will update it below.
    new_cat_file = cat_file

    if mag_cut > 0:
        mags = numpy.sort(data['MAG_AUTO'])
        min_star = numpy.median(mags[0:nbright_stars])
        print '   min mag = ',mags[0]
        print '   median of brightest %d is '%nbright_stars, min_star

        mag_mask = data['MAG_AUTO'] > (min_star+mag_cut)
        print '   select stars dimmer than ',min_star+mag_cut
        print '   which includes %d stars'%numpy.count_nonzero(mag_mask)
            
        data = data[mag_mask]
        print '   after exclude bright: len(data) = ',len(data)
        new_cat_file = new_cat_file.replace('psfcat','psfcat_magcut_%0.1f'%mag_cut)

    if use_tapebumps:
        data = exclude_tapebumps(tbdata[ccdnum], data, tapebump_extra * fwhm)
        print '   after exclude tapebumps: len(data) = ',len(data)
        new_cat_file = new_cat_file.replace('psfcat','psfcat_tb')

    # create new catalog file with only these entries
    with pyfits.open(cat_file) as pyf:
        hdu1 = copy.copy(pyf[0])
        hdu2 = copy.copy(pyf[1])
        hdu3 = pyfits.BinTableHDU(data)
        hdu3.name = 'LDAC_OBJECTS'
        list = pyfits.HDUList([hdu1, hdu2, hdu3])
        # Apparently pyf still needs to be open when this command occurs in order to 
        # be able to handle the hdu1 and hdu2 objects correctly.
        # Hence, we don't read those earlier when we read data.
        list.writeto(new_cat_file,clobber=True)

    return new_cat_file, len(data)


def get_fwhm(cat_file):
    """Get the fwhm from the SExtractor FLUX_RADIUS estimates.
    """
    import pyfits
    import numpy

    # get the brightest 10 stars that have flags=0 and take the median just in case some
    # strange magnitudes were selected
    with pyfits.open(cat_file) as pyf:
        data = pyf[2].data
        flux_radius = data['FLUX_RADIUS']
    return numpy.median(flux_radius)


def run_psfex(odir, root, cat_file, psf_file, used_file, logfile, psfex_dir, psfex_config):
    """Run PSFEx
    """
    print '   running psfex'
    psf_cmd = '{psfex_dir}/psfex {cat_file} -c {config} -OUTCAT_TYPE FITS_LDAC -OUTCAT_NAME {used_file} >> {log} 2>&1'.format(
            psfex_dir=psfex_dir, cat_file=cat_file, config=psfex_config, used_file=used_file,
            log=logfile)
    os.system(psf_cmd)

    # PSFEx generates its output filename from the input catalog name.  If this doesn't match
    # our target name, then rename it.
    actual_psf_file = cat_file.replace('.fits','.psf')
    if psf_file != actual_psf_file:
        os.rename(actual_psf_file, psf_file)

def remove_temp_files(odir, root, *args):
    """Remove odir/root* except for any files listed in the args
    """
    import glob
    files = glob.glob('%s/%s*'%(odir,root))
    for save in args:
        if save in files:
            files.remove(save)

    print '   Removing the following files from ',odir
    for f in files:
        print '       ',os.path.split(f)[1]
        os.remove(f)
    print '   Done'


def make_symlinks(odir, link_dir, *args):
    """Make symlinks for files given in args from link_dir to odir
    """
    print '   Make symlinks from link_dir', link_dir
    if not os.path.isdir(link_dir):
        try: 
            os.makedirs(link_dir, mode=0775)
            print '   made directory ',link_dir
        except OSError:
            if not os.path.isdir(link_dir): raise
    for file in args:
        link = os.path.join(link_dir,os.path.basename(file))
        print '   make link: ',link,' to ',file
        try:
            os.remove(link)
        except OSError:
            pass
        os.symlink(file,link)

def main():
    import glob
    args = parse_args()
    if args.use_tapebumps:
        tbdata = read_tapebump_file(args.tapebump_file)
    blacklist_file = os.path.join(args.output,'psfex_blacklist')

    # Make the output directory if it does not exist yet.
    try:
        os.mkdir(args.output)
    except OSError:
        pass

    datadir = '/astro/u/astrodat/data/DES'

    for run,exp in zip(args.runs,args.exps):

        print 'Start work on run, exp = ',run,exp

        # Make the output directory for this exposure and clear it if necessary.
        odir = os.path.join(args.output,exp)
        if args.clear_output:
            import shutil
            try:
                shutil.rmtree(odir)
            except OSError:
                pass
        try:
            os.makedirs(odir)
        except:
            if not os.path.exists(odir): raise

        # Outputs of shell commands will be logged here.
        logfile = os.path.join(odir, 'log.'+exp)

        # The input directory from the main DESDM reduction location.
        input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))

        # Get the file names in that directory.
        files = glob.glob('%s/%s'%(input_dir,args.exp_match))

        for file_name in files:
            print '\nProcessing ', file_name
            flag = 0

            try:
                root, ccdnum = parse_file_name(file_name)
            except:
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue
            print '   root, ccdnum = ',root,ccdnum

            try:

                if args.run_psfex or args.use_findstars or args.mag_cut>0 or args.use_tapebumps:
                    # Unpack the image file if necessary
                    img_file = unpack_file(file_name, odir, logfile)

                    # extract the saturation level, this is how desdm runs sextractor
                    # we need the fwhm for class star
                    # Also need the fwhm for doing the tape bumps.
                    sat, fwhm = read_image_header(img_file)
                    print '   fwhm = ',fwhm

                    cat_file = run_sextractor(odir, root, img_file, sat, fwhm, logfile,
                                              args.sex_dir, args.sex_config, args.sex_params, 
                                              args.sex_filter, args.sex_nnw)

                # if we want to use only the stars selected by findstars
                if args.use_findstars:
                    tmp = cat_file
                    cat_file, nstars = run_findstars(odir, root, cat_file, logfile,
                                                     args.findstars_dir, args.findstars_config)
                    # Check if there are few or many staras.
                    if nstars < FEW_STARS:
                        print '     -- flag for too few stars: ',nstars
                        flag |= TOO_FEW_STARS_FLAG
                    if nstars > MANY_STARS:
                        print '     -- flag for too many stars: ',nstars
                        flag |= TOO_MANY_STARS_FLAG


                # If we want to cut the brighest magnitudes
                if args.mag_cut>0 or args.use_tapebumps:
                    cat_file, nstars = remove_bad_stars(
                            odir, root, ccdnum, cat_file, tbdata,
                            args.mag_cut, args.nbright_stars,
                            args.use_tapebumps, args.tapebump_extra, fwhm)
                    # Recheck this.
                    if nstars < FEW_STARS:
                        print '     -- flag for too few stars: ',nstars
                        flag |= TOO_FEW_STARS_FLAG
                    if nstars == 0:
                        raise NoStarsException()

                if args.run_psfex or args.use_findstars or args.mag_cut>0 or args.use_tapebumps:
                    # Get the median fwhm of the given stars
                    star_fwhm = get_fwhm(cat_file)
                    print '   fwhm of stars = ',star_fwhm
                    if star_fwhm > HIGH_FWHM:
                        print '     -- flag for too high fwhm'
                        flag |= TOO_HIGH_FWHM_FLAG
                    if star_fwhm > 1.5 * fwhm:
                        print '     -- flag for too high fwhm compared to fwhm from fits header'
                        flag |= TOO_HIGH_FWHM_FLAG
    
                # Need these names even if not running psfex, since we may be updating symlinks.
                psf_file = os.path.join(odir,root+'_psfcat.psf')
                used_file = os.path.join(odir,root+'_psfcat.used.fits')

                if args.run_psfex:
                    run_psfex(odir, root, cat_file, psf_file, used_file, logfile,
                              args.psfex_dir, args.psfex_config)

                if args.rm_files:
                    remove_temp_files(odir, root, psf_file, used_file)

                if args.make_symlinks:
                    link_dir = os.path.join(datadir,'EXTRA/red/%s/psfex-rerun/%s/'%(run,exp))
                    make_symlinks(odir, link_dir, psf_file, used_file)

            except NoStarsException:
                print 'No stars.  Log this in the blacklist and continue.'
                flag |= NO_STARS_FLAG
            except Exception as e:
                print 'Caught exception: ',e
                traceback.print_exc()
                print 'Log this in the blacklist and continue.'
                flag |= ERROR_FLAG

            if flag:
                log_blacklist(blacklist_file,run,exp,ccdnum,flag)

            if args.single_ccd:
                break

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
