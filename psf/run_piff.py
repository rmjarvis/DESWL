#! /usr/bin/env python
# Run PSFEx for a set of exposures, including making any necessarily input files.
# It also logs errors into a psf blacklist file.

from __future__ import print_function
import os
import sys
import shutil
import logging
import datetime
import traceback
import numpy as np
import copy
import glob
import time
import fitsio
import pixmappy
import pandas
import galsim
import galsim.des
import piff
import ngmix

import matplotlib
matplotlib.use('Agg') # needs to be done before import pyplot
import matplotlib.pyplot as plt
plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')

# Don't skip columns in describe output  (default is 20, which is a bit too small)
pandas.options.display.max_columns = 200

# Define the parameters for the blacklist

# How many stars are too few or too many?
FEW_STARS = 25
MANY_STARS_FRAC = 0.3
# How high is a high FWHM?  3.6 arcsec / 0.26 arcsec/pixel = 13.8 pixels
#HIGH_FWHM = 13.8
HIGH_FWHM = 3.6  # (We switched to measuring this in arcsec)
NSIG_T_OUTLIER = 4   # How many sigma for a chip to be an outlier in <T>.

# flag values for blacklist
NO_STARS_FLAG = 1
TOO_FEW_STARS_FLAG = 2
TOO_MANY_STARS_FLAG = 4
TOO_HIGH_FWHM_FLAG = 8
FINDSTARS_FAILURE = 16
PSF_FAILURE = 32
ERROR_FLAG = 64
LARGE_SIZE_SPREAD = 128
PIFF_BAD_SOLUTION = 256
OUTLIER_SIZE = 512

# flag values for psf catalog

MAX_CENTROID_SHIFT = 1.0

NOT_USED = 1
BAD_MEASUREMENT = 2
CENTROID_SHIFT = 4
OUTLIER = 8
FAILURE = 32
RESERVED = 64
NOT_STAR = 128
BLACK_FLAG_FACTOR = 512 # blacklist flags are this times the original exposure blacklist flag
                        # blacklist flags go up to 64,

rng = galsim.BaseDeviate(1234)

# array to convert ccdnum to detpos
detpos = [None,'S29','S30','S31','S25','S26','S27','S28','S20','S21','S22',
          'S23','S24','S14','S15','S16','S17','S18','S19','S8','S9',
          'S10','S11','S12','S13','S1','S2','S3','S4','S5','S6','S7',
          'N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11',
          'N12','N13','N14','N15','N16','N17','N18','N19','N20','N21',
          'N22','N23','N24','N25','N26','N27','N28','N29','N30','N31']


class CatastrophicFailure(Exception):
    pass

def ps():
    home = os.path.expanduser('~')
    with open(os.path.join(home,'.desar2')) as f:
        s = f.read().strip()
    return s

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description='Run PSFEx on a set of exposures')

    # Directory arguments
    parser.add_argument('--sex_dir', default='/astro/u/mjarvis/bin/',
                        help='location of sextrator executable')
    parser.add_argument('--piff_exe', default='/astro/u/mjarvis/.conda/envs/py2.7/bin/piffify',
                        help='location of piffify executable')
    parser.add_argument('--findstars_dir', default='/astro/u/mjarvis/bin',
                        help='location wl executables')
    parser.add_argument('--work', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of intermediate outputs')
    parser.add_argument('--scratch', default='/data/mjarvis/y3_piff',
                        help='location of intermediate outputs')
    parser.add_argument('--pixmappy_dir', default='/astro/u/mjarvis/work/y3_piff/astro',
                        help='location of pixmappy astrometric solutions')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--base_exposures',
                        default='/astro/u/mjarvis/work/y3_piff/exposures-ccds-Y3A1_COADD.fits',
                        help='The base file with information about the DES exposures')
    parser.add_argument('--file', default='',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Configuration files
    parser.add_argument('--sex_config',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/y3.sex',
                        help='sextractor config file')
    parser.add_argument('--piff_config',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/piff.yaml',
                        help='piff config file')
    parser.add_argument('--findstars_config',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/y3.config',
                        help='findstars config file')
    parser.add_argument('--sex_params',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/sex.param_piff',
                        help='sextractor param file')
    parser.add_argument('--sex_filter',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/sex.conv',
                        help='name of sextractor filter file')
    parser.add_argument('--sex_nnw',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/sex.nnw',
                        help='name of sextractor star file')
    parser.add_argument('--tapebump_file',
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/mask_ccdnum.txt',
                        help='name of tape bump file')
    parser.add_argument('--make_symlinks', default=0, type=int,
                        help='make symlinks in output dir, rather than move files')
    parser.add_argument('--noweight', default=False, action='store_const', const=True,
                        help='do not try to use a weight image.')


    # Options
    parser.add_argument('--clear_output', default=0, type=int,
                        help='should the output directory be cleared before writing new files?')
    parser.add_argument('--rm_files', default=1, type=int,
                        help='remove unpacked files after finished')
    parser.add_argument('--use_existing', default=0, type=int,
                        help='use previously downloaded files if they exist')
    parser.add_argument('--blacklist', default=1, type=int,
                        help='add failed CCDs to the blacklist')
    parser.add_argument('--run_piff', default=1, type=int,
                        help='run piff on files')
    parser.add_argument('--run_sextractor', default=1, type=int,
                        help='run sextractor to remake input catalog')
    parser.add_argument('--run_findstars', default=1, type=int,
                        help='force a run of findstars to get input star catalog')
    parser.add_argument('--mag_cut', default=-1, type=float,
                        help='remove the top mags using mag_auto')
    parser.add_argument('--min_mag', default=-1, type=float,
                        help='remove stars brighter than this mag')
    parser.add_argument('--nbright_stars', default=1, type=int,
                        help='use median of this many brightest stars for min mag')
    parser.add_argument('--max_mag', default=0, type=float,
                        help='only use stars brighter than this mag')
    parser.add_argument('--use_tapebumps', default=1, type=int,
                        help='avoid stars in or near tape bumps')
    parser.add_argument('--tapebump_extra', default=2, type=float,
                        help='How much extra room around tape bumps to exclude stars in units of FWHM')
    parser.add_argument('--single_ccd', default=0, type=int,
                        help='Only do the specified ccd (used for debugging)')
    parser.add_argument('--reserve', default=0, type=float,
                        help='Reserve some fraction of the good stars for testing')
    parser.add_argument('--get_psfex', default=False, action='store_const', const=True,
                        help='Download the PSFEx files along the way')
    parser.add_argument('--plot_fs', default=False, action='store_const', const=True,
                        help='Make a size-magnitude plot of the findstars output')
    parser.add_argument('--use_ngmix', default=False, action='store_const', const=True,
                        help='Use ngmix rather than hsm for the measurements')

    args = parser.parse_args()
    return args


def read_tapebump_file(file_name, logger):
    """Read and parse the tapebump file if we are going to need it.
    """
    raw_tbdata = np.genfromtxt(file_name, delimiter=',')
    # repackage this as  dict (key = ccdnum) of lists of tuples (ymin, xmin, ymax, xmax)
    tbdata = {}
    for d in raw_tbdata:
        ccdnum = int(d[0])
        if ccdnum not in tbdata:
            tbdata[ccdnum] = []
        tbdata[ccdnum].append( (int(d[1]), int(d[2]), int(d[3]), int(d[4])) )
    logger.info('read in tape bump file.  %d bumps for %d chips'%(len(raw_tbdata),len(tbdata)))
    return tbdata

def exclude_tapebumps(df, tbd, extra):
    """
    Process the tapebump data for a particular chip.

    tbd = tbdata for a particular ccdnum,
    data = the input data for the stars
    extra = how much extra distance around the tape bumps to exclude stars in pixels
    """

def log_blacklist(blacklist_file, exp, ccdnum, flag, logger):
    try:
        with open(blacklist_file,'a') as f:
            f.write("%s %d %d\n"%(exp,ccdnum,flag))
    except OSError as e:
        logger.info(e)
        logger.info('Error opening blacklist.  Wait and try again.')
        time.sleep(1)
        return log_blacklist(blacklist_file,exp,ccdnum,flag)


def unpack_file(file_name, logger):
    """Create the unpacked file in the work directory if necessary.

    If the unpacked file already exists, then it is removed.
    Then funpack is run, outputting the result into the work directory.
    """
    logger.info('unpack %s',file_name)

    img_file = os.path.splitext(file_name)[0]
    if os.path.lexists(img_file):
        logger.info('   %s exists already.  Removing.',img_file)
        os.remove(img_file)
    logger.info('   unpacking fz file')
    cmd = 'funpack -O {outf} {inf}'.format(outf=img_file, inf=file_name)
    logger.info(cmd)
    for itry in range(5):
        run_with_timeout(cmd, 120)
        if os.path.lexists(img_file): break
        logger.info('%s was not properly made.  Retrying.',img_file)
        time.sleep(10)

    if not os.path.lexists(img_file):
        logger.info('Unable to create %s.  Skip this file.',img_file)
        return None

    return img_file

def read_image_header(row, img_file, logger):
    """Read some information from the image header and write into the df row.
    """
    hdu = 0

    # Note: The next line usually works, but fitsio doesn't support CONTINUE lines, which DES
    #       image headers sometimes include.
    #h = fitsio.read_header(img_file, hdu)
    # I don't  care about any of the lines the sometimes use CONITNUE (e.g. OBSERVER), so I
    # just remove them and make the header with the rest of the entries.
    f = fitsio.FITS(img_file)
    header_list = f[hdu].read_header_list()
    header_list = [ d for d in header_list if 'CONTINUE' not in d['name'] ]
    h = fitsio.FITSHDR(header_list)
    try:
        date = h['DATE-OBS']
        date, time = date.strip().split('T',1)

        filter = h['FILTER']
        filter = filter.split()[0]

        sat = h['SATURATE']
        fwhm = h['FWHM']

        ccdnum = int(h['CCDNUM'])
        detpos = h['DETPOS'].strip()

        telra = h['TELRA']
        teldec = h['TELDEC']
        telha = h['HA']
        if galsim.__version__ >= '1.5.1':
            telra = galsim.Angle.from_hms(telra) / galsim.degrees
            teldec = galsim.Angle.from_dms(teldec) / galsim.degrees
            telha = galsim.Angle.from_hms(telha) / galsim.degrees
        else:
            telra = galsim.HMS_Angle(telra) / galsim.degrees
            teldec = galsim.DMS_Angle(teldec) / galsim.degrees
            telha = galsim.HMS_Angle(telha) / galsim.degrees

        airmass = float(h.get('AIRMASS',-999))
        sky = float(h.get('SKYBRITE',-999))
        sigsky = float(h.get('SKYSIGMA',-999))

        tiling = int(h.get('TILING',0))
        hex = int(h.get('HEX',0))

    except Exception as e:
        logger.info("Caught %s",e)
        logger.info("Cannot read header information from %s", img_file)
        raise

    row['date'] = date
    row['time'] = time
    row['sat'] = sat
    row['fits_filter'] = filter
    row['fits_fwhm'] = fwhm
    row['fits_ccdnum'] = ccdnum
    row['telra'] = telra
    row['teldec'] = teldec
    row['telha'] = telha
    row['airmass'] = airmass
    row['sky'] = sky
    row['sigsky'] = sigsky
    row['tiling'] = tiling
    row['hex'] = hex


def run_sextractor(wdir, root, img_file, sat, fwhm, noweight,
                   sex_dir, sex_config, sex_params, sex_filter, sex_nnw, logger):
    """Run sextractor, but only if the output file does not exist yet.
    """
    cat_file = os.path.join(wdir,root+'_cat.fits')

    logger.info('   running sextractor')
    cat_cmd = "{sex_dir}/sex {img_file}[0] -c {config} -CATALOG_NAME {cat_file} -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME {params} -FILTER_NAME {filter}  -STARNNW_NAME {nnw} -DETECT_MINAREA 3".format(
            sex_dir=sex_dir, img_file=img_file, config=sex_config,
            cat_file=cat_file, params=sex_params, filter=sex_filter,
            nnw=sex_nnw, fwhm=fwhm, sat=sat)
    if not noweight:
        cat_cmd += " -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {img_file}[2]".format(img_file=img_file)
    if fwhm != 0:
        cat_cmd += " -SEEING_FWHM {fwhm}".format(fwhm=fwhm)
    if sat != -1:
        cat_cmd += " -SATUR_LEVEL {sat}".format(sat=sat)
    logger.info(cat_cmd)
    run_with_timeout(cat_cmd, 120)

    if not os.path.exists(cat_file) or os.path.getsize(cat_file) == 0:
        logger.info('   Error running SExtractor.  No ouput file was written.')
        logger.info('   Try again, in case it was a fluke.')
        run_with_timeout(cat_cmd, 120)
        if os.path.getsize(cat_file) == 0:
            os.remove(cat_file)
            logger.info('   Error running SExtractor (again).')
            return None
    return cat_file

def run_findstars(row, wdir, fs_dir, fs_config, logger):
    """Run findstars
    """
    star_file = row['star_file']

    # run find stars
    logger.info('   running findstars')
    findstars_cmd = '{fs_dir}/findstars {fs_config} root={root} cat_file={cat_file} stars_file={star_file} input_prefix={wdir}/'.format(
            fs_dir=fs_dir, fs_config=fs_config, root=row['root'], cat_file=row['cat_file'],
            star_file=star_file, wdir=wdir)
    logger.info(findstars_cmd)
    run_with_timeout(findstars_cmd, 120)

    if not os.path.exists(star_file) or os.path.getsize(star_file) == 0:
        logger.info('   Error running findstars.  Rerun with verbose=2.')
        debug_file = star_file.replace('.fits','_fs.debug')
        findstars_cmd = findstars_cmd + ' verbose=2 debug_file=' + debug_file
        logger.info(findstars_cmd)
        run_with_timeout(findstars_cmd, 240)
        logger.info('   The debug file is %s',debug_file)
        if not os.path.exists(star_file) or os.path.getsize(star_file) == 0:
            logger.info('   Error running findstars (again).')
            return None
    return star_file

def read_findstars(star_file, cat_file, magzp, logger):
    """Read the findstars output file
    """
    if not os.path.exists(star_file):
        return None

    # Read the output and make a DataFrome with the contents
    data = fitsio.read(star_file)
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)
    # This has the following columns:
    # id: The original id from the SExtractor catalog
    # x: The x position
    # y: The y position
    # sky: The local sky value
    # noise: The estimated noise.  But these are all 0, so I think this isn't being calculated.
    # size_flags: Error flags that occurred when estimating the size
    # mag: The magnitude from SExtractor
    # sg: SExtractor's star/galaxy estimate.  Currently SPREAD_MODEL.  (Actually, currently none)
    # sigma0: The shapelet sigma that results in a b_11 = 0 shapelet parameter.
    # star_flag: 1 if findstars thought this was a star, 0 otherwise.

    ntot = len(df)
    nstars = df['star_flag'].sum()
    logger.info('   found %d stars',nstars)

    #print('mag range = ',np.min(df['mag']), np.max(df['mag']))
    is_star = df['star_flag'] == 1
    #print('star mag range = ',np.min(df['mag'][is_star]), np.max(df['mag'][is_star]))
    #print('zero point = ',magzp)
    df['mag'] += magzp - 25.
    #print('star mag range => ',np.min(df['mag'][is_star]), np.max(df['mag'][is_star]))

    # Add on some extra information from the sextractor catalog
    sdata = fitsio.read(cat_file, 2)
    assert len(data) == len(sdata)
    df['ra'] = sdata['ALPHAWIN_J2000']
    df['dec'] = sdata['DELTAWIN_J2000']

    return df


def remove_bad_stars(df, ccdnum, tbdata,
                     mag_cut, nbright_stars, min_mag, max_mag,
                     use_tapebumps, tapebump_extra, reserve, fwhm, logger):
    """Remove stars that are considered bad for some reason.

    Currently these reasons include:
    - Magnitude indicates that the star is significantly contaminated by the brighter/fatter
      effect.
    - Star falls in or near the tape bumps.
    """
    use = df['star_flag'] == 1
    mags = np.sort(df['mag'][use])
    if use.sum() == 0:
        df['use'] = use
        if reserve:
            df['reserve'] = use  # All False, so just copy use.
        return use.sum(), len(use)

    min_star = mags[0]
    if max_mag < 0:
        max_mag = min_star + abs(max_mag)

    if mag_cut > 0:
        logger.info('   min star mag = %s',min_star)
        if nbright_stars > 1:
            min_star = np.median(mags[0:nbright_stars])
            logger.info('   median of brightest %d is %s',nbright_stars, min_star)
        if min_star + mag_cut > min_mag:
            min_mag = min_star + mag_cut

    if min_mag > 0:
        use = use & (df['mag'] > min_mag)
        logger.info('   select stars dimmer than %s',min_mag)
        logger.info('   which brings star count to %s',use.sum())

    if max_mag > 0:
        use = use & (df['mag'] < max_mag)
        logger.info('   also select stars brighter than %s',max_mag)
        logger.info('   which brings star count to %s',use.sum())

    if use_tapebumps:
        # I'm sure it doesn't matter, but add an extra 0.5 pixel to the slop because the tape bump
        # values are in terms of which pixels to include as part of the tape bump.  So the edges
        # are an extra 0.5 pixel outside of that.
        extra = tapebump_extra * fwhm + 0.5
        x = df['x']
        y = df['y']
        tbd = tbdata[ccdnum]
        mask = [(y>tb[0]-extra) & (x>tb[1]-extra) & (y<tb[2]+extra) & (x<tb[3]+extra) for tb in tbd]
        mask = np.any(mask, axis=0)
        use = use & ~mask
        logger.info('   excluding tape bumps brings star count to %s',use.sum())

    if reserve:
        #print('   reserve ',reserve)
        n = use.sum()
        perm = np.random.permutation(n)
        n1 = int(reserve * n)
        logger.info('   reserving %s',n1)
        #print('   initial ids = ',df[use]['id'].values)
        # There is surely a more efficient way to do this, but I kept getting confused about
        # when numpy makes a copy vs a view.  This step-by-step calculation works.
        r1 = np.zeros((n), dtype=bool)
        r1[:n1] = True
        r2 = np.zeros((n), dtype=bool)
        r2[perm] = r1
        reserve = np.zeros_like(use, dtype=bool)
        reserve[use] = r2
        #print('   nreserve = ',np.sum(reserve))
        df['reserve'] = reserve
        #print('   reserve ids = ',df['id'][df['reserve']].values)
        use = use & ~reserve
        #print('   final ids = ',df['id'][use].values)
        logger.info('   after reserve: nstars = %s',use.sum())

    df['use'] = use
    return use.sum(), len(use)


def get_fwhm(df):
    """Get the fwhm from the SExtractor FLUX_RADIUS estimates.
    """

    # get the brightest 10 stars that have flags=0 and take the median just in case some
    # strange magnitudes were selected
    use = df['use']
    fwhm = 2.35482 * df[use]['sigma0']  # For Gaussian, fwhm = 2.35482 * sigma
    if len(fwhm) >= 1:
        stats = ( np.min(fwhm), np.max(fwhm),
                  np.mean(fwhm), np.median(fwhm) )
    else:
        stats = (-999., -999., -999., -999.)
    return stats


def plot_fs(df, filename, logger):
    fig, ax = plt.subplots(1,1)

    size = 2. * df['sigma0']**2
    mag = df['mag']

    # Make a simpler flag to use below:
    # 0 = good detection
    # 1 = star candidate
    # 2 = final star
    flag = df['star_flag'].copy()
    if 'piff_flag' in df:
        flag[df['piff_flag']==0] = 2
        flag[df['piff_flag']==RESERVED] = 2  # Also count the reserve stars as psf stars
    logger.info('num detections: %s',np.sum(flag == 0))
    logger.info('num candidates: %s',np.sum(flag == 1))
    logger.info('num final psf stars: %s',np.sum(flag == 2))

    logger.info('mag range = %s .. %s',np.min(mag),np.max(mag))
    logger.info('mag range for size < 1.2 = %s .. %s',np.min(mag[size<1.2]),np.max(mag[size<1.2]))

    ax.set_xlim(np.floor(np.min(mag[size<1.2])), np.ceil(np.max(mag[size<1.2])))
    #ax.set_xlim(14, 24.)
    ax.set_ylim(0., 1.2)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel(r'$T = 2\sigma^2$ (arcsec${}^2$)')

    index = np.argsort(mag)
    n = np.arange(len(mag))

    for i in range(20):
        s = index[n%20==i]

        det = ax.scatter(mag[s][flag[s] == 0], size[s][flag[s] == 0],
                         color='black', marker='o', s=2.)
        cand = ax.scatter(mag[s][flag[s] == 1], size[s][flag[s] == 1],
                          color='magenta', marker='o', s=8.)
        psf = ax.scatter(mag[s][flag[s] == 2], size[s][flag[s] == 2],
                         color='green', marker='*', s=40.)

        legend_items = [det, cand, psf]

    ax.legend(legend_items,
              ['Detected Object', 'Candidate Star', 'PSF Star', 'Modest Star'],
              loc='upper left', frameon=True)

    plt.xlim(right= np.max(mag[mag < 30]))
    plt.tight_layout()
    plt.savefig(filename)
    logger.info('Wrote fs plot to %s',filename)
    plt.close(fig)

    np.savetxt(os.path.splitext(filename)[0] + '.dat',
               np.array(zip(size, mag, flag), dtype='f8, f8, i2'), fmt='%r',
               header='size  mag  flag (0=detected, 1=candidate, 2=psf, 4=bad measurement)')


def run_with_timeout(cmd, timeout_sec):
    # cf. https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    import subprocess, shlex
    from threading import Timer

    proc = subprocess.Popen(shlex.split(cmd))
    kill_proc = lambda p: p.kill()
    timer = Timer(timeout_sec, kill_proc, [proc])
    try:
        timer.start()
        proc.communicate()
    finally:
        timer.cancel()

def run_piff(df, img_file, cat_file, psf_file, piff_exe, piff_config,
             pixmappy, exp, ccdnum, logger):
    """Run Piffify

    Returns True if successful, False if there was a catastrophic failure and no output
    file was written.
    """
    if not os.path.exists(piff_config) and '/' not in piff_config:
        piff_config = os.path.join('/astro/u/mjarvis/rmjarvis/DESWL/psf/',piff_config)
    if os.path.lexists(psf_file):
        logger.info('   deleting existing %s',psf_file)
        os.unlink(psf_file)

    logger.info('   making cat file for piff')
    stars_df = df[df['use']]
    logger.debug('stars_df = \n%s',stars_df.describe())
    piff_cat_file = cat_file.replace('stars','use_stars')
    fitsio.write(piff_cat_file, stars_df.to_records(index=False), clobber=True)

    logger.info('   running piff')
    piff_cmd = '{piff_exe} {config} input.image_file_name={image} input.cat_file_name={cat} output.file_name={psf}'.format(
        piff_exe=piff_exe, config=piff_config, image=img_file, cat=piff_cat_file, psf=psf_file)
    piff_cmd += ' input.wcs.file_name={pixmappy} input.wcs.exp={exp} input.wcs.ccdnum={ccdnum}'.format(
        pixmappy=pixmappy, exp=exp, ccdnum=ccdnum)
    logger.info(piff_cmd)
    if True:
        config = piff.read_config(piff_config)
        config['input']['image_file_name'] = img_file
        config['input']['cat_file_name'] = piff_cat_file
        config['output']['file_name'] = psf_file
        config['input']['wcs']['file_name'] = pixmappy
        config['input']['wcs']['exp'] = exp
        config['input']['wcs']['ccdnum'] = ccdnum
        piff.piffify(config, logger)
    else:
        # Old way using piffify executable.
        run_with_timeout(piff_cmd, 300)  # 5 minutes should be way more than plenty!

        if not os.path.exists(psf_file) or os.path.getsize(psf_file) == 0:
            logger.info('   Error running Piff.  No ouput file was written.')
            logger.info('   Try again, in case it was a fluke.')
            piff_cmd += ' verbose=3'  # Add more debugging so we can see what might have gone wrong.
            run_with_timeout(piff_cmd, 600)  # And double the time just in case that was the problem.
            if not os.path.exists(psf_file) or os.path.getsize(psf_file) == 0:
                logger.info('   Error running Piff (again).')
                return False
    return True

def remove_temp_files(sdir, root, keep_files, logger):
    """Remove sdir/root* except for any files listed in the keep_files list
    """
    files = sorted(glob.glob('%s/%s*'%(sdir,root)))
    for save in keep_files:
        if save in files:
            files.remove(save)

    logger.info('   Removing the following files from %s',sdir)
    for f in files:
        logger.info('       %s',os.path.split(f)[1])
        os.remove(f)
    logger.info('   Done')


def wget(url_base, path, wdir, file, logger):
    url = url_base + path + file
    full_file = os.path.join(wdir,file)

    if not os.path.isfile(full_file):
        logger.info('Downloading %s',full_file)
        # Sometimes this fails with an "http protocol error, bad status line".
        # Maybe from too many requests at once or something.  So we retry up to 5 times.
        nattempts = 5
        cmd = 'wget -q --no-check-certificate %s -O %s'%(url, full_file)
        for attempt in range(1,nattempts+1):
            logger.info('%s  (attempt %d)',cmd, attempt)
            run_with_timeout(cmd, 300)
            if os.path.exists(full_file):
                break
    return full_file


def hsm(im, wt, logger):
    #print('im stats: ',im.array.min(),im.array.max(),im.array.mean(),np.median(im.array))
    #print('wt = ',wt)
    #if wt:
        #print('im stats: ',wt.array.min(),wt.array.max(),wt.array.mean(),np.median(wt.array))
    flag = 0
    try:
        shape_data = im.FindAdaptiveMom(weight=wt, strict=False)
        #print('shape_data = ',shape_data)
    except Exception as e:
        logger.info(e)
        logger.info(' *** Bad measurement (caught exception).  Mask this one.')
        flag |= BAD_MEASUREMENT

    if shape_data.moments_status != 0:
        logger.info('status = %s',shape_data.moments_status)
        logger.info(' *** Bad measurement (hsm status).  Mask this one.')
        flag |= BAD_MEASUREMENT

    if galsim.__version__ >= '1.5.1':
        dx = shape_data.moments_centroid.x - im.true_center.x
        dy = shape_data.moments_centroid.y - im.true_center.y
    else:
        dx = shape_data.moments_centroid.x - im.trueCenter().x
        dy = shape_data.moments_centroid.y - im.trueCenter().y
    #print('dx, dy = ',dx,dy)
    if dx**2 + dy**2 > MAX_CENTROID_SHIFT**2:
        logger.info(' *** Centroid shifted by %f,%f in hsm.  Mask this one.',dx,dy)
        flag |= CENTROID_SHIFT

    flux = shape_data.moments_amp
    #print('flux = ',flux)

    # Account for the image wcs
    if im.wcs.isPixelScale():
        g1 = shape_data.observed_shape.g1
        g2 = shape_data.observed_shape.g2
        T = 2 * shape_data.moments_sigma**2 * im.scale**2
        #print('simple shape = ',g1,g2,T)
    else:
        e1 = shape_data.observed_shape.e1
        e2 = shape_data.observed_shape.e2
        s = shape_data.moments_sigma
        #print('simple shape = ',e1,e2,s)

        if galsim.__version__ >= '1.5.1':
            jac = im.wcs.jacobian(im.true_center)
        else:
            jac = im.wcs.jacobian(im.trueCenter())
        M = np.matrix( [[ 1 + e1, e2 ], [ e2, 1 - e1 ]] ) * s*s
        J = jac.getMatrix()
        M = J * M * J.T

        e1 = (M[0,0] - M[1,1]) / (M[0,0] + M[1,1])
        e2 = (2.*M[0,1]) / (M[0,0] + M[1,1])
        T = M[0,0] + M[1,1]

        shear = galsim.Shear(e1=e1, e2=e2)
        g1 = shear.g1
        g2 = shear.g2
        #print('distorted shape = ',g1,g2,T)

    return dx, dy, g1, g2, T, flux, flag

def make_ngmix_prior(T, pixel_scale):
    from ngmix import priors, joint_prior

    # centroid is 1 pixel gaussian in each direction
    cen_prior=priors.CenPrior(0.0, 0.0, pixel_scale, pixel_scale)

    # g is Bernstein & Armstrong prior with sigma = 0.1
    gprior=priors.GPriorBA(0.1)

    # T is log normal with width 0.2
    Tprior=priors.LogNormal(T, 0.2)

    # flux is the only uninformative prior
    Fprior=priors.FlatPrior(-10.0, 1.e10)

    prior=joint_prior.PriorSimpleSep(cen_prior, gprior, Tprior, Fprior)
    return prior


def ngmix_fit(im, wt, fwhm, x, y, logger):
    flag = 0
    dx, dy, g1, g2, flux = 0., 0., 0., 0., 0.
    T_guess = (fwhm / 2.35482)**2 * 2.
    T = T_guess
    #print('fwhm = %s, T_guess = %s'%(fwhm, T_guess))
    try:
        #hsm_dx,hsm_dy,hsm_g1,hsm_g2,hsm_T,hsm_flux,hsm_flag = hsm(im, None, logger)
        #logger.info('hsm: %s, %s, %s, %s, %s, %s, %s',hsm_dx,hsm_dy,hsm_g1,hsm_g2,hsm_T,hsm_flux,hsm_flag)
        #if hsm_flag != 0:
            #print('hsm: ',g1,g2,T,flux,hsm_flag)
            #print('Bad hsm measurement.  Reverting to g=(0,0) and T=T_guess = %s'%(T_guess))
            #T = T_guess
        #elif np.abs(np.log(T/T_guess)) > 0.5:
            #print('hsm: ',g1,g2,T,flux,hsm_flag)
            #print('T = %s is not near T_guess = %s.  Reverting to T_guess'%(T,T_guess))
            #T = T_guess
        if galsim.__version__ >= '1.5.1':
            wcs = im.wcs.local(im.center)
        else:
            wcs = im.wcs.local(im.center())

        prior = make_ngmix_prior(T, wcs.minLinearScale())

        if galsim.__version__ >= '1.5.1':
            cen = im.true_center - im.origin
        else:
            cen = im.trueCenter() - im.origin()
        jac = ngmix.Jacobian(wcs=wcs, x=cen.x + x - int(x+0.5), y=cen.y + y - int(y+0.5))
        if wt is None:
            obs = ngmix.Observation(image=im.array, jacobian=jac)
        else:
            obs = ngmix.Observation(image=im.array, weight=wt.array, jacobian=jac)

        lm_pars = {'maxfev':4000}
        runner=ngmix.bootstrap.PSFRunner(obs, 'gauss', T, lm_pars, prior=prior)
        runner.go(ntry=3)

        ngmix_flag = runner.fitter.get_result()['flags']
        gmix = runner.fitter.get_gmix()
    except Exception as e:
        logger.info(e)
        logger.info(' *** Bad measurement (caught exception).  Mask this one.')
        flag |= BAD_MEASUREMENT
        return dx,dy,g1,g2,T,flux,flag

    if ngmix_flag != 0:
        logger.info(' *** Bad measurement (ngmix flag = %d).  Mask this one.',ngmix_flag)
        flag |= BAD_MEASUREMENT

    dx, dy = gmix.get_cen()
    if dx**2 + dy**2 > MAX_CENTROID_SHIFT**2:
        logger.info(' *** Centroid shifted by %f,%f in ngmix.  Mask this one.',dx,dy)
        flag |= CENTROID_SHIFT

    g1, g2, T = gmix.get_g1g2T()
    if abs(g1) > 0.5 or abs(g2) > 0.5:
        logger.info(' *** Bad shape measurement (%f,%f).  Mask this one.',g1,g2)
        flag |= BAD_MEASUREMENT

    flux = gmix.get_flux() / wcs.pixelArea()  # flux is in ADU.  Should ~ match sum of pixels
    #logger.info('ngmix: %s %s %s %s %s %s %s',dx,dy,g1,g2,T,flux,flag)
    return dx, dy, g1, g2, T, flux, flag

def measure_star_shapes(df, image_file, noweight, wcs, use_ngmix, fwhm, logger):
    """Measure shapes of the raw stellar images at each location.
    """
    logger.info('Read in stars in file: %s',image_file)

    ind = df.index[df['star_flag'] == 1]
    logger.info('ind = %s',ind)
    n_psf = len(ind)
    logger.info('n_psf = %s',n_psf)

    df['obs_dx'] = [ -999. ] * len(df)
    df['obs_dy'] = [ -999. ] * len(df)
    df['obs_e1'] = [ -999. ] * len(df)
    df['obs_e2'] = [ -999. ] * len(df)
    df['obs_T'] = [ -999. ] * len(df)
    df['obs_flux'] = [ -999. ] * len(df)
    df['obs_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'obs_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'obs_flag'] |= RESERVED
        df.loc[~df['use'] & ~df['reserve'], 'obs_flag'] |= NOT_USED
    else:
        df.loc[~df['use'], 'obs_flag'] |= NOT_USED

    full_image = galsim.fits.read(image_file, hdu=0)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        full_weight = galsim.fits.read(image_file, hdu=2)
        full_weight.array[full_weight.array < 0] = 0.

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]

        #print('Measure shape for star at ',x,y)
        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2,
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        if noweight:
            wt = None
        else:
            wt = full_weight[b]

        if use_ngmix:
            dx, dy, e1, e2, T, flux, flag = ngmix_fit(im, wt, fwhm, x, y, logger)
        else:
            dx, dy, e1, e2, T, flux, flag = hsm(im, wt, logger)
        #logger.info('ngmix measurement: (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
        if np.any(np.isnan([dx,dy,e1,e2,T,flux])):
            logger.info(' *** NaN detected (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
            flag |= BAD_MEASUREMENT
        else:
            df.loc[i, 'obs_dx'] = dx
            df.loc[i, 'obs_dy'] = dy
            df.loc[i, 'obs_e1'] = e1
            df.loc[i, 'obs_e2'] = e2
            df.loc[i, 'obs_T'] = T
            df.loc[i, 'obs_flux'] = flux
        df.loc[i, 'obs_flag'] |= flag
    logger.info('final obs_flag = %s',df['obs_flag'][ind].values)
    #print('df[ind] = ',df.loc[ind].describe())
    flag_outliers(df, ind, 'obs', 4., logger)

    # Any stars that weren't measurable here, don't use for PSF fitting.
    df.loc[df['obs_flag']!=0, 'use'] = False

def flag_outliers(df, ind, prefix, nsig, logger):
    logger.info('ind = %s',ind)
    mask = df[prefix + '_flag'][ind] == 0
    ind = ind[mask]
    if len(ind) == 0:
        return
    logger.info('ind => %s',ind)
    mean_e1 = np.mean(df[prefix + '_e1'][ind])
    mean_e2 = np.mean(df[prefix + '_e2'][ind])
    mean_T = np.mean(df[prefix + '_T'][ind])
    std_e1 = np.std(df[prefix + '_e1'][ind])
    std_e2 = np.std(df[prefix + '_e2'][ind])
    std_T = np.std(df[prefix + '_T'][ind])
    logger.info('e1 = %s +- %s',mean_e1,std_e1)
    logger.info('e2 = %s +- %s',mean_e2,std_e2)
    logger.info('T = %s +- %s',mean_T,std_T)
    outlier = np.abs(df[prefix + '_e1'][ind] - mean_e1) > nsig*std_e1
    outlier |= np.abs(df[prefix + '_e2'][ind] - mean_e2) > nsig*std_e2
    outlier |= np.abs(df[prefix + '_T'][ind] - mean_T) > nsig*std_T
    logger.info('n outliers = %s = %s',np.sum(outlier),len(ind[outlier]))
    if len(ind[outlier]) > 0:
        logger.info('outlier = %s',ind[outlier])
        logger.info('outlier x = %s',df['x'][ind[outlier]])
        logger.info('outlier y = %s',df['y'][ind[outlier]])
        logger.info('outlier e1 = %s',df[prefix + '_e1'][ind[outlier]])
        logger.info('outlier e2 = %s',df[prefix + '_e2'][ind[outlier]])
        logger.info('outlier T = %s',df[prefix + '_T'][ind[outlier]])
        df.loc[ind[outlier], prefix + '_flag'] |= OUTLIER
        logger.info('Repeat with the new subset')
        # Increase nsig by 20% so we don't just keep chipping away at the edge of a normal
        # distribution.
        return flag_outliers(df, ind, prefix, nsig*1.2, logger)

def find_index(x1, y1, x2, y2, logger):
    """Find the index of the closest point in (x2,y2) to each (x1,y1)

    Any points that do not have a corresponding point within 1 arcsec gets index = -1.
    """
    index = np.zeros(len(x1),dtype=int)
    for i in range(len(x1)):
        close = np.where( (x2 > x1[i]-1.) &
                          (x2 < x1[i]+1.) &
                          (y2 > y1[i]-1.) &
                          (y2 < y1[i]+1.) )[0]
        if len(close) == 0:
            #print('Could not find object near x,y = (%f,%f)'%(x,y))
            index[i] = -1
        elif len(close) == 1:
            index[i] = close[0]
        else:
            logger.info('Multiple objects found near x,y = (%f,%f)',x1[i],y1[i])
            amin = np.argmin((x2[close] - x1[i])**2 + (y2[close] - y1[i])**2)
            logger.info('Found minimum at %s, where (x,y) = (%f,%f)',
                    close[amin],x2[close[amin]], y2[close[amin]])
            index[i] = close[amin]
    return index


def measure_piff_shapes(df, psf_file, image_file, noweight, wcs, use_ngmix, fwhm, row, logger):
    """Measure shapes of the Piff solution at each location.
    """
    logger.info('Read in Piff file: %s',psf_file)

    ind = df.index[df['star_flag'] == 1]
    logger.info('ind = %s',ind)
    n_psf = len(ind)
    logger.info('n_psf = %s',n_psf)

    df['piff_dx'] = [ -999. ] * len(df)
    df['piff_dy'] = [ -999. ] * len(df)
    df['piff_e1'] = [ -999. ] * len(df)
    df['piff_e2'] = [ -999. ] * len(df)
    df['piff_T'] = [ -999. ] * len(df)
    df['piff_flux'] = [ -999. ] * len(df)
    df['piff_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'piff_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'piff_flag'] |= RESERVED
    df.loc[~df['use'], 'piff_flag'] |= NOT_USED

    # Flag the stars not used by Piff to make the solution.
    used_data = fitsio.read(psf_file, ext='psf_stars')
    logger.info('Piff used %d stars for the solution',len(used_data))
    used_xmin = used_data['x'].min()
    used_xmax = used_data['x'].max()
    used_ymin = used_data['y'].min()
    used_ymax = used_data['y'].max()
    logger.info('   final bounds of used stars = %s .. %s, %s .. %s',
                used_xmin,used_xmax,used_ymin,used_ymax)
    used_area = (used_xmax-used_xmin)*(used_ymax-used_ymin)
    logger.info('   area = %s',used_area)
    used_index = find_index(df['x'], df['y'], used_data['x'], used_data['y'], logger)
    df.loc[used_index < 0, 'piff_flag'] |= NOT_USED

    try:
        psf = piff.read(psf_file)
    except Exception as e:
        logger.info('Caught %s',e)
        df.loc[ind, 'piff_flag'] = FAILURE
        return PSF_FAILURE

    full_image = galsim.fits.read(image_file, hdu=0)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        full_weight = galsim.fits.read(image_file, hdu=2)
        full_weight.array[full_weight.array < 0] = 0.

    stamp_size = 48

    nbad = 0

    good_bounds = galsim.BoundsD()
    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]
        #print('Measure Piff model shape at ',x,y)

        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2,
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        im = psf.draw(x=x, y=y, image=im)
        im *= df['obs_flux'].iloc[i]

        if noweight:
            wt = None
        else:
            wt = full_weight[b]
            var = wt.copy()
            var.invertSelf()
            im.addNoise(galsim.VariableGaussianNoise(rng, var))

        if use_ngmix:
            dx, dy, e1, e2, T, flux, flag = ngmix_fit(im, wt, fwhm, x, y, logger)
        else:
            dx, dy, e1, e2, T, flux, flag = hsm(im, wt, logger)

        #if flag != 0:
            #im1_name = 'failed_image_%d.fits'%i
            #im2_name = 'failed_model_%d.fits'%i
            #print('Writing images to ',im1_name,im2_name)
            #local_wcs = wcs.local(image_pos=galsim.PositionD(x,y))
            #full_image[b].view(wcs=local_wcs).write(im1_name)
            #im.view(wcs=local_wcs).write(im2_name)
        if np.any(np.isnan([dx,dy,e1,e2,T,flux])):
            logger.info(' *** NaN detected (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
            flag |= BAD_MEASUREMENT
            nbad += 1
        else:
            df.loc[i, 'piff_dx'] = dx
            df.loc[i, 'piff_dy'] = dy
            df.loc[i, 'piff_e1'] = e1
            df.loc[i, 'piff_e2'] = e2
            df.loc[i, 'piff_T'] = T
            df.loc[i, 'piff_flux'] = flux
            good_bounds += galsim.PositionD(x,y)

        df.loc[i, 'piff_flag'] |= flag
    logger.info('final piff_flag = %s',df['piff_flag'][ind].values)
    #print('df[ind] = ',df.loc[ind].describe())
    flag_outliers(df, ind, 'piff', 4., logger)

    row['piff_xmin'] = good_bounds.xmin
    row['piff_xmax'] = good_bounds.xmax
    row['piff_ymin'] = good_bounds.ymin
    row['piff_ymax'] = good_bounds.ymax
    row['piff_chisq'] = psf.chisq
    row['piff_dof'] = psf.dof

    logger.info('nbad = %d/%d',nbad,len(ind))
    logger.info('chisq = %f, dof = %f',psf.chisq,psf.dof)
    if nbad > len(ind)/3 or psf.chisq > 1.5*psf.dof:
        return PIFF_BAD_SOLUTION
    else:
        return 0


def measure_psfex_shapes(df, psfex_file, image_file, noweight, wcs, use_ngmix, fwhm, logger):
    """Measure shapes of the PSFEx solution at each location.
    """
    logger.info('Read in PSFEx file: %s',psfex_file)

    ind = df.index[df['star_flag'] == 1]
    logger.info('ind = %s',ind)
    n_psf = len(ind)
    logger.info('n_psf = %s',n_psf)

    df['psfex_dx'] = [ -999. ] * len(df)
    df['psfex_dy'] = [ -999. ] * len(df)
    df['psfex_e1'] = [ -999. ] * len(df)
    df['psfex_e2'] = [ -999. ] * len(df)
    df['psfex_T'] = [ -999. ] * len(df)
    df['psfex_flux'] = [ -999. ] * len(df)
    df['psfex_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'psfex_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'psfex_flag'] |= RESERVED
    df.loc[~df['use'], 'psfex_flag'] |= NOT_USED

    try:
        psf = galsim.des.DES_PSFEx(psfex_file, image_file)
    except Exception as e:
        logger.info('Caught %s',e)
        df.loc[ind, 'psfex_flag'] = FAILURE
        return

    full_image = galsim.fits.read(image_file, hdu=0)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        full_weight = galsim.fits.read(image_file, hdu=2)
        full_weight.array[full_weight.array < 0] = 0.

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]
        #print('Measure PSFEx model shape at ',x,y)
        image_pos = galsim.PositionD(x,y)
        psf_i = psf.getPSF(image_pos)

        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2,
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        im = psf_i.drawImage(image=im, method='no_pixel')
        im *= df['obs_flux'].iloc[i]

        if noweight:
            wt = None
        else:
            wt = full_weight[b]
            var = wt.copy()
            var.invertSelf()
            im.addNoise(galsim.VariableGaussianNoise(rng, var))

        if use_ngmix:
            dx, dy, e1, e2, T, flux, flag = ngmix_fit(im, wt, fwhm, x, y, logger)
        else:
            dx, dy, e1, e2, T, flux, flag = hsm(im, wt, logger)
        if np.any(np.isnan([dx,dy,e1,e2,T,flux])):
            logger.info(' *** NaN detected (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
            flag |= BAD_MEASUREMENT
        else:
            df.loc[i, 'psfex_dx'] = dx
            df.loc[i, 'psfex_dy'] = dy
            df.loc[i, 'psfex_e1'] = e1
            df.loc[i, 'psfex_e2'] = e2
            df.loc[i, 'psfex_T'] = T
            df.loc[i, 'psfex_flux'] = flux
        df.loc[i, 'psfex_flag'] |= flag
    logger.info('final psfex_flag = %s',df['psfex_flag'][ind].values)
    #print('df[ind] = ',df.loc[ind].describe())
    flag_outliers(df, ind, 'psfex', 4., logger)

def run_single_ccd(row, args, wdir, sdir, tbdata, which_zone, logger):

    # The url to use up to just before OPS
    url_base = 'https://rmjarvis:%s@desar2.cosmology.illinois.edu/DESFiles/desarchive/'%ps()

    key, expnum, ccdnum, band = row['key'], row['expnum'], row['ccdnum'], row['band']
    magzp = row['magzp']
    path = row['path'].strip()
    logger.info('path = %s',path)

    # Use a well-defined seed so results are repeatable if we see a problem.
    np.random.seed(((expnum+76876876) * (ccdnum+23424524) * 8675309) % 2**32)

    # Store all information about the stars in a pandas data frame.
    df = None
    flag = 0

    keep_files = []

    try:

        # Download the files we need:
        base_path, _, _, image_file_name = path.rsplit('/',3)
        root, ext = image_file_name.rsplit('_',1)
        logger.info('root, ext = |%s| |%s|',root,ext)

        image_file = os.path.join(sdir, root + '_' + ext)
        if not (args.use_existing and os.path.exists(image_file)):
            image_file = wget(url_base, base_path + '/red/immask/', sdir, root + '_' + ext, logger)
        logger.info('image_file = %s',image_file)
        row['root'] = root
        row['image_file'] = image_file

        bkg_file = os.path.join(sdir, root + '_bkg.fits.fz')
        if not (args.use_existing and os.path.exists(bkg_file)):
            bkg_file = wget(url_base, base_path + '/red/bkg/', sdir, root + '_bkg.fits.fz', logger)
        logger.info('bkg_file = %s',bkg_file)
        row['bkg_file'] = bkg_file

        # Unpack the image file if necessary
        unpack_image_file = unpack_file(image_file, logger)
        if unpack_image_file is None:
            flag |= ERROR_FLAG
            raise CatastrophicFailure()
        logger.info('unpacked to %s',unpack_image_file)

        # Subtract off the background right from the start
        with fitsio.FITS(unpack_image_file, 'rw') as f:
            bkg = fitsio.read(bkg_file)
            #print('after read bkg')
            img = f[0].read()
            img -= bkg
            #print('after subtract bkg')
            f[0].write(img)
            #print('after write new img')
        logger.info('subtracted off background image')

        if args.get_psfex:
            psfex_file = os.path.join(sdir, root + '_psfexcat.psf')
            if not (args.use_existing and os.path.exists(psfex_File)):
                psfex_file = wget(url_base, base_path + '/psf/', wdir, root + '_psfexcat.psf', logger)
            logger.info('psfex_file = %s',psfex_file)
            row['psfex_file'] = psfex_file
            keep_files.append(psfex_file)

        read_image_header(row, unpack_image_file, logger)
        #print('read image header.  row = ',row)
        sat = row['sat']
        fits_fwhm = row['fits_fwhm']

        # Run Sextractor
        cat_file = os.path.join(sdir, root + '_cat.fits')
        # Unfortunately, the desdm catalog in the cat directory doesn't seem to be
        # complete.  Looks like it only has stars maybe?  Which is nominally fine, but
        # we're not sure if we trust their star selection.  So we run sextractor ourself.
        if args.run_sextractor or not os.path.isfile(cat_file):
            # Also need the fwhm for doing the tape bumps.
            cat_file = run_sextractor(sdir, root, unpack_image_file, sat, fits_fwhm,
                    args.noweight, args.sex_dir, args.sex_config, args.sex_params,
                    args.sex_filter, args.sex_nnw, logger)
            if cat_file is None:
                flag |= ERROR_FLAG
                raise CatastrophicFailure()
        logger.info('cat_file = %s',cat_file)
        row['cat_file'] = cat_file

        # Run findstars
        star_file = os.path.join(sdir, root + '_stars.fits')
        row['star_file'] = star_file
        if args.run_findstars or not os.path.isfile(star_file):
            star_file = run_findstars(row, sdir, args.findstars_dir, args.findstars_config, logger)
            if star_file is None:
                flag |= FINDSTARS_FAILURE
                raise CatastrophicFailure()

        # Read the findstars catalog (need to do this even if star_file exists)
        df = read_findstars(star_file, cat_file, magzp, logger)
        if df is None:
            logger.info('     -- flag for findstars failure')
            flag |= FINDSTARS_FAILURE
            raise CatastrophicFailure()

        # Cut the brighest magnitudes or other exclusions/reservations
        nstars, ntot = remove_bad_stars(
                df, ccdnum, tbdata,
                args.mag_cut, args.nbright_stars, args.min_mag, args.max_mag,
                args.use_tapebumps, args.tapebump_extra, args.reserve, fits_fwhm, logger)

        # Check if there are few or many staras.
        if nstars == 0:
            flag |= NO_STARS_FLAG
        if nstars < FEW_STARS:
            logger.info('     -- flag for too few stars: %s',nstars)
            flag |= TOO_FEW_STARS_FLAG
        if nstars > MANY_STARS_FRAC * ntot:
            logger.info('     -- flag for too many stars: %d/%d',nstars,ntot)
            flag |= TOO_MANY_STARS_FLAG
        row['ntot'] = ntot
        row['nstars'] = df['star_flag'].sum()
        row['nstars_piff'] = nstars
        row['nreserve'] = np.sum(df['reserve'])

        # Get the median fwhm of the given stars
        # Returns min, max, mean, median.  We use median, which is index 3.
        star_fwhm = get_fwhm(df)
        logger.info('   fwhm of stars = %s',star_fwhm)
        logger.info('   cf. header fwhm = %s',fits_fwhm * 0.26)
        if star_fwhm[3] > HIGH_FWHM:
            logger.info('     -- flag for too high fwhm')
            flag |= TOO_HIGH_FWHM_FLAG
        if star_fwhm[3] > 1.5 * fits_fwhm * 0.26:
            logger.info('     -- flag for too high fwhm compared to fwhm from fits header')
            flag |= TOO_HIGH_FWHM_FLAG
        row['star_fwhm'] = star_fwhm[3]
        fwhm = star_fwhm[3]

        # Get the pixmappy wcs for this ccd for correcting the shapes
        dp = detpos[ccdnum]
        wz = np.where((which_zone['expnum'] == expnum) & (which_zone['detpos'] == dp))[0][0]
        logger.info('row in which_zone is %s',wz)
        #print('  ',which_zone[wz])
        zone = which_zone['zone'][wz]
        logger.info('zone = %s',zone)
        row['zone'] = zone
        pixmappy_file = os.path.join(args.pixmappy_dir, 'zone%03d.astro'%zone)
        logger.info('pixmappy_file = %s',pixmappy_file)
        wcs = pixmappy.GalSimWCS(pixmappy_file, exp=expnum, ccdnum=ccdnum, default_color=0)
        wcs._color = 0.  # For now.  Revisit when doing color-dependent PSF.

        # Measure the shpes and sizes of the stars
        measure_star_shapes(df, unpack_image_file, args.noweight, wcs,
                            args.use_ngmix, fwhm, logger)

        # Another check is that the spread in star sizes isn't too large
        obs_T = df.loc[df['obs_flag']==0, 'obs_T']
        if np.std(obs_T) > 0.2 * np.mean(obs_T):
            logger.info('     -- flag for large size spread: %f > 0.15 * %f',
                        np.std(obs_T), np.mean(obs_T))
            flag |= LARGE_SIZE_SPREAD

        # Run piff
        psf_file = os.path.join(wdir, root + '_piff.fits')
        row['piff_file'] = psf_file
        keep_files.append(psf_file)
        if args.run_piff and (flag & NO_STARS_FLAG) == 0:
            success = run_piff(df, unpack_image_file, star_file, psf_file,
                               args.piff_exe, args.piff_config,
                               pixmappy_file, expnum, ccdnum, logger)
            if not success:
                flag |= PSF_FAILURE

        # Measure the psf shapes
        if args.run_piff and (flag & (NO_STARS_FLAG | PSF_FAILURE) == 0):
            piff_flag = measure_piff_shapes(df, psf_file, unpack_image_file, args.noweight, wcs,
                                            args.use_ngmix, fwhm, row, logger)
            flag |= piff_flag
            good = (df['piff_flag'] == 0) & (df['obs_flag'] == 0)
            ngood = np.sum(good)
            row['ngood'] = ngood
            logger.info('Number of detected objects = %s',row['ntot'])
            logger.info('Number of stars = %s',row['nstars'])
            logger.info('Number of stars given to Piff = %s',row['nstars_piff'])
            logger.info('Number of good stars = %s',row['ngood'])
            logger.info('Number of reserve stars = %s',row['nreserve'])
            if ngood < FEW_STARS:
                logger.info('     -- flag for too few good stars: %s',ngood)
                flag |= TOO_FEW_STARS_FLAG

            de1 = df.loc[good, 'piff_e1'] - df.loc[good, 'obs_e1']
            de2 = df.loc[good, 'piff_e2'] - df.loc[good, 'obs_e2']
            dT = df.loc[good, 'piff_T'] - df.loc[good, 'obs_T']
            logger.info('de1 = %s, %s, %s, %s',np.min(de1),np.max(de1),np.mean(de1),np.std(de1))
            logger.info('de2 = %s, %s, %s, %s',np.min(de2),np.max(de2),np.mean(de2),np.std(de2))
            logger.info('dT = %s, %s, %s, %s',np.min(dT),np.max(dT),np.mean(dT),np.std(dT))
        else:
            de1 = -999.
            de2 = -999.
            dT = -999.

        # After Piff rejects outliers, the spread should be even smaller.  Don't go too small
        # thought, since some chips do have a mild spread in sizes across CCD.
        if args.run_piff and flag == 0:
            obs_T = df.loc[df['piff_flag']==0, 'obs_T']
            if np.std(obs_T) > 0.15 * np.mean(obs_T):
                logger.info('     -- flag for large size spread: %f > 0.15 * %f',
                            np.std(obs_T), np.mean(obs_T))
                flag |= LARGE_SIZE_SPREAD

        if args.plot_fs:
            fs_plot_file = os.path.join(wdir, root + '_fs.pdf')
            plot_fs(df, fs_plot_file, logger)
            keep_files.append(fs_plot_file)

        if args.get_psfex:
            measure_psfex_shapes(df, psfex_file, unpack_image_file, args.noweight, wcs,
                                 args.use_ngmix, fwhm, logger)
            xgood = (df['psfex_flag'] == 0) & (df['obs_flag'] == 0)
            xde1 = df.loc[xgood, 'psfex_e1'] - df.loc[xgood, 'obs_e1']
            xde2 = df.loc[xgood, 'psfex_e2'] - df.loc[xgood, 'obs_e2']
            xdT = df.loc[xgood, 'psfex_T'] - df.loc[xgood, 'obs_T']
            logger.info('de1 = %s, %s, %s, %s',np.min(xde1),np.max(xde1),np.mean(xde1),np.std(xde1))
            logger.info('de2 = %s, %s, %s, %s',np.min(xde2),np.max(xde2),np.mean(xde2),np.std(xde2))
            logger.info('dT = %s, %s, %s, %s',np.min(xdT),np.max(xdT),np.mean(xdT),np.std(xdT))

        if args.rm_files:
            logger.info('removing temp files')
            remove_temp_files(sdir, root, keep_files, logger)

    except CatastrophicFailure:
        logger.info('Catastrophic failure')
        logger.info('Log this in the blacklist and continue.')
        # A appropriate flag has already been set.
    except Exception as e:
        logger.info('Caught exception: %s',e)
        logger.info(traceback.format_exc())
        logger.info('Log this in the blacklist and continue.')
        flag |= ERROR_FLAG
    else:
        row['obs_mean_e1'] = np.mean(df.loc[df['obs_flag']==0, 'obs_e1'])
        row['obs_mean_e2'] = np.mean(df.loc[df['obs_flag']==0, 'obs_e2'])
        row['obs_mean_T'] = np.mean(df.loc[df['obs_flag']==0, 'obs_T'])
        row['obs_std_e1'] = np.std(df.loc[df['obs_flag']==0, 'obs_e1'])
        row['obs_std_e2'] = np.std(df.loc[df['obs_flag']==0, 'obs_e2'])
        row['obs_std_T'] = np.std(df.loc[df['obs_flag']==0, 'obs_T'])
        if args.run_piff and (flag & (NO_STARS_FLAG | PSF_FAILURE) == 0):
            row['piff_mean_e1'] = np.mean(df.loc[df['piff_flag']==0, 'piff_e1'])
            row['piff_mean_e2'] = np.mean(df.loc[df['piff_flag']==0, 'piff_e2'])
            row['piff_mean_T'] = np.mean(df.loc[df['piff_flag']==0, 'piff_T'])
            row['piff_std_e1'] = np.std(df.loc[df['piff_flag']==0, 'piff_e1'])
            row['piff_std_e2'] = np.std(df.loc[df['piff_flag']==0, 'piff_e2'])
            row['piff_std_T'] = np.std(df.loc[df['piff_flag']==0, 'piff_T'])
            row['piff_mean_de1'] = np.mean(de1)
            row['piff_mean_de2'] = np.mean(de2)
            row['piff_mean_dT'] = np.mean(dT)
            row['piff_std_de1'] = np.std(de1)
            row['piff_std_de2'] = np.std(de2)
            row['piff_std_dT'] = np.std(dT)
        if args.get_psfex:
            row['psfex_mean_e1'] = np.mean(df.loc[df['psfex_flag']==0, 'psfex_e1'])
            row['psfex_mean_e2'] = np.mean(df.loc[df['psfex_flag']==0, 'psfex_e2'])
            row['psfex_mean_T'] = np.mean(df.loc[df['psfex_flag']==0, 'psfex_T'])
            row['psfex_std_e1'] = np.std(df.loc[df['psfex_flag']==0, 'psfex_e1'])
            row['psfex_std_e2'] = np.std(df.loc[df['psfex_flag']==0, 'psfex_e2'])
            row['psfex_std_T'] = np.std(df.loc[df['psfex_flag']==0, 'psfex_T'])
            row['psfex_mean_de1'] = np.mean(xde1)
            row['psfex_mean_de2'] = np.mean(xde2)
            row['psfex_mean_dT'] = np.mean(xdT)
            row['psfex_std_de1'] = np.std(xde1)
            row['psfex_std_de2'] = np.std(xde2)
            row['psfex_std_dT'] = np.std(xdT)

    row['flag'] = flag
    return df, row


def add_stream_handler(logger, level):
    formatter = logging.Formatter('%(message)s')
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    logger.setLevel(level)  # Not sure if both this and sh.setLevel are required...

def add_file_handler(logger, level, sdir, root):
    log_filename = root + datetime.datetime.now().strftime('%Y-%m-%d-%H%M') + '.log'
    formatter = logging.Formatter('%(asctime)s : %(message)s')
    fh = logging.FileHandler(filename=os.path.join(sdir, log_filename))
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return log_filename

def remove_file_handler(logger, sdir, wdir):
    sdir = os.path.abspath(sdir)
    wdir = os.path.abspath(wdir)
    if sdir != wdir:
        for h in logger.handlers[:]:  # Note slice = copy, so removing items is safe.
            if isinstance(h, logging.FileHandler):
                h.flush()
                h.close()
                time.sleep(1)  # Occasionally this isn't quite written out immediately.
                if sdir in h.baseFilename:
                    final_log_file = h.baseFilename.replace(sdir, wdir)
                    shutil.move(h.baseFilename, final_log_file)
                    logger.removeHandler(h)
                    logger.info('Moved logfile to %s',final_log_file)
                else:
                    logger.removeHandler(h)

def check_T_outliers(df, logger):
    logger.info('Finished processing all CCDs in this exposure.')
    logger.info('flags = %s',df['flag'].values)
    logger.info('mean T = %s',df['obs_mean_T'].values)

    while np.sum(df['flag'] == 0) > 1:
        mask = df['flag'] == 0
        logger.info('ngood = %d',np.sum(mask))
        T = df['piff_mean_T'][mask]  # Use piff_mean_T so we have removed Piff outliers.
        meanT = T.mean()
        stdT = T.std()
        logger.info('mean T = %s,  std T = %s', meanT, stdT)

        i = (T-meanT).abs().idxmax()
        nsig = abs(T[i] - meanT) / stdT
        logger.info('max diff at i=%d.  T=%s,  delta = %f sigma', i, T[i], nsig)

        if nsig < NSIG_T_OUTLIER:
            break
        else:
            df.loc[i, 'flag'] = OUTLIER_SIZE

    logger.info('Done: flags = %s',df['flag'].values)


def main():

    # TODO: Everything is info now.  Should allow different logging levels with -v arg.
    logging_level = logging.INFO
    logger = logging.getLogger('run_piff')
    add_stream_handler(logger, logging_level)

    args = parse_args()
    if args.use_tapebumps:
        tbdata = read_tapebump_file(args.tapebump_file, logger)
    blacklist_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psf'
    if args.tag:
        blacklist_file += '-' + args.tag
    blacklist_file += '.txt'
    if args.blacklist:
        logger.info('Logging blacklisted chips to %s',blacklist_file)

    # Make the work directory if it does not exist yet.
    work = os.path.expanduser(args.work)
    logger.info('work dir = %s',work)
    try:
        if not os.path.exists(work):
            os.makedirs(work)
    except OSError:
        if not os.path.exists(work): raise
    scratch = os.path.expanduser(args.scratch)
    logger.info('scratch dir = %s',scratch)
    try:
        if not os.path.exists(scratch):
            os.makedirs(scratch)
    except OSError:
        if not os.path.exists(scratch): raise

    # A listing Erin made of all the exposures in Y3 used in meds files
    all_exp = fitsio.read(args.base_exposures)
    # Switch to native endians, so pandas doesn't complain.
    all_exp = all_exp.astype(all_exp.dtype.newbyteorder('='))

    if args.file != '':
        logger.info('Read file %s',args.file)
        with open(args.file) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        logger.info('File includes %d exposures',len(exps))
    elif args.exps is not None:
        exps = args.exps
        logger.info('Explicit listing of %d exposures',len(exps))
    else:
        exps = list(set(all_exp['expnum']))
        logger.info('There are a total of %d exposures',len(exps))

    which_zone = fitsio.read(os.path.join(args.pixmappy_dir, 'which_zone.fits'))
    which_zone = pandas.DataFrame(which_zone)
    if sys.version_info >= (3,):
        which_zone['detpos'] = which_zone['detpos'].str.decode("utf-8")

    exps = sorted(exps)
    logger.info('exps = %s',exps)

    for exp in exps:
        exp = int(exp)
        data = all_exp[all_exp['expnum'] == exp]

        logger.info('Start work on exp = %s',exp)

        logger.info('%d images for this exposure',len(data))
        if len(data) == 0:
            logger.info("This exposure is not in Erin's file.  Skip it.")
            continue

        # Make the work and scratch directories for this exposure and clear it if necessary.
        wdir = os.path.join(work,str(exp))
        if args.clear_output:
            if os.path.exists(wdir):
                for f in os.listdir(wdir):
                    try:
                        os.remove(os.path.join(wdir, f))
                    except OSError as e:
                        pass
        try:
            os.makedirs(wdir)
        except OSError:
            if not os.path.exists(wdir): raise
        logger.info('wdir = %s',wdir)

        sdir = os.path.join(scratch,str(exp))
        if (args.clear_output or (args.rm_files and sdir != wdir)) and not args.use_existing:
            if os.path.exists(sdir):
                for f in os.listdir(sdir):
                    try:
                        os.remove(os.path.join(sdir, f))
                    except OSError as e:
                        pass
        try:
            os.makedirs(sdir)
        except OSError:
            if not os.path.exists(sdir): raise
        logger.info('sdir = %s',sdir)

        # Store all information about the exposure in a pandas data frame.
        # Copy over the information from Erin's file.
        exp_info_df = pandas.DataFrame(data)

        if sys.version_info >= (3,):
            for str_key in ['key', 'band', 'path']:
                exp_info_df[str_key] = exp_info_df[str_key].str.decode("utf-8")
        # Add some blank columns to be filled in below.
        for k in ['root', 'date', 'time', 'fits_filter',
                  'image_file', 'bkg_file', 'cat_file', 'star_file', 'piff_file']:
            exp_info_df[k] = np.array([''] * len(data), dtype=str)
        for k in ['zone', 'fits_ccdnum', 'tiling', 'hex',
                  'flag', 'ntot', 'nstars', 'nstars_piff', 'nreserve', 'ngood',]:
            exp_info_df[k] = np.array([-999] * len(data), dtype=int)
        for k in ['sat', 'fits_fwhm', 'star_fwhm',
                  'telra', 'teldec', 'telha', 'airmass', 'sky', 'sigsky',
                  'obs_mean_e1', 'obs_mean_e2', 'obs_mean_T',
                  'piff_mean_e1', 'piff_mean_e2', 'piff_mean_T',
                  'piff_std_e1', 'piff_std_e2', 'piff_std_T',
                  'piff_mean_de1', 'piff_mean_de2', 'piff_mean_dT',
                  'piff_std_de1', 'piff_std_de2', 'piff_std_dT',
                  'piff_chisq', 'piff_dof',
                  'piff_xmin', 'piff_xmax', 'piff_ymin', 'piff_ymax',
                 ]:
            exp_info_df[k] = np.array([-999.] * len(data), dtype=float)
        if args.get_psfex:
            exp_info_df['psfex_file'] = [''] * len(data)
            for k in ['psfex_mean_e1', 'psfex_mean_e2', 'psfex_mean_T',
                      'psfex_std_e1', 'psfex_std_e2', 'psfex_std_T',
                      'psfex_mean_de1', 'psfex_mean_de2', 'psfex_mean_dT',
                      'psfex_std_de1', 'psfex_std_de2', 'psfex_std_dT',
                     ]:
                exp_info_df[k] = np.array([-999.] * len(data), dtype=float)
        for k in ['flag']:
            exp_info_df[k] = np.array([0] * len(data), dtype=int)
        logger.info('dtype = %s',exp_info_df.dtypes)

        exp_info_df.sort_values('ccdnum', inplace=True)
        exp_stars_df = pandas.DataFrame()

        for k, row in exp_info_df.iterrows():

            key, expnum, ccdnum, band = row['key'], row['expnum'], row['ccdnum'], row['band']
            if args.single_ccd and ccdnum != args.single_ccd: continue

            log_file = add_file_handler(logger, logging_level, sdir, '%s_%s_'%(expnum,ccdnum))
            logger.info('\nProcessing %s %s %s %s', key, expnum, ccdnum, band)
            logger.info('Logging for this ccd being written to %s',log_file)

            done = False
            psf_cat_file = os.path.join(wdir, 'psf_cat_%d_%d.fits'%(exp,ccdnum))
            if not args.clear_output and os.path.exists(psf_cat_file):
                logger.info('%s already exists.  Reading the existing file.',psf_cat_file)
                try:
                    with fitsio.FITS(psf_cat_file,'r') as f:
                        all_obj = f['all_obj'].read()
                        stars = f['stars'].read()
                        info = f['info'].read()
                    all_obj = all_obj.astype(all_obj.dtype.newbyteorder('='))
                    stars = stars.astype(stars.dtype.newbyteorder('='))
                    info = info.astype(info.dtype.newbyteorder('='))
                    df = pandas.DataFrame(all_obj)
                    stars = pandas.DataFrame(stars)
                    row = pandas.DataFrame(info).iloc[0]
                    logger.info('Read chip-level psf information from %s',psf_cat_file)
                    logger.info('df = \n%s',df.describe())
                    logger.info('stars = \n%s',stars.describe())
                    logger.info('row = %s',row)
                    done = True
                except Exception as e:
                    logger.error("Caught %r",e)
                    # Leave done = False

            if not done:
                df, row = run_single_ccd(row, args, wdir, sdir, tbdata, which_zone, logger)

                logger.info('row = %s', row)
                # This construction keeps the dtypes of the columns in exp_info_df.
                rowdf = exp_info_df.iloc[0:0].append(row)

                if df is None:
                    logger.info('Catastrophic error for %s, %s', expnum, ccdnum)
                    stars = []
                else:
                    logger.info('all_obj = \n%s', df.describe())
                    star_mask = df['star_flag'] == 1
                    stars = df.loc[star_mask].copy()
                    logger.info('stars = \n%s', stars.describe())
                    with fitsio.FITS(psf_cat_file,'rw',clobber=True) as f:
                        f.write_table(df.to_records(index=False), extname='all_obj')
                        f.write_table(stars.to_records(index=False), extname='stars')
                        f.write_table(rowdf.to_records(index=False), extname='info')
                    logger.info('Wrote chip-level psf information to %s',psf_cat_file)

            # Do this immediately after writing the psf_cat file, so it usually exists in the
            # wdir directory if the psf_cat file exists. (If job times out and is restarted, any
            # log file still in sdir is lost.)
            remove_file_handler(logger, sdir, wdir)

            # row is a copy of the row in the original, so make sure to copy it back.
            exp_info_df.iloc[k] = row

            # Add this chip's stars to the exposure stars df
            flag = row['flag']
            if len(stars) > 0:
                stars.loc[:, 'ccdnum'] = ccdnum
                if 'obs_flag' in stars:
                    stars.loc[:, 'obs_flag'] |= flag * BLACK_FLAG_FACTOR
                if 'piff_flag' in stars:
                    stars.loc[:, 'piff_flag'] |= flag * BLACK_FLAG_FACTOR
                if 'psfex_flag' in stars:
                    stars.loc[:, 'psfex_flag'] |= flag * BLACK_FLAG_FACTOR
                exp_stars_df = exp_stars_df.append(stars, ignore_index=True)

            if flag and args.blacklist:
                log_blacklist(blacklist_file,exp,ccdnum,flag, logger)
                logger.info('Logged flag %d in blacklist',flag)

            if args.single_ccd:
                sys.exit()

        if len(exp_stars_df) == 0:
            logger.error('All CCDs failed for exposure %s',exp)
            continue

        check_T_outliers(exp_info_df, logger)
        logger.info('Done with exposure %s',exp)

        logger.info('exp_stars = \n%s',exp_stars_df.describe())
        logger.info('exp_info = \n%s',exp_info_df.describe())
        exp_cat_file = os.path.join(wdir, 'exp_psf_cat_%d.fits'%exp)
        with fitsio.FITS(exp_cat_file,'rw',clobber=True) as f:
            f.write_table(exp_stars_df.to_records(index=False), extname='stars')
            f.write_table(exp_info_df.to_records(index=False), extname='info')

        logger.info('Wrote exposure information to %s',exp_cat_file)
        logger.info('Done with exposure %s',exp)

    logger.info('\nFinished processing all exposures')


if __name__ == "__main__":
    main()
