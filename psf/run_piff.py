#! /usr/bin/env python
# Run PSFEx for a set of exposures, including making any necessarily input files.
# It also logs errors into a psf blacklist file.

from __future__ import print_function
import os
import sys
import traceback
import numpy
import copy
import glob
import time
import fitsio
import pixmappy
import pandas
import galsim
import galsim.des
import piff

import matplotlib
matplotlib.use('Agg') # needs to be done before import pyplot
import matplotlib.pyplot as plt
plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')

# Don't skip columns in describe output  (default is 20, which is a bit too small)
pandas.options.display.max_columns = 200

# Define the parameters for the blacklist

# How many stars are too few or too many?
FEW_STARS = 20
MANY_STARS_FRAC = 0.5
# How high is a high FWHM?  3.6 arcsec / 0.26 arcsec/pixel = 13.8 pixels
#HIGH_FWHM = 13.8
HIGH_FWHM = 3.6  # (We switched to measuring this in arcsec)

# flag values for blacklist
NO_STARS_FLAG = 1
TOO_FEW_STARS_FLAG = 2
TOO_MANY_STARS_FLAG = 4
TOO_HIGH_FWHM_FLAG = 8
FINDSTARS_FAILURE = 16
PSF_FAILURE = 32
ERROR_FLAG = 64
LARGE_SIZE_SPREAD = 128

# flag values for psf catalog

MAX_CENTROID_SHIFT = 1.0

NOT_USED = 1
BAD_MEASUREMENT = 2
CENTROID_SHIFT = 4
FAILURE = 32
RESERVED = 64
NOT_STAR = 128
BLACK_FLAG_FACTOR = 512 # blacklist flags are this times the original exposure blacklist flag
                        # blacklist flags go up to 64, 


class NoStarsException(Exception):
    pass

def obfuscate(s):
    mask = 'I Have a Dream'
    lmask = len(mask)
    nmask = [ord(c) for c in mask]
    return ''.join([chr(ord(c) ^ nmask[i % lmask]) for i, c in enumerate(s)])

def ps():
    return obfuscate('~\x10+\t\x1f\x15S\x07T3')

def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of exposures')

    # Directory arguments
    parser.add_argument('--sex_dir', default='/astro/u/rarmst/soft/bin/',
                        help='location of sextrator executable')
    parser.add_argument('--piff_exe', default='/astro/u/mjarvis/bin/piffify',
                        help='location of piffify executable')
    parser.add_argument('--findstars_dir', default='/astro/u/mjarvis/bin',
                        help='location wl executables')
    parser.add_argument('--work', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of intermediate outputs')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')
    parser.add_argument('--clear_output', default=False, action='store_const', const=True,
                        help='should the output directory be cleared before writing new files?')

    # Exposure inputs
    parser.add_argument('--file', default='',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--pixmappy', default='zone029', type=str,
                        help='Use the given Pixmappy WCS solution')
    parser.add_argument('--bands', default='grizY', type=str,
                        help='Limit to the given bands')

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
                        default='/astro/u/mjarvis/rmjarvis/DESWL/psf/sex.param_psfex',
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
    parser.add_argument('--rm_files', default=1, type=int,
                        help='remove unpacked files after finished')
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
    parser.add_argument('--nbright_stars', default=10, type=int,
                        help='use median of this many brightest stars for min mag')
    parser.add_argument('--max_mag', default=-1, type=float,
                        help='only use stars brighter than this mag')
    parser.add_argument('--use_tapebumps', default=1, type=int,
                        help='avoid stars in or near tape bumps')
    parser.add_argument('--tapebump_extra', default=2, type=float,
                        help='How much extra room around tape bumps to exclude stars in units of FWHM')
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')
    parser.add_argument('--reserve', default=0, type=float,
                        help='Reserve some fraction of the good stars for testing')
    parser.add_argument('--get_psfex', default=False, action='store_const', const=True,
                        help='Download the PSFEx files along the way')
    parser.add_argument('--plot_fs', default=False, action='store_const', const=True,
                        help='Make a size-magnitude plot of the findstars output')

    args = parser.parse_args()
    return args


def read_tapebump_file(file_name):
    """Read and parse the tapebump file if we are going to need it.
    """
    raw_tbdata = numpy.genfromtxt(file_name, delimiter=',')
    # repackage this as  dict (key = ccdnum) of lists of tuples (ymin, xmin, ymax, xmax)
    tbdata = {}
    for d in raw_tbdata:
        ccdnum = int(d[0])
        if ccdnum not in tbdata:
            tbdata[ccdnum] = []
        tbdata[ccdnum].append( (int(d[1]), int(d[2]), int(d[3]), int(d[4])) )
    print('read in tape bump file.  %d bumps for %d chips'%(len(raw_tbdata),len(tbdata)))
    return tbdata

def exclude_tapebumps(df, tbd, extra):
    """
    Process the tapebump data for a particular chip.

    tbd = tbdata for a particular ccdnum,
    data = the input data for the stars
    extra = how much extra distance around the tape bumps to exclude stars in pixels
    """

def log_blacklist(blacklist_file, exp, ccdnum, flag):
    try:
        with open(blacklist_file,'a') as f:
            f.write("%s %d %d\n"%(exp,ccdnum,flag))
    except OSError as e:
        print(e)
        print('Error opening blacklist.  Wait and try again.')
        import time
        time.sleep(1)
        return log_blacklist(blacklist_file,exp,ccdnum,flag)


def unpack_file(file_name):
    """Create the unpacked file in the work directory if necessary.

    If the unpacked file already exists, then a link is made.
    Otherwise funpack is run, outputting the result into the work directory.
    """
    print('unpack ',file_name)

    img_file = os.path.splitext(file_name)[0]
    if os.path.lexists(img_file):
        print('   %s exists already.  Removing.'%img_file)
        os.remove(img_file)
    print('   unpacking fz file')
    cmd = 'funpack -O {outf} {inf}'.format(outf=img_file, inf=file_name)
    print(cmd)
    for itry in range(5):
        run_with_timeout(cmd, 120)
        if os.path.lexists(img_file): break
        print('%s was not properly made.  Retrying.'%img_file)
        time.sleep(10)

    if not os.path.lexists(img_file):
        print('Unable to create %s.  Skip this file.'%img_file)
        return None

    return img_file

def read_image_header(row, img_file):
    """Read some information from the image header and write into the df row.
    """
    hdu = 1

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
        telra = galsim.HMS_Angle(telra) / galsim.degrees
        teldec = galsim.DMS_Angle(teldec) / galsim.degrees
        telha = galsim.HMS_Angle(telha) / galsim.degrees

        airmass = float(h.get('AIRMASS',-999))
        sky = float(h.get('SKYBRITE',-999))
        sigsky = float(h.get('SKYSIGMA',-999))

        tiling = int(h.get('TILING',0))
        hex = int(h.get('HEX',0))

    except:
        print("Cannot read header information from " + img_file)
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
                   sex_dir, sex_config, sex_params, sex_filter, sex_nnw):
    """Run sextractor, but only if the output file does not exist yet.
    """
    cat_file = os.path.join(wdir,root+'_cat.fits')

    print('   running sextractor')
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
    print(cat_cmd)
    run_with_timeout(cat_cmd, 120)

    if not os.path.exists(cat_file) or os.path.getsize(cat_file) == 0:
        print('   Error running SExtractor.  No ouput file was written.')
        print('   Try again, in case it was a fluke.')
        run_with_timeout(cat_cmd, 120)
        if os.path.getsize(psf_file) == 0:
            os.remove(cat_file)
            print('   Error running SExtractor (again).')
            return None
    return cat_file

def run_findstars(row, wdir, fs_dir, fs_config):
    """Run findstars
    """
    star_file = row['star_file']

    # run find stars
    print('   running findstars')
    findstars_cmd = '{fs_dir}/findstars {fs_config} root={root} cat_file={cat_file} stars_file={star_file} input_prefix={wdir}/'.format(
            fs_dir=fs_dir, fs_config=fs_config, root=row['root'], cat_file=row['cat_file'],
            star_file=star_file, wdir=wdir)
    print(findstars_cmd)
    run_with_timeout(findstars_cmd, 120)

    if not os.path.exists(star_file) or os.path.getsize(star_file) == 0:
        print('   Error running findstars.  Rerun with verbose=2.')
        findstars_cmd = findstars_cmd + ' verbose=2 debug_ext=_fs.debug'
        print(findstars_cmd)
        run_with_timeout(findstars_cmd, 240)
        print('   The debug file is',row['root'] + '_fs.debug')
        if not os.path.exists(star_file) or os.path.getsize(star_file) == 0:
            return None
    return star_file

def read_findstars(star_file, cat_file):
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
    print('   found %d stars'%nstars)
    if nstars == 0:
        # Can't really do any of the rest of this, so skip out to the end.
        raise NoStarsException()

    # Add on some extra information from the sextractor catalog
    sdata = fitsio.read(cat_file, 2)
    assert len(data) == len(sdata)
    df['ra'] = sdata['ALPHAWIN_J2000']
    df['dec'] = sdata['DELTAWIN_J2000']

    return df


def remove_bad_stars(df, ccdnum, tbdata,
                     mag_cut, nbright_stars, max_mag,
                     use_tapebumps, tapebump_extra, reserve, fwhm):
    """Remove stars that are considered bad for some reason.

    Currently these reasons include:
    - Magnitude indicates that the star is significantly contaminated by the brighter/fatter
      effect.
    - Star falls in or near the tape bumps.
    """

    use = df['star_flag'] == 1

    if mag_cut > 0:
        mags = numpy.sort(df['mag'])
        min_star = numpy.median(mags[0:nbright_stars])
        print('   min mag = ',mags[0])
        print('   median of brightest %d is '%nbright_stars, min_star)

        use = use & (df['mag'] > min_star + mag_cut)
        print('   select stars dimmer than',min_star+mag_cut)
        print('   which brings star count to ',use.sum())

    if max_mag > 0:
        use = use & (df['mag'] < max_mag)
        print('   also select stars brighter than',max_mag)
        print('   which brings star count to ',use.sum())

    if use_tapebumps:
        # I'm sure it doesn't matter, but add an extra 0.5 pixel to the slop because the tape bump
        # values are in terms of which pixels to include as part of the tape bump.  So the edges
        # are an extra 0.5 pixel outside of that.
        extra = tapebump_extra * fwhm + 0.5
        x = df['x']
        y = df['y']
        tbd = tbdata[ccdnum]
        mask = [(y>tb[0]-extra) & (x>tb[1]-extra) & (y<tb[2]+extra) & (x<tb[3]+extra) for tb in tbd]
        mask = numpy.any(mask, axis=0)
        use = use & ~mask
        print('   excluding tape bumps brings star count to ',use.sum())

    if reserve:
        print('   reserve ',reserve)
        n = use.sum()
        perm = numpy.random.permutation(n)
        n1 = int(reserve * n)
        print('   reserving',n1)
        print('   initial ids = ',df[use]['id'].values)
        # There is surely a more efficient way to do this, but I kept getting confused about
        # when numpy makes a copy vs a view.  This step-by-step calculation works.
        r1 = numpy.zeros((n), dtype=bool)
        r1[:n1] = True
        r2 = numpy.zeros((n), dtype=bool)
        r2[perm] = r1
        reserve = numpy.zeros_like(use, dtype=bool)
        reserve[use] = r2
        print('   nreserve = ',numpy.sum(reserve))
        df['reserve'] = reserve
        print('   reserve ids = ',df['id'][df['reserve']].values)
        print('   final ids = ',df['id'][~df['reserve']].values)
        use = use & ~reserve
        print('   after reserve: nstars = ',use.sum())

    df['use'] = use
    return use.sum(), len(use)


def get_fwhm(df):
    """Get the fwhm from the SExtractor FLUX_RADIUS estimates.
    """

    # get the brightest 10 stars that have flags=0 and take the median just in case some
    # strange magnitudes were selected
    use = df['use']
    fwhm = 2. * df[use]['sigma0']  # 2 * sigma is approx fwhm
    stats = ( numpy.min(fwhm), numpy.max(fwhm),
              numpy.mean(fwhm), numpy.median(fwhm) )
    return stats


def plot_fs(df, filename):
    import numpy
    fig, ax = plt.subplots(1,1)

    # Overall zeropoint adjustment from Eli.
    zp = 5.3

    size = 2. * df['sigma0']**2
    mag = df['mag']

    # Make a simpler flag to use below:
    # 0 = good detection
    # 1 = star candidate
    # 2 = final star
    flag = df['star_flag'].copy()
    flag[df['piff_flag']==0] = 2
    flag[df['piff_flag']==RESERVED] = 2  # Also count the reserve stars as psf stars
    print('num detections: ',numpy.sum(flag == 0))
    print('num candidates: ',numpy.sum(flag == 1))
    print('num final psf stars: ',numpy.sum(flag == 2))

    # Mag adjustment from Eli.
    mag = mag + zp

    ax.set_xlim(10+zp, 24.)
    ax.set_ylim(0., 1.2)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel(r'$T = 2\sigma^2$ (arcsec${}^2$)')

    index = numpy.argsort(mag)
    n = numpy.arange(len(mag))

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

    plt.tight_layout()
    plt.savefig(filename)
    print('Wrote fs plot to ',filename)
    plt.close(fig)

    numpy.savetxt(os.path.splitext(filename)[0] + '.dat',
                  numpy.array(zip(size, mag, flag), dtype='f8, f8, i2'), fmt='%r',
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
             pixmappy, exp, ccdnum):
    """Run Piffify

    Returns True if successful, False if there was a catastrophic failure and no output 
    file was written.
    """
    if not os.path.exists(piff_config) and '/' not in piff_config:
        piff_config = os.path.join('/astro/u/mjarvis/rmjarvis/DESWL/psf/',piff_config)
    if os.path.lexists(psf_file):
        print('   deleting existing',psf_file)
        os.unlink(psf_file)

    print('   making cat file for piff')
    stars_df = df[df['use']]
    print('stars_df = \n',stars_df.describe())
    piff_cat_file = cat_file.replace('stars','use_stars')
    fitsio.write(piff_cat_file, stars_df.to_records(index=False), clobber=True)

    print('   running piff')
    piff_cmd = '{piff_exe} {config} input.image_file_name={image} input.cat_file_name={cat} output.file_name={psf}'.format(
            piff_exe=piff_exe, config=piff_config, image=img_file, cat=piff_cat_file, psf=psf_file)
    piff_cmd += ' input.wcs.file_name={pixmappy} input.wcs.exp={exp} input.wcs.ccdnum={ccdnum}'.format(
            pixmappy=pixmappy, exp=exp, ccdnum=ccdnum)
    print(piff_cmd)
    run_with_timeout(piff_cmd, 300)  # 5 minutes should be way more than plenty!

    if not os.path.exists(psf_file) or os.path.getsize(psf_file) == 0:
        print('   Error running Piff.  No ouput file was written.')
        print('   Try again, in case it was a fluke.')
        piff_cmd += ' verbose=3'  # Add more debugging so we can see what might have gone wrong.
        run_with_timeout(piff_cmd, 600)  # And double the time just in case that was the problem.
        if not os.path.exists(psf_file) or os.path.getsize(psf_file) == 0:
            print('   Error running Piff (again).')
            return False
    return True

def remove_temp_files(wdir, root, keep_files):
    """Remove wdir/root* except for any files listed in the keep_files list
    """
    files = sorted(glob.glob('%s/%s*'%(wdir,root)))
    for save in keep_files:
        if save in files:
            files.remove(save)
        else:
            print('WARNING: %s not found in %s'%(save,wdir))

    print('   Removing the following files from ',wdir)
    for f in files:
        print('       ',os.path.split(f)[1])
        os.remove(f)
    print('   Done')


def wget(url_base, path, wdir, file):
    url = url_base + path + file
    full_file = os.path.join(wdir,file)

    import wget
    if not os.path.isfile(full_file):
        print('Downloading ',full_file)
        # Sometimes this fails with an "http protocol error, bad status line".
        # Maybe from too many requests at once or something.  So we retry up to 5 times.
        nattempts = 5
        for attempt in range(1,nattempts+1):
            print('wget %s  (attempt %d)'%(url, attempt))
            try:
                wget.download(url, bar=None, out=full_file)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print('Caught ',e)
                if attempt < nattempts:
                    print('Try again.')
                    import time
                    time.sleep(2)
                continue
            else:
                break
    return full_file


def hsm(im, wt=None):
    flag = 0
    try:
        shape_data = im.FindAdaptiveMom(weight=wt, strict=False)
    except:
        print(' *** Bad measurement (caught exception).  Mask this one.')
        flag |= BAD_MEASUREMENT

    if shape_data.moments_status != 0:
        print('status = ',shape_data.moments_status)
        print(' *** Bad measurement.  Mask this one.')
        flag |= BAD_MEASUREMENT

    dx = shape_data.moments_centroid.x - im.trueCenter().x
    dy = shape_data.moments_centroid.y - im.trueCenter().y
    if dx**2 + dy**2 > MAX_CENTROID_SHIFT**2:
        print(' *** Centroid shifted by ',dx,dy,'.  Mask this one.')
        flag |= CENTROID_SHIFT

    # Account for the image wcs
    if im.wcs.isPixelScale():
        g1 = shape_data.observed_shape.g1
        g2 = shape_data.observed_shape.g2
        T = 2 * shape_data.moments_sigma**2 * im.scale**2
    else:
        e1 = shape_data.observed_shape.e1
        e2 = shape_data.observed_shape.e2
        s = shape_data.moments_sigma

        jac = im.wcs.jacobian(im.trueCenter())
        M = numpy.matrix( [[ 1 + e1, e2 ], [ e2, 1 - e1 ]] ) * s*s
        J = jac.getMatrix()
        M = J * M * J.T

        e1 = (M[0,0] - M[1,1]) / (M[0,0] + M[1,1])
        e2 = (2.*M[0,1]) / (M[0,0] + M[1,1])
        T = M[0,0] + M[1,1]

        shear = galsim.Shear(e1=e1, e2=e2)
        g1 = shear.g1
        g2 = shear.g2

    return g1, g2, T, flag

def measure_star_shapes(df, image_file, bkg_file, noweight, wcs):
    """Measure shapes of the raw stellar images at each location.
    """
    print('Read in stars in file: ',image_file)

    ind = df.index[df['star_flag'] == 1]
    print('ind = ',ind)
    n_psf = len(ind)
    print('n_psf = ',n_psf)

    df['obs_e1'] = [ -999. ] * len(df)
    df['obs_e2'] = [ -999. ] * len(df)
    df['obs_T'] = [ -999. ] * len(df)
    df['obs_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'obs_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'obs_flag'] |= RESERVED
    df.loc[~df['use'], 'obs_flag'] |= NOT_USED

    full_image = galsim.fits.read(image_file, hdu=1)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        full_weight = galsim.fits.read(image_file, hdu=3)
        full_weight.array[full_weight.array < 0] = 0.

    # Subtract off the background
    if bkg_file is not None:
        bkg = galsim.fits.read(bkg_file)
        print('median = ',numpy.median(bkg.array))
        full_image -= bkg

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]

        print('Measure shape for star at ',x,y)
        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2, 
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds

        im = full_image[b]
        if noweight:
            wt = None
        else:
            wt = full_weight[b]

        # If don't have a bkg file, this is probably ok.
        if bkg_file is None:
            im -= df['sky']

        e1, e2, T, flag = hsm(im, wt)
        df.loc[i, 'obs_e1'] = e1
        df.loc[i, 'obs_e2'] = e2
        df.loc[i, 'obs_T'] = T
        df.loc[i, 'obs_flag'] |= flag
    print('final obs_flag = ',df['obs_flag'][ind].values)
    print('df[ind] = ',df.loc[ind].describe())

def find_index(x1, y1, x2, y2):
    """Find the index of the closest point in (x2,y2) to each (x1,y1)

    Any points that do not have a corresponding point within 1 arcsec gets index = -1.
    """
    index = numpy.zeros(len(x1),dtype=int)
    for i in range(len(x1)):
        close = numpy.where( (x2 > x1[i]-1.) &
                             (x2 < x1[i]+1.) &
                             (y2 > y1[i]-1.) &
                             (y2 < y1[i]+1.) )[0]
        if len(close) == 0:
            #print('Could not find object near x,y = (%f,%f)'%(x,y))
            index[i] = -1
        elif len(close) == 1:
            index[i] = close[0]
        else:
            print('Multiple objects found near x,y = (%f,%f)'%(x1[i],y1[i]))
            amin = numpy.argmin((x2[close] - x1[i])**2 + (y2[close] - y1[i])**2)
            print('Found minimum at ',close[amin],', where (x,y) = (%f,%f)'%(
                    x2[close[amin]], y2[close[amin]]))
            index[i] = close[amin]
    return index


def measure_piff_shapes(df, psf_file, image_file, wcs):
    """Measure shapes of the Piff solution at each location.
    """
    print('Read in Piff file: ',psf_file)

    ind = df.index[df['star_flag'] == 1]
    print('ind = ',ind)
    n_psf = len(ind)
    print('n_psf = ',n_psf)

    df['piff_e1'] = [ -999. ] * len(df)
    df['piff_e2'] = [ -999. ] * len(df)
    df['piff_T'] = [ -999. ] * len(df)
    df['piff_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'piff_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'piff_flag'] |= RESERVED
    df.loc[~df['use'], 'piff_flag'] |= NOT_USED

    # Flag the stars not used by Piff to make the solution.
    used_data = fitsio.read(psf_file, ext='psf_stars')
    print('Piff used %d stars for the solution'%len(used_data))
    used_xmin = used_data['x'].min()
    used_xmax = used_data['x'].max()
    used_ymin = used_data['y'].min()
    used_ymax = used_data['y'].max()
    print('   final bounds of used stars = ',used_xmin,used_xmax,used_ymin,used_ymax)
    used_area = (used_xmax-used_xmin)*(used_ymax-used_ymin)
    print('   area = ',used_area)
    used_index = find_index(df['x'], df['y'], used_data['x'], used_data['y'])
    df.loc[used_index < 0, 'piff_flag'] |= NOT_USED

    try:
        psf = piff.read(psf_file)
    except Exception as e:
        print('Caught ',e)
        df.loc[ind, 'piff_flag'] = FAILURE
        return

    full_image = galsim.fits.read(image_file, hdu=1)

    if wcs is not None:
        full_image.wcs = wcs

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]
        print('Measure Piff model shape at ',x,y)

        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2, 
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        im = psf.draw(x=x, y=y, image=im)

        e1, e2, T, flag = hsm(im)
        df.loc[i, 'piff_e1'] = e1
        df.loc[i, 'piff_e2'] = e2
        df.loc[i, 'piff_T'] = T
        df.loc[i, 'piff_flag'] |= flag
    print('final piff_flag = ',df['piff_flag'][ind].values)
    print('df[ind] = ',df.loc[ind].describe())


def measure_psfex_shapes(df, psfex_file, image_file, wcs):
    """Measure shapes of the PSFEx solution at each location.
    """
    print('Read in PSFEx file: ',psfex_file)

    ind = df.index[df['star_flag'] == 1]
    print('ind = ',ind)
    n_psf = len(ind)
    print('n_psf = ',n_psf)

    df['psfex_e1'] = [ -999. ] * len(df)
    df['psfex_e2'] = [ -999. ] * len(df)
    df['psfex_T'] = [ -999. ] * len(df)
    df['psfex_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'psfex_flag'] = 0

    if 'reserve' in df:
        df.loc[df['reserve'], 'psfex_flag'] |= RESERVED
    df.loc[~df['use'], 'psfex_flag'] |= NOT_USED

    try:
        psf = galsim.des.DES_PSFEx(psfex_file, image_file)
    except Exception as e:
        print('Caught ',e)
        df.loc[ind, 'psfex_flag'] = FAILURE
        return

    full_image = galsim.fits.read(image_file, hdu=1)

    if wcs is not None:
        full_image.wcs = wcs

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]
        print('Measure PSFEx model shape at ',x,y)
        image_pos = galsim.PositionD(x,y)
        psf_i = psf.getPSF(image_pos)

        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2, 
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        im = psf_i.drawImage(image=im, method='no_pixel')

        e1, e2, T, flag = hsm(im)
        df.loc[i, 'psfex_e1'] = e1
        df.loc[i, 'psfex_e2'] = e2
        df.loc[i, 'psfex_T'] = T
        df.loc[i, 'psfex_flag'] |= flag
    print('final psfex_flag = ',df['psfex_flag'][ind].values)
    print('df[ind] = ',df.loc[ind].describe())


def main():
    args = parse_args()
    if args.use_tapebumps:
        tbdata = read_tapebump_file(args.tapebump_file)
    blacklist_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psf'
    if args.tag:
        blacklist_file += '-' + args.tag
    blacklist_file += '.txt'
    if args.blacklist:
        print('Logging blacklisted chips to',blacklist_file)


    # Make the work directory if it does not exist yet.
    work = os.path.expanduser(args.work)
    print('work dir = ',work)
    try:
        if not os.path.exists(work):
            os.makedirs(work)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    
    # The url to use up to just before OPS
    url_base = 'https://rmjarvis:%s@desar2.cosmology.illinois.edu/DESFiles/desarchive/'%ps()

    # A listing Erin made of all the exposures in Y3 used in meds files
    all_exp = fitsio.read('/astro/u/mjarvis/work/y3_piff/exposures-ccds-Y3A1_COADD.fits')
    # Switch to native endians, so pandas doesn't complain.
    all_exp = all_exp.astype(all_exp.dtype.newbyteorder('='))

    if args.file != '':
        print('Read file ',args.file)
        with open(args.file) as fin:
            exps = [ line for line in fin if line[0] != '#' ]
        print('File includes %d exposures'%len(exps))
    elif args.exps is not None:
        exps = args.exps
        print('Explicit listing of %d exposures'%len(exps))
    else:
        exps = list(set(all_exp['expnum']))
        print('There are a total of %d exposures'%len(exps))

    if args.pixmappy is not None:
        pixmappy_wcs = pixmappy.PixelMapCollection(args.pixmappy).wcs
        wcs_exps = [k[1:7] for k in pixmappy_wcs.keys() if k[0]=='D']
        print('Pixmappy file has %d exposures'%len(wcs_exps))
        exps = [exp for exp in exps if exp in wcs_exps]
        print('Limiting to the %d exposures in pixmappy wcs file'%len(exps))

    exps = sorted(exps)
    print('exps = ',exps)

    for exp in exps:
        exp = int(exp)
        data = all_exp[all_exp['expnum'] == exp]

        # TODO: trim to the listed bands above rather than here.
        if data['band'][0] not in args.bands:  # (All bands should be equal for same exp)
            print('Skipping %d because band %s not in %s'%(exp,data['band'][0],args.bands))
            continue

        print('Start work on exp = ',exp)
        print('%d images for this exposure'%len(data))

        # Store all information about the exposure in a pandas data frame.
        # Copy over the information from Erin's file.
        exp_df = pandas.DataFrame(data)
        # Add some blank columns to be filled in below.
        for k in ['root', 'image_file', 'bkg_file', 'psfex_file', 'cat_file', 'star_file',
                  'piff_file']:
            exp_df[k] = [''] * len(data)
        for k in ['sat', 'fits_fwhm',
                  'obs_mean_e1', 'obs_mean_e2', 'obs_mean_T',
                  'piff_mean_de1', 'piff_mean_de2', 'piff_mean_dT',
                  'piff_std_de1', 'piff_std_de2', 'piff_std_dT',
                 ]:
            exp_df[k] = [-999.] * len(data)
        if args.get_psfex:
            for k in ['psfex_mean_de1', 'psfex_mean_de2', 'psfex_mean_dT',
                      'psfex_std_de1', 'psfex_std_de2', 'psfex_std_dT',
                     ]:
                exp_df[k] = [-999.] * len(data)
        for k in ['flag']:
            exp_df[k] = [0] * len(data)

        # Make the work directory for this exposure and clear it if necessary.
        wdir = os.path.join(work,str(exp))
        if args.clear_output:
            import shutil
            if os.path.exists(wdir):
                for f in os.listdir(wdir):
                    try:
                        os.remove(os.path.join(wdir, f))
                    except OSError as e:
                        print("Ignore OSError from remove(wdir/f):")
                        print(e)
                        pass
        try:
            os.makedirs(wdir)
        except:
            if not os.path.exists(wdir): raise
        print('wdir = ',wdir)

        exp_df.sort_values('ccdnum', inplace=True)

        for k, row in exp_df.iterrows():
            key, expnum, ccdnum, band = row['key'], row['expnum'], row['ccdnum'], row['band']
            print('\nProcessing ', key, expnum, ccdnum, band)
            path = row['path'].strip()
            print('path = ',path)

            # Store all information about the stars in a pandas data frame.
            df = pandas.DataFrame()
            flag = 0

            keep_files = []

            try:

                # Download the files we need:
                base_path, _, _, image_file_name = path.rsplit('/',3)
                root, ext = image_file_name.rsplit('_',1)
                print('root, ext = |%s| |%s|'%(root,ext))
                image_file = wget(url_base, base_path + '/red/immask/', wdir, root + '_' + ext)
                print('image_file = ',image_file)
                row['root'] = root
                row['image_file'] = image_file

                bkg_file = wget(url_base, base_path + '/red/bkg/', wdir, root + '_bkg.fits.fz')
                row['bkg_file'] = bkg_file

                if args.get_psfex:
                    psfex_file = wget(url_base, base_path + '/psf/', wdir, root + '_psfexcat.psf')
                    row['psfex_file'] = psfex_file
                    keep_files.append(psfex_file)

                read_image_header(row, image_file)
                sat = row['sat']
                fits_fwhm = row['fits_fwhm']

                # Run Sextractor
                cat_file = os.path.join(wdir, root + '_cat.fits')
                # Unfortunately, the desdm catalog in the cat directory doesn't seem to be
                # complete.  Looks like it only has stars maybe?  Which is nominally fine, but
                # we're not sure if we trust their star selection.  So we run sextractor ourself.
                if args.run_sextractor or not os.path.isfile(cat_file):
                    # Unpack the image file if necessary
                    unpack_image_file = unpack_file(image_file)
                    if unpack_image_file is None:
                        # This was our signal to skip this without blacklisting.  Just continue.
                        flag |= ERROR_FLAG
                        raise NoStarsException()
                    # Also need the fwhm for doing the tape bumps.
                    cat_file = run_sextractor(wdir, root, unpack_image_file, sat, fits_fwhm,
                            args.noweight, args.sex_dir, args.sex_config, args.sex_params,
                            args.sex_filter, args.sex_nnw)
                    if cat_file is None:
                        raise NoStarsException()
                print('cat_file = ',cat_file)
                row['cat_file'] = cat_file

                # Run findstars
                star_file = os.path.join(wdir, root + '_stars.fits')
                row['star_file'] = star_file
                if args.run_findstars or not os.path.isfile(star_file):
                    star_file = run_findstars(row, wdir, args.findstars_dir, args.findstars_config)

                # Read the findstars catalog (need to do this even if star_file exists)
                df = read_findstars(star_file, cat_file)
                if df is None:
                    print('     -- flag for findstars failure')
                    flag |= FINDSTARS_FAILURE
                    raise NoStarsException()

                # Cut the brighest magnitudes or other exclusions/reservations
                nstars, ntot = remove_bad_stars(
                        df, ccdnum, tbdata,
                        args.mag_cut, args.nbright_stars, args.max_mag,
                        args.use_tapebumps, args.tapebump_extra, args.reserve, fits_fwhm)

                # Check if there are few or many staras.
                if nstars < FEW_STARS:
                    print('     -- flag for too few stars: ',nstars)
                    flag |= TOO_FEW_STARS_FLAG
                if nstars <= 1:
                    raise NoStarsException()
                if nstars > MANY_STARS_FRAC * ntot:
                    print('     -- flag for too many stars: %d/%d'%(nstars,ntot))
                    flag |= TOO_MANY_STARS_FLAG

                # Get the median fwhm of the given stars
                # Returns min, max, mean, median.  We use median, which is index 3.
                star_fwhm = get_fwhm(df) 
                print('   fwhm of stars = ',star_fwhm)
                print('   cf. header fwhm = ',fits_fwhm * 0.26)
                if star_fwhm[3] > HIGH_FWHM:
                    print('     -- flag for too high fwhm')
                    flag |= TOO_HIGH_FWHM_FLAG
                if star_fwhm[3] > 1.5 * fits_fwhm * 0.26:
                    print('     -- flag for too high fwhm compared to fwhm from fits header')
                    flag |= TOO_HIGH_FWHM_FLAG
                row['star_fwhm'] = star_fwhm[3]

                # Get the pixmappy wcs for this ccd for correcting the shapes
                if args.pixmappy is not None:
                    wcs = pixmappy.GalSimWCS(args.pixmappy, exp=exp, ccdnum=ccdnum)
                    wcs._color = 0.  # For now.  Revisit when doing color-dependent PSF.
                else:
                    wcs = None

                # Measure the shpes and sizes of the stars
                measure_star_shapes(df, image_file, bkg_file, args.noweight, wcs)

                # Another check is that the spread in star sizes isn't too large
                obs_T = df.loc[df['obs_flag']==0, 'obs_T']
                if numpy.std(obs_T) > 0.15 * numpy.mean(obs_T):
                    print('     -- flag for large size spread: %f > 0.15 * %f'%(
                            numpy.std(obs_T), numpy.mean(obs_T)))
                    flag |= LARGE_SIZE_SPREAD

                # Run piff
                psf_file = os.path.join(wdir, root + '_piff.fits')
                keep_files.append(psf_file)
                if args.run_piff:
                    success = run_piff(df, image_file, star_file, psf_file,
                                       args.piff_exe, args.piff_config,
                                       args.pixmappy, exp, ccdnum)
                    if not success:
                        flag |= PSF_FAILURE

                # Measure the psf shapes
                if not (flag & PSF_FAILURE):
                    measure_piff_shapes(df, psf_file, image_file, wcs)
                    good = (df['piff_flag'] == 0) & (df['obs_flag'] == 0)
                    de1 = df.loc[good, 'piff_e1'] - df.loc[good, 'obs_e1']
                    de2 = df.loc[good, 'piff_e2'] - df.loc[good, 'obs_e2']
                    dT = df.loc[good, 'piff_T'] - df.loc[good, 'obs_T']
                    print('de1 = ',numpy.min(de1),numpy.max(de1),numpy.mean(de1),numpy.std(de1))
                    print('de2 = ',numpy.min(de2),numpy.max(de2),numpy.mean(de2),numpy.std(de2))
                    print('dT = ',numpy.min(dT),numpy.max(dT),numpy.mean(dT),numpy.std(dT))
                else:
                    de1 = -999.
                    de2 = -999.
                    dT = -999.

                if args.plot_fs:
                    fs_plot_file = os.path.join(wdir, root + '_fs.pdf')
                    plot_fs(df, fs_plot_file)
                    keep_files.append(fs_plot_file)

                if args.get_psfex:
                    measure_psfex_shapes(df, psfex_file, image_file, wcs)
                    xgood = (df['psfex_flag'] == 0) & (df['obs_flag'] == 0)
                    xde1 = df.loc[xgood, 'psfex_e1'] - df.loc[xgood, 'obs_e1']
                    xde2 = df.loc[xgood, 'psfex_e2'] - df.loc[xgood, 'obs_e2']
                    xdT = df.loc[xgood, 'psfex_T'] - df.loc[xgood, 'obs_T']
                    print('de1 = ',numpy.min(xde1),numpy.max(xde1),numpy.mean(xde1),numpy.std(xde1))
                    print('de2 = ',numpy.min(xde2),numpy.max(xde2),numpy.mean(xde2),numpy.std(xde2))
                    print('dT = ',numpy.min(xdT),numpy.max(xdT),numpy.mean(xdT),numpy.std(xdT))

                if args.rm_files:
                    print('removing temp files')
                    remove_temp_files(wdir, root, keep_files)


            except NoStarsException:
                print('No stars.  Log this in the blacklist and continue.')
                flag |= NO_STARS_FLAG
            except Exception as e:
                print('Caught exception: ',e)
                traceback.print_exc()
                print('Log this in the blacklist and continue.')
                flag |= ERROR_FLAG
            else:
                print('df = \n',df.describe())
                psf_info_file = os.path.join(wdir, 'psf_cat_%d_%d.fits'%(exp,ccdnum))
                fitsio.write(psf_info_file, df.to_records(index=False), clobber=True)
                print('Wrote psf information to ',psf_info_file)
                row['obs_mean_e1'] = numpy.mean(df.loc[df['obs_flag']==0, 'obs_e1'])
                row['obs_mean_e2'] = numpy.mean(df.loc[df['obs_flag']==0, 'obs_e2'])
                row['obs_mean_T'] = numpy.mean(df.loc[df['obs_flag']==0, 'obs_T'])
                row['obs_std_e1'] = numpy.std(df.loc[df['obs_flag']==0, 'obs_e1'])
                row['obs_std_e2'] = numpy.std(df.loc[df['obs_flag']==0, 'obs_e2'])
                row['obs_std_T'] = numpy.std(df.loc[df['obs_flag']==0, 'obs_T'])
                row['piff_mean_de1'] = numpy.mean(de1)
                row['piff_mean_de2'] = numpy.mean(de2)
                row['piff_mean_dT'] = numpy.mean(dT)
                row['piff_std_de1'] = numpy.std(de1)
                row['piff_std_de2'] = numpy.std(de2)
                row['piff_std_dT'] = numpy.std(dT)
                if args.get_psfex:
                    row['psfex_mean_de1'] = numpy.mean(xde1)
                    row['psfex_mean_de2'] = numpy.mean(xde2)
                    row['psfex_mean_dT'] = numpy.mean(xdT)
                    row['psfex_std_de1'] = numpy.std(xde1)
                    row['psfex_std_de2'] = numpy.std(xde2)
                    row['psfex_std_dT'] = numpy.std(xdT)

            row['flag'] = flag
            # row is a copy of the row in the original, so make sure to copy it back.
            exp_df.iloc[k] = row

            if flag and args.blacklist:
                log_blacklist(blacklist_file,exp,ccdnum,flag)
                print('Logged flag %d in blacklist'%flag)

            if args.single_ccd:
                sys.exit()

        print('Done with exposure ',exp)
        print('exp_df = \n',exp_df.describe())
        exp_info_file = os.path.join(wdir, 'exp_info_%d.fits'%exp)
        fitsio.write(exp_info_file, exp_df.to_records(index=False), clobber=True)
        print('Wrote exposure information to ',exp_info_file)

    print('\nFinished processing all exposures')


if __name__ == "__main__":
    main()
