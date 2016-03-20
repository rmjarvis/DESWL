#! /usr/bin/env python
# Calculate shapes of stars and the shapes of the PSFEx measurements of the stars.

# Define the flag values:

NOT_USED = 1
MEAS_BAD_MEASUREMENT = 2
MEAS_CENTROID_SHIFT = 4
PSFEX_BAD_MEASUREMENT = 8
PSFEX_CENTROID_SHIFT = 16
PSFEX_FAILURE = 32
DESDM_BAD_MEASUREMENT = 64
DESDM_CENTROID_SHIFT = 128
DESDM_FAILURE = 256
DESDM_FLAG_FACTOR = DESDM_BAD_MEASUREMENT / PSFEX_BAD_MEASUREMENT
BLACK_FLAG_FACTOR = 512 # blacklist flags are this times the original exposure blacklist flag
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Build PSF catalogs for a set of runs/exposures')

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


def get_wcs(img_file):
    """Read the wcs from the image header
    """
    import pyfits
    import galsim

    if img_file.endswith('fz'):
        hdu = 1
    else:
        hdu = 0
    wcs = galsim.FitsWCS(img_file, hdu=hdu)
    return wcs
 
def parse_file_name(file_name):
    """Parse the PSFEx file name to get the root name and the chip number
    """
    import os

    dir, base_file = os.path.split(file_name)
    if os.path.splitext(base_file)[1] == '.fz':
        base_file=os.path.splitext(base_file)[0]
    if os.path.splitext(base_file)[1] != '.fits':
        raise ValueError("Invalid file name "+file)
    root = os.path.splitext(base_file)[0]

    ccdnum = int(root.split('_')[-1])
    return dir, root, ccdnum


def read_used(exp_dir, root):
    """Read in the .used.fits file that PSFEx generates with the list of stars that actually
    got used in making the PSFEx file.
    """
    import os
    import astropy.io.fits as pyfits
    import copy

    file_name = os.path.join(exp_dir, root + '_psfcat.used.fits')
    if not os.path.isfile(file_name):
        return None
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
    import os
    import astropy.io.fits as pyfits
    import copy

    file_name = os.path.join(exp_dir, root + '_findstars.fits')
    print 'file_name = ',file_name
    if not os.path.isfile(file_name):
        return None
    try:
        with pyfits.open(file_name, memmap=False) as fp:
            data = copy.copy(fp[1].data)
        # This has the following columns:
        # id: The original id from the SExtractor catalog
        # x: The x position
        # y: The y position
        # sky: The local sky value
        # noise: The estimated noise.  But these are all 0, so I think this isn't being calculated.
        # size_flags: Error flags that occurred when estimating the size
        # mag: The magnitude from SExtractor
        # sg: SExtractor's star/galaxy estimate.  Currently SPREAD_MODEL
        # sigma0: The shapelet sigma that results in a b_11 = 0 shapelet parameter.
        # star_flag: 1 if findstars thought this was a star, 0 otherwise.
        return data
    except Exception as e:
        print 'Caught exception:'
        print e
        return None
 
def find_index(x1, y1, x2, y2):
    """Find the index of the closest point in (x2,y2) to each (x1,y1)

    Any points that do not have a corresponding point within 1 arcsec gets index = -1.
    """
    import numpy
    index = numpy.zeros(len(x1),dtype=int)
    for i in range(len(x1)):
        close = numpy.where( (x2 > x1[i]-1.) &
                             (x2 < x1[i]+1.) &
                             (y2 > y1[i]-1.) &
                             (y2 < y1[i]+1.) )[0]
        if len(close) == 0:
            #print 'Could not find object near x,y = (%f,%f)'%(x,y)
            index[i] = -1
        elif len(close) == 1:
            index[i] = close[0]
        else:
            print 'Multiple objects found near x,y = (%f,%f)'%(x1[i],y1[i])
            amin = numpy.argmin((x2[close] - x1[i])**2 + (y2[close] - y1[i])**2)
            print 'Found minimum at ',close[amin],', where (x,y) = (%f,%f)'%(
                    x2[close[amin]], y2[close[amin]])
            index[i] = close[amin]
    return index

 
def find_fs_index(used_data, fs_data):
    """Find the index in the fs_data records corresponding to each star in used_data.
    """
    used_x = used_data['X_IMAGE']
    used_y = used_data['Y_IMAGE']
    fs_x = fs_data['x']
    fs_y = fs_data['y']
    return find_index(used_x, used_y, fs_x, fs_y)

def find_used_index(fs_data, used_data):
    """Find the index in the used_data records corresponding to each star in fs_data.
    """
    fs_x = fs_data['x']
    fs_y = fs_data['y']
    used_x = used_data['X_IMAGE']
    used_y = used_data['Y_IMAGE']
    return find_index(fs_x, fs_y, used_x, used_y)


def measure_shapes(xlist, ylist, file_name, wcs):
    """Given x,y positions, an image file, and the wcs, measure shapes and sizes.

    We use the HSM module from GalSim to do this.

    Returns e1, e2, size, flag.
    """
    import galsim
    import numpy
    import astropy.io.fits as pyfits

    #im = galsim.fits.read(file_name)
    #bp_im = galsim.fits.read(file_name, hdu=2)
    #wt_im = galsim.fits.read(file_name, hdu=3)
    with pyfits.open(file_name) as f:
        # Some DES images have bad cards.  Fix them with verify('fix') before sending to GalSim.
        f[1].verify('fix')
        f[2].verify('fix')
        f[3].verify('fix')
        print 'after verify, f = ',f
        print 'f[1] = ',f[1]
        print 'f[2] = ',f[2]
        print 'f[3] = ',f[3]
        im = galsim.fits.read(hdu_list=f, hdu=1, compression='rice')
        print 'im = ',im
        bp_im = galsim.fits.read(hdu_list=f, hdu=2, compression='rice')
        wt_im = galsim.fits.read(hdu_list=f, hdu=3, compression='rice')

    # The badpix image is offset by 32768 from the true value.  Subtract it off.
    if numpy.any(bp_im.array > 32767):
        bp_im -= 32768
    # Also, convert to int16, since it isn't by default weirdly.  I think this is
    # an error in astropy's RICE algorith, since fitsio converts it correctly to uint16.
    bp_im = galsim.ImageS(bp_im)

    # Also, it seems that the weight image has negative values where it should be 0.
    # Make them 0.
    wt_im.array[wt_im.array < 0] = 0.
    print 'file_name = ',file_name

    # Read the background image as well.
    bkg_file_name = file_name[:-8] + '_bkg.fits.fz'
    print 'bkg_file_name = ',bkg_file_name
    #bkg_im = galsim.fits.read(bkg_file_name)
    with pyfits.open(bkg_file_name) as f:
        f[1].verify('fix')
        bkg_im = galsim.fits.read(hdu_list=f, hdu=1, compression='rice')
    im -= bkg_im # Subtract off the sky background.

    stamp_size = 48

    n_psf = len(xlist)
    e1_list = [ 999. ] * n_psf
    e2_list = [ 999. ] * n_psf
    s_list = [ 999. ] * n_psf
    flag_list = [ 0 ] * n_psf
    print 'len(xlist) = ',len(xlist)

    for i in range(n_psf):
        x = xlist[i]
        y = ylist[i]
        print 'Measure shape for star at ',x,y
        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2, 
                           int(y)-stamp_size/2, int(y)+stamp_size/2)

        try:
            subim = im[b]
            subbp = bp_im[b]
            subwt = wt_im[b]
            #print 'subim = ',subim.array
            #print 'subwt = ',subwt.array
            #print 'subbp = ',subbp.array
            #shape_data = subim.FindAdaptiveMom(weight=subwt, badpix=subbp, strict=False)
            shape_data = subim.FindAdaptiveMom(weight=subwt, strict=False)
        except Exception as e:
            print 'Caught ',e
            print ' *** Bad measurement (caught exception).  Mask this one.'
            flag_list[i] = MEAS_BAD_MEASUREMENT
            continue

        #print 'shape_data = ',shape_data
        #print 'image_bounds = ',shape_data.image_bounds
        print 'shape = ',shape_data.observed_shape
        print 'sigma = ',shape_data.moments_sigma
        #print 'amp = ',shape_data.moments_amp
        #print 'centroid = ',shape_data.moments_centroid
        #print 'rho4 = ',shape_data.moments_rho4
        #print 'niter = ',shape_data.moments_n_iter

        if shape_data.moments_status != 0:
            print 'status = ',shape_data.moments_status
            print ' *** Bad measurement.  Mask this one.'
            flag_list[i] = MEAS_BAD_MEASUREMENT
            continue

        dx = shape_data.moments_centroid.x - x
        dy = shape_data.moments_centroid.y - y
        #print 'dcentroid = ',dx,dy
        if dx**2 + dy**2 > 0.1**2:
            print ' *** Centroid shifted by ',dx,dy,'.  Mask this one.'
            flag_list[i] = MEAS_CENTROID_SHIFT
            continue

        e1 = shape_data.observed_shape.e1
        e2 = shape_data.observed_shape.e2
        s = shape_data.moments_sigma
        # Note: this is (det M)^1/4, not ((Ixx+Iyy)/2)^1/2.
        # For reference, the latter is size * (1-e^2)^-1/4
        # So, not all that different, especially for stars with e ~= 0.

        # Account for the WCS:
        jac = wcs.jacobian(galsim.PositionD(x,y))
        # ( Iuu  Iuv ) = ( dudx  dudy ) ( Ixx  Ixy ) ( dudx  dvdx )
        # ( Iuv  Ivv )   ( dvdx  dvdy ) ( Ixy  Iyy ) ( dudy  dvdy )
        M = numpy.matrix( [[ 1+e1, e2 ], [ e2, 1-e1 ]] )
        #print 'M = ',M
        #print 'det(M) = ',numpy.linalg.det(M)
        M2 = jac.getMatrix() * M * jac.getMatrix().T
        #print 'M2 = ',M2
        #print 'det(M2) = ',numpy.linalg.det(M2)
        e1 = (M2[0,0] - M2[1,1]) / (M2[0,0] + M2[1,1])
        e2 = (2. * M2[0,1]) / (M2[0,0] + M2[1,1])
        #print 's = ',s
        s *= abs(numpy.linalg.det(jac.getMatrix()))**0.5
        #print 's -> ',s
        # Now convert back to a more normal shear definition, rather than distortion.
        shear = galsim.Shear(e1=e1,e2=e2)
        e1_list[i] = shear.g1
        e2_list[i] = shear.g2
        s_list[i] = s

    return e1_list,e2_list,s_list,flag_list

def measure_psfex_shapes(xlist, ylist, psfex_file_name, file_name):
    """Given x,y positions, a psfex solution file, and the wcs, measure shapes and sizes
    of the PSF model.

    We use the HSM module from GalSim to do this.

    Returns e1, e2, size, flag.
    """
    import galsim
    import galsim.des
    import numpy
    print 'Read in PSFEx file: ',psfex_file_name

    n_psf = len(xlist)
    e1_list = [ 999. ] * n_psf
    e2_list = [ 999. ] * n_psf
    s_list = [ 999. ] * n_psf
    flag_list = [ 0 ] * n_psf

    try:
        psfex = galsim.des.DES_PSFEx(psfex_file_name, file_name)
    except Exception as e:
        print 'Caught ',e
        flag_list = [ PSFEX_FAILURE ] * n_psf
        return e1_list,e2_list,s_list,flag_list

    stamp_size = 64
    pixel_scale = 0.2

    im = galsim.Image(stamp_size, stamp_size, scale=pixel_scale)

    for i in range(n_psf):
        x = xlist[i]
        y = ylist[i]
        print 'Measure PSFEx model shape at ',x,y
        image_pos = galsim.PositionD(x,y)
        psf = psfex.getPSF(image_pos)
        im = psf.drawImage(image=im, method='no_pixel')
        #print 'im = ',im.array

        try:
            shape_data = im.FindAdaptiveMom(strict=False)
        except:
            print ' *** Bad measurement (caught exception).  Mask this one.'
            flag_list[i] = PSFEX_BAD_MEASUREMENT
            continue

        if shape_data.moments_status != 0:
            print 'status = ',shape_data.moments_status
            print ' *** Bad measurement.  Mask this one.'
            flag_list[i] = PSFEX_BAD_MEASUREMENT
            continue

        dx = shape_data.moments_centroid.x - im.trueCenter().x
        dy = shape_data.moments_centroid.y - im.trueCenter().y
        #print 'centroid = ',shape_data.moments_centroid
        #print 'trueCenter = ',im.trueCenter()
        #print 'dcentroid = ',dx,dy
        if dx**2 + dy**2 > 0.1**2:
            print ' *** Centroid shifted by ',dx,dy,'.  Mask this one.'
            flag_list[i] = PSFEX_CENTROID_SHIFT
            continue

        g1 = shape_data.observed_shape.g1
        g2 = shape_data.observed_shape.g2
        s = shape_data.moments_sigma * pixel_scale

        #print 'g1,g2,s = ',g1,g2,s

        e1_list[i] = g1
        e2_list[i] = g2
        s_list[i] = s

    return e1_list,e2_list,s_list,flag_list


def apply_wcs(wcs, g1, g2, s):
    import numpy

    scale = 2./(1.+g1*g1+g2*g2)
    e1 = g1 * scale
    e2 = g2 * scale
    I = numpy.matrix( [[ 1 + e1, e2 ], [ e2, 1 - e1 ]] ) * s*s

    J = wcs.getMatrix()
    I = J * I * J.transpose()

    e1 = (I[0,0] - I[1,1]) / (I[0,0] + I[1,1])
    e2 = (2.*I[0,1]) / (I[0,0] + I[1,1])
    s = numpy.sqrt( (I[0,0] + I[1,1]) / 2.)

    esq = e1*e1 + e2*e2
    scale = (1.-numpy.sqrt(1.-esq))/esq
    g1 = e1 * scale
    g2 = e2 * scale

    return g1, g2, s

def measure_psfex_shapes_erin(xlist, ylist, psfex_file_name, file_name):
    """Given x,y positions, a psfex solution file, and the wcs, measure shapes and sizes
    of the PSF model.

    We use the HSM module from GalSim to do this.

    Also, this uses Erin's psfex module to render the images rather than the GalSim module.

    Returns e1, e2, size, flag.
    """
    import galsim
    import numpy
    import psfex
    print 'Read in PSFEx file: ',psfex_file_name

    n_psf = len(xlist)
    e1_list = [ 999. ] * n_psf
    e2_list = [ 999. ] * n_psf
    s_list = [ 999. ] * n_psf
    flag_list = [ 0 ] * n_psf

    try:
        psf = psfex.PSFEx(psfex_file_name)
    except Exception as e:
        print 'Caught ',e
        flag_list = [ PSFEX_FAILURE ] * n_psf
        return e1_list,e2_list,s_list,flag_list

    if psf._psfex is None:
        # Erin doesn't throw an exception for errors.
        # The _psfex attribute just ends up as None, so check for that.
        print 'psf._psfex is None'
        flag_list = [ PSFEX_FAILURE ] * n_psf
        return e1_list,e2_list,s_list,flag_list

    wcs = galsim.FitsWCS(file_name)

    for i in range(n_psf):
        x = xlist[i]
        y = ylist[i]
        print 'Measure PSFEx model shape at ',x,y,' with Erin\'s code.'

        # Note that this code renders the image in the original coordinate system,
        # rather than in RA/Dec oriented coordinates.  So we'll need to correct for that.
        im_ar = psf.get_rec(y,x)
        local_wcs = wcs.jacobian(galsim.PositionD(x,y))
        #print 'local wcs = ',local_wcs
        #print 'pixel scale = ',numpy.sqrt(local_wcs.pixelArea())
        #print 'psf center = ',psf.get_center(y,x)
        pixel_scale= numpy.sqrt(local_wcs.pixelArea())
        im = galsim.Image(array=im_ar, scale=pixel_scale)

        try:
            shape_data = im.FindAdaptiveMom(strict=False)
        except:
            print ' *** Bad measurement (caught exception).  Mask this one.'
            flag_list[i] = PSFEX_BAD_MEASUREMENT
            continue

        if shape_data.moments_status != 0:
            print 'status = ',shape_data.moments_status
            print ' *** Bad measurement.  Mask this one.'
            flag_list[i] = PSFEX_BAD_MEASUREMENT
            continue

        cen = psf.get_center(y,x)
        true_center = galsim.PositionD( cen[1]+1, cen[0]+1 )
        dx = shape_data.moments_centroid.x - true_center.x
        dy = shape_data.moments_centroid.y - true_center.y
        #print 'centroid = ',shape_data.moments_centroid
        #print 'trueCenter = ',true_center
        #print 'dcentroid = ',dx,dy
        if dx**2 + dy**2 > 0.5**2:
            print ' *** Centroid shifted by ',dx,dy,'.  Mask this one.'
            flag_list[i] = PSFEX_CENTROID_SHIFT
            continue
        #print 'shape = ',shape_data.observed_shape
        #print 'sigma = ',shape_data.moments_sigma * pixel_scale
        g1,g2,s = apply_wcs(local_wcs, shape_data.observed_shape.g1, shape_data.observed_shape.g2,
                            shape_data.moments_sigma)
        #print 'after wcs: ',g1,g2,s

        e1_list[i] = g1
        e2_list[i] = g2
        s_list[i] = s

    return e1_list,e2_list,s_list,flag_list

def read_blacklists(tag):
    """Read the psfex blacklist file and the other blacklists.

    Returns a dict indexed by the tuple (expnum, ccdnum) with the bitmask value.
    """
    import astropy.io.fits as pyfits
    import numpy
    d = {}  # The dict will be indexed by (expnum, ccdnum)
    print 'reading blacklists'

    if False:
        # First read Eli's astrometry flags
        # cf. https://github.com/esheldon/deswl/blob/master/deswl/desmeds/genfiles.py#L498
        eli_file = '/astro/u/astrodat/data/DES/EXTRA/astrorerun/sva1_astrom_run1.0.1_stats_flagged_sheldon.fit'
        with pyfits.open(eli_file) as pyf:
            data = pyf[1].data
            for expnum, ccdnum, flag in zip(data['EXPNUM'],data['CCDNUM'],data['ASTROM_FLAG']):
                key = (int(expnum), int(ccdnum))
                d[key] = int(flag)
        print 'after astrom, len(d) = ',len(d)

    # Then Alex and Steve's blacklists
    # cf. https://github.com/esheldon/deswl/blob/master/deswl/desmeds/genfiles.py#L588)
    ghost_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/ghost-scatter-y1-uniq.txt'
    streak_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/streak-y1-uniq.txt'
    noise_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/noise-y1-uniq.txt'
    with open(ghost_file) as f:
        for line in f:
            expnum, ccdnum = line.split()
            key = (int(expnum), int(ccdnum))
            if key in d:
                d[key] |= (1 << 10)
            else:
                d[key] = (1 << 10)
    with open(noise_file) as f:
        for line in f:
            expnum, ccdnum = line.split()
            key = (int(expnum), int(ccdnum))
            if key in d:
                d[key] |= (1 << 11)
            else:
                d[key] = (1 << 11)
    with open(streak_file) as f:
        for line in f:
            expnum, ccdnum = line.split()
            key = (int(expnum), int(ccdnum))
            if key in d:
                d[key] |= (1 << 13)
            else:
                d[key] = (1 << 13)
    print 'after ghost, streak, len(d) = ',len(d)

    # And finally the PSFEx blacklist file.
    psfex_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psfex'
    if tag:
        psfex_file += '-' + tag
    psfex_file += '.txt'
    with open(psfex_file) as f:
        for line in f:
            run, exp, ccdnum, flag = line.split()
            expnum = exp[6:]
            key = (int(expnum), int(ccdnum))
            flag = int(flag)
            if key in d:
                d[key] |= (flag << 15)
            else:
                d[key] = (flag << 15)
    print 'after psfex, len(d) = ',len(d)

    return d


def main():
    import os
    import glob
    import galsim
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

    flag_dict = read_blacklists(args.tag)

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    # Directory to put output files.
    cat_dir = os.path.join(work,'psf_cats')
    if not os.path.exists(cat_dir):
        os.makedirs(cat_dir)

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        # The input directory from the main DESDM reduction location.
        input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))

        # Get the file names in that directory.
        files = glob.glob('%s/%s'%(input_dir,args.exp_match))

        # Setup the columns for the output catalog:
        ccdnum_col = []
        x_col = []
        y_col = []
        ra_col = []
        dec_col = []
        mag_col = []
        flag_col = []
        e1_col = []
        e2_col = []
        size_col = []
        psfex_e1_col = []
        psfex_e2_col = []
        psfex_size_col = []
        erin_e1_col = []
        erin_e2_col = []
        erin_size_col = []
        desdm_e1_col = []
        desdm_e2_col = []
        desdm_size_col = []

        for file_name in files:
            print '\nProcessing ', file_name

            try:
                desdm_dir, root, ccdnum = parse_file_name(file_name)
            except:
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue
            print '   root, ccdnum = ',root,ccdnum
            print '   desdm_dir = ',desdm_dir

            key = (expnum, ccdnum)
            if key in flag_dict:
                black_flag = flag_dict[key]
                print '   blacklist flag = ',black_flag
                if black_flag & (113 << 15):
                    print '   Catastrophic flag.  Skipping this file.'
                    if args.single_ccd:
                        break
                    continue
            else:
                black_flag = 0

            # Read the star data.  From both findstars and the PSFEx used file.
            fs_data = read_findstars(exp_dir, root)
            if fs_data is None:
                print '   No _findstars.fits file found'
                if args.single_ccd:
                    break
                continue
            n_tot = len(fs_data)
            n_fs = fs_data['star_flag'].sum()
            print '   n_tot = ',n_tot
            print '   n_fs = ',n_fs
            mask = fs_data['star_flag'] == 1

            used_data = read_used(exp_dir, root)
            if used_data is None:
                print '   No .used.fits file found'
                continue
            n_used = len(used_data)
            print '   n_used = ',n_used
            if n_used == 0:
                print '   No stars were used.'
                continue

            tot_xmin = fs_data['x'].min()
            tot_xmax = fs_data['x'].max()
            tot_ymin = fs_data['y'].min()
            tot_ymax = fs_data['y'].max()
            tot_area = (tot_xmax-tot_xmin)*(tot_ymax-tot_ymin)
            print '   bounds from sextractor = ',tot_xmin,tot_xmax,tot_ymin,tot_ymax
            print '   area = ',tot_area

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
            print '   fraction used = ',float(used_area) / tot_area

            # Figure out which fs objects go with which used objects.
            fs_index = find_fs_index(used_data, fs_data)
            used_index = find_used_index(fs_data[mask], used_data)
            print '   fs_index = ',fs_index
            print '   used_index = ',used_index

            # Check: This should be the same as the used bounds
            alt_used_xmin = fs_data['x'][fs_index].min()
            alt_used_xmax = fs_data['x'][fs_index].max()
            alt_used_ymin = fs_data['y'][fs_index].min()
            alt_used_ymax = fs_data['y'][fs_index].max()
            print '   bounds from findstars[fs_index] = ',
            print alt_used_xmin,alt_used_xmax,alt_used_ymin,alt_used_ymax
 
            # Get the magnitude range for each catalog.
            tot_magmin = fs_data['mag'].min()
            tot_magmax = fs_data['mag'].max()
            print '   magnitude range of full catalog = ',tot_magmin,tot_magmax
            fs_magmin = fs_data['mag'][mask].min()
            fs_magmax = fs_data['mag'][mask].max()
            print '   magnitude range of fs stars = ',fs_magmin,fs_magmax
            used_magmin = fs_data['mag'][fs_index].min()
            used_magmax = fs_data['mag'][fs_index].max()
            print '   magnitude range of used stars = ',used_magmin,used_magmax

            try:
                # Get the wcs from the image file
                wcs = get_wcs(file_name)

                # Measure the shpes and sizes of the stars used by PSFEx.
                x = fs_data['x'][mask]
                y = fs_data['y'][mask]
                mag = fs_data['mag'][mask]
                e1, e2, size, meas_flag = measure_shapes(x, y, file_name, wcs)

                # Measure the model shapes, sizes.
                psfex_file_name = os.path.join(exp_dir, root + '_psfcat.psf')
                psfex_e1, psfex_e2, psfex_size, psfex_flag = measure_psfex_shapes(
                        x, y, psfex_file_name, file_name)
                erin_e1, erin_e2, erin_size, erin_flag = measure_psfex_shapes_erin(
                        x, y, psfex_file_name, file_name)

                # Measure the desdm model shapes, sizes.
                desdm_file_name = os.path.join(desdm_dir, root + '_psfcat.psf')
                desdm_e1, desdm_e2, desdm_size, desdm_flag = measure_psfex_shapes(
                        x, y, desdm_file_name, file_name)
                desdm_flag = [ f * DESDM_FLAG_FACTOR for f in desdm_flag ]
            except Exception as e:
                print 'Catastrophic error trying to measure the shapes:'
                print e
                print 'Skip this file'
                continue

            # Put all the flags together:
            flag = [ m | p | e | d for m,p,e,d in zip(meas_flag,psfex_flag,erin_flag,desdm_flag) ]
            print 'meas_flag = ',meas_flag
            print 'psfex_flag = ',psfex_flag
            print 'erin_flag = ',erin_flag
            print 'desdm_flag = ',desdm_flag
            print 'flag = ',flag

            # Add in flags for bad indices
            bad_index = numpy.where(used_index < 0)[0]
            print 'bad_index = ',bad_index
            for i in bad_index:
                flag[i] |= NOT_USED
            print 'flag => ',flag

            # If the ccd is blacklisted, everything gets the blacklist flag
            if black_flag:
                print 'black_flag = ',black_flag
                print 'type(black_flag) = ',type(black_flag)
                print 'type(flag[0]) = ',type(flag[0])
                print 'type(flag[0] | black_flag) = ',type(flag[0] | black_flag)
                black_flag *= BLACK_FLAG_FACTOR
                print 'black_flag => ',black_flag
                flag = [ f | black_flag for f in flag ]
                print 'flag => ',flag

            # Compute ra,dec from the wcs:
            coord = [ wcs.toWorld(galsim.PositionD(xx,yy)) for xx,yy in zip(x,y) ]
            ra = [ c.ra / galsim.degrees for c in coord ]
            dec = [ c.dec / galsim.degrees for c in coord ]

            # Extend the column arrays with this chip's data.
            ccdnum_col.extend([ccdnum] * n_fs)
            x_col.extend(x)
            y_col.extend(y)
            ra_col.extend(ra)
            dec_col.extend(dec)
            mag_col.extend(mag)
            flag_col.extend(flag)
            e1_col.extend(e1)
            e2_col.extend(e2)
            size_col.extend(size)
            psfex_e1_col.extend(psfex_e1)
            psfex_e2_col.extend(psfex_e2)
            psfex_size_col.extend(psfex_size)
            erin_e1_col.extend(erin_e1)
            erin_e2_col.extend(erin_e2)
            erin_size_col.extend(erin_size)
            desdm_e1_col.extend(desdm_e1)
            desdm_e2_col.extend(desdm_e2)
            desdm_size_col.extend(desdm_size)
            assert len(ccdnum_col) == len(x_col)
            assert len(ccdnum_col) == len(y_col)
            assert len(ccdnum_col) == len(ra_col)
            assert len(ccdnum_col) == len(dec_col)
            assert len(ccdnum_col) == len(mag_col)
            assert len(ccdnum_col) == len(flag_col)
            assert len(ccdnum_col) == len(e1_col)
            assert len(ccdnum_col) == len(e2_col)
            assert len(ccdnum_col) == len(size_col)
            assert len(ccdnum_col) == len(psfex_e1_col)
            assert len(ccdnum_col) == len(psfex_e2_col)
            assert len(ccdnum_col) == len(psfex_size_col)
            assert len(ccdnum_col) == len(erin_e1_col)
            assert len(ccdnum_col) == len(erin_e2_col)
            assert len(ccdnum_col) == len(erin_size_col)
            assert len(ccdnum_col) == len(desdm_e1_col)
            assert len(ccdnum_col) == len(desdm_e2_col)
            assert len(ccdnum_col) == len(desdm_size_col)

            if args.single_ccd:
                break

        cols = pyfits.ColDefs([
            pyfits.Column(name='ccdnum', format='I', array=ccdnum_col),
            pyfits.Column(name='x', format='E', array=x_col),
            pyfits.Column(name='y', format='E', array=y_col),
            pyfits.Column(name='ra', format='E', array=ra_col),
            pyfits.Column(name='dec', format='E', array=dec_col),
            pyfits.Column(name='mag', format='E', array=mag_col),
            pyfits.Column(name='flag', format='J', array=flag_col),
            pyfits.Column(name='e1', format='E', array=e1_col),
            pyfits.Column(name='e2', format='E', array=e2_col),
            pyfits.Column(name='size', format='E', array=size_col),
            pyfits.Column(name='psfex_e1', format='E', array=psfex_e1_col),
            pyfits.Column(name='psfex_e2', format='E', array=psfex_e2_col),
            pyfits.Column(name='psfex_size', format='E', array=psfex_size_col),
            pyfits.Column(name='erin_e1', format='E', array=erin_e1_col),
            pyfits.Column(name='erin_e2', format='E', array=erin_e2_col),
            pyfits.Column(name='erin_size', format='E', array=erin_size_col),
            pyfits.Column(name='desdm_e1', format='E', array=desdm_e1_col),
            pyfits.Column(name='desdm_e2', format='E', array=desdm_e2_col),
            pyfits.Column(name='desdm_size', format='E', array=desdm_size_col),
            ])

        # Depending on the version of pyfits, one of these should work:
        try:
            tbhdu = pyfits.BinTableHDU.from_columns(cols)
        except:
            tbhdu = pyfits.new_table(cols)
        exp_root = root.rsplit('_',1)[0]
        print 'exp_root = ',exp_root
        cat_file = os.path.join(cat_dir, exp_root + "_psf.fits")
        tbhdu.writeto(cat_file, clobber=True)
        print 'wrote cat_file = ',cat_file

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
