#! /usr/bin/env python
# Program to computer rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='',
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

def bump_center(bump):
    """A quick helper function that returns the center of a tape bump
    """
    import galsim
    return galsim.PositionD( (bump[1]+bump[3])/2., (bump[0]+bump[2])/2. )

def read_blacklists(tag):
    """Read the psfex blacklist file and the other blacklists.

    Returns a dict indexed by the tuple (expnum, ccdnum) with the bitmask value.
    """
    import pyfits
    import numpy
    d = {}  # The dict will be indexed by (expnum, ccdnum)
    print 'reading blacklists'

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
    ghost_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/ghost-scatter-sv-uniq.txt'
    streak_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/streak-sv-uniq.txt'
    with open(ghost_file) as f:
        for line in f:
            expnum, ccdnum = line.split()
            key = (int(expnum), int(ccdnum))
            if key in d:
                d[key] |= (1 << 10)
            else:
                d[key] = (1 << 10)
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
    psfex_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psfex-sv'
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

def read_image_header(img_file):
    """Read some information from the image header.

    Returns date, time, filter, ccdnum, detpos, telra, teldec, ha, airmass, wcs
    """
    import pyfits
    import galsim

    if img_file.endswith('fz'):
        hdu = 1
    else:
        hdu = 0

    with pyfits.open(img_file) as pyf:
        # DATE-OBS looks like '2012-12-03T07:38:54.174780', so split on T.
        date = pyf[hdu].header['DATE-OBS']
        date, time = date.strip().split('T',1)

        # FILTER looks like 'z DECam SDSS c0004 9260.0 1520.0', so split on white space
        filter = pyf[hdu].header['FILTER']
        filter = filter.split()[0]

        # CCDNUM is 1-62.  DETPOS is a string such as 'S29     '.  Strip off the whitespace.
        ccdnum = pyf[hdu].header['CCDNUM']
        ccdnum = int(ccdnum)
        detpos = pyf[hdu].header['DETPOS']
        detpos = detpos.strip()

        # TELRA, TELDEC look like '-62:30:22.320'.  Use GalSim to convert to decimal degrees.
        telra = pyf[hdu].header['TELRA']
        teldec = pyf[hdu].header['TELDEC']
        telra = galsim.HMS_Angle(telra) / galsim.degrees
        teldec = galsim.DMS_Angle(teldec) / galsim.degrees

        # HA looks like '03:12:02.70''.  Use GalSim to convert to decimal degrees.
        ha = pyf[hdu].header['HA']
        ha = galsim.HMS_Angle(ha) / galsim.degrees

        # AIRMASS should already be a float I think, but just make sure.
        airmass = pyf[hdu].header['AIRMASS']
        airmass = float(airmass)

        # Use Galsim to read WCS
        wcs = galsim.FitsWCS(header=pyf[hdu].header)

    return date, time, filter, ccdnum, detpos, telra, teldec, ha, airmass, wcs
 
def convert_to_year(date, time):
    """Given string representations of the date and time, convert to a decimal year.

    e.g. 2012.89142
    """
    # cf. http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    from datetime import datetime as dt

    def toYearFraction(date):
        """Input is a datetime object.  Returns the year with a decimal fraction.
        """
        import time

        def s(date): # returns seconds from canonical epoch to the given date
            return time.mktime(date.timetuple())

        year = date.year
        startOfThisYear = dt(year=year, month=1, day=1)
        startOfNextYear = dt(year=year+1, month=1, day=1)
        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = yearElapsed/yearDuration

        return date.year + fraction

    # Date looks like 2012-12-03
    year, month, day = date.split('-')
    year = int(year)
    month = int(month)
    day = int(day)

    # Time looks like 07:38:54.174780
    hour, min, sec = time.split(':')
    hour = int(hour)
    min = int(min)
    sec = float(sec)

    d = dt(year, month, day, hour, min, int(sec), int((sec-int(sec)) * 1.e6))
    print 'd = ',d
    y = toYearFraction(d)
    print 'y = ',y
    return y

def parse_file_name(file_name):
    """Parse the PSFEx file name to get the root name and the chip number
    """
    import os

    base_file = os.path.split(file_name)[1]
    if os.path.splitext(base_file)[1] == '.fz':
        base_file=os.path.splitext(base_file)[0]
    if os.path.splitext(base_file)[1] != '.fits':
        raise ValueError("Invalid file name "+file)
    root = os.path.splitext(base_file)[0]

    ccdnum = int(root.split('_')[-1])
    return root, ccdnum


def read_used(exp_dir, root):
    """Read in the .used.fits file that PSFEx generates with the list of stars that actually
    got used in making the PSFEx file.
    """
    import os
    import astropy.io.fits as pyfits
    import copy

    file_name = os.path.join(exp_dir, root + '_psfcat.used.fits')
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


def measure_shapes(xlist, ylist, file_name, wcs):
    """Given x,y positions, an image file, and the wcs, measure shapes and sizes.

    We use the HSM module from GalSim to do this.

    Returns e1, e2, size, mask.
    """
    import galsim
    import numpy

    im = galsim.fits.read(file_name)
    bp_im = galsim.fits.read(file_name, hdu=2)
    wt_im = galsim.fits.read(file_name, hdu=3)
    wt_im.array[wt_im.array < 0] = 0.
    print 'file_name = ',file_name
    bkg_file_name = file_name[:-8] + '_bkg.fits.fz'
    print 'bgk_file_name = ',bkg_file_name
    bkg_im = galsim.fits.read(bkg_file_name)
    im -= bkg_im # Subtract off the sky background.

    stamp_size = 32

    e1_list = numpy.zeros_like(xlist)
    e2_list = numpy.zeros_like(xlist)
    s_list = numpy.zeros_like(xlist)
    mask = numpy.ones_like(xlist, dtype=bool)

    for i in range(len(xlist)):
        x = xlist[i]
        y = ylist[i]
        print 'Measure shape for star at ',x,y
        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2, 
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        subim = im[b]
        subbp = bp_im[b]
        subwt = wt_im[b]
        shape_data = subim.FindAdaptiveMom(weight=subwt, badpix=subbp, strict=False)
        #print 'shape_data = ',shape_data
        #print 'image_bounds = ',shape_data.image_bounds
        #print 'shape = ',shape_data.observed_shape
        #print 'sigma = ',shape_data.moments_sigma
        #print 'amp = ',shape_data.moments_amp
        #print 'centroid = ',shape_data.moments_centroid
        #print 'rho4 = ',shape_data.moments_rho4
        #print 'niter = ',shape_data.moments_n_iter

        if shape_data.moments_status != 0:
            print 'status = ',shape_data.moments_status
            print ' *** Bad measurement.  Mask this one.'
            mask[i] = 0
            continue

        dx = shape_data.moments_centroid.x - x
        dy = shape_data.moments_centroid.y - y
        #print 'dcentroid = ',dx,dy
        if dx**2 + dy**2 > 0.1**2:
            print ' *** Centroid shifted by ',dx,dy,'.  Mask this one.'
            mask[i] = 0
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
        e1_list[i] = e1
        e2_list[i] = e2
        s_list[i] = s

    return e1_list,e2_list,s_list,mask


def measure_psfex_shapes(xlist, ylist, psfex_file_name, file_name):
    """Given x,y positions, a psfex solution file, and the wcs, measure shapes and sizes
    of the PSF model.

    We use the HSM module from GalSim to do this.

    Returns e1, e2, size, mask.
    """
    import galsim
    import galsim.des
    import numpy

    psfex = galsim.des.DES_PSFEx(psfex_file_name, file_name)

    stamp_size = 32
    pixel_scale = 0.1

    e1 = numpy.zeros_like(xlist)
    e2 = numpy.zeros_like(xlist)
    size = numpy.zeros_like(xlist)
    mask = numpy.ones_like(xlist, dtype=bool)
    im = galsim.Image(stamp_size, stamp_size, scale=pixel_scale)

    for i in range(len(xlist)):
        x = xlist[i]
        y = ylist[i]
        print 'Measure PSFEx model shape at ',x,y
        image_pos = galsim.PositionD(x,y)
        psf = psfex.getPSF(image_pos)
        im = psf.drawImage(image=im, method='no_pixel')
        #print 'im = ',im.array

        shape_data = im.FindAdaptiveMom(strict=False)

        if shape_data.moments_status != 0:
            print 'status = ',shape_data.moments_status
            print ' *** Bad measurement.  Mask this one.'
            mask[i] = 0
            continue

        dx = shape_data.moments_centroid.x - im.trueCenter().x
        dy = shape_data.moments_centroid.y - im.trueCenter().y
        #print 'centroid = ',shape_data.moments_centroid
        #print 'trueCenter = ',im.trueCenter()
        #print 'dcentroid = ',dx,dy
        if dx**2 + dy**2 > 0.1**2:
            print ' *** Centroid shifted by ',dx,dy,'.  Mask this one.'
            mask[i] = 0
            continue

        #print 'shape = ',shape_data.observed_shape
        e1[i] = shape_data.observed_shape.e1
        e2[i] = shape_data.observed_shape.e2
        size[i] = shape_data.moments_sigma * pixel_scale

    return e1,e2,size,mask

def measure_rho(e1,e2,s,m_e1,m_e2,m_s,mask):
    """Compute the rho statistics
    """
    import numpy
    import treecorr

    print 'e1 = ',e1
    print 'm_e1 = ',m_e1
    print 'de1 = ',e1-m_e1
    print 'mean de1 = ',numpy.mean(e1[mask]-m_e1[mask])
    print 'e2 = ',e2
    print 'm_e2 = ',m_e2
    print 'de2 = ',e2-m_e2
    print 'mean de2 = ',numpy.mean(e2[mask]-m_e2[mask])
    print 's = ',s
    print 'm_s = ',m_s
    print 'ds = ',s-m_s
    print 'mean ds = ',numpy.mean(s[mask]-m_s[mask])


def main():
    import os
    import glob
    import galsim
    import json

    args = parse_args()

    tbdata = read_tapebump_file('/astro/u/mjarvis/rmjarvis/DESWL/psfex/mask_ccdnum.txt')

    flag_dict = read_blacklists(args.tag)

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
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        # The input directory from the main DESDM reduction location.
        input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))

        # Get the file names in that directory.
        files = glob.glob('%s/%s'%(input_dir,args.exp_match))

        for file_name in files:
            print '\nProcessing ', file_name

            # Start by getting some basic information about the exposure / chip
            # to put in the output file
            try:
                root, ccdnum = parse_file_name(file_name)
            except:
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue
            print '   root, ccdnum = ',root,ccdnum

            date, time, filter, ccdnum2, detpos, telra, teldec, ha, airmass, wcs = \
                read_image_header(file_name)
            print '   date, time = ',date,time
            print '   filter, ccdnum, detpos = ',filter,ccdnum,detpos
            print '   telra, teldec, ha, airmass = ',telra, teldec, ha, airmass
            if ccdnum != ccdnum2:
                raise ValueError("CCDNUM from FITS header doesn't match ccdnum from file name.")

            year = convert_to_year(date, time)
            print '   year = ',year

            # Check if we have a blacklist flag for this chip
            key = (expnum, ccdnum)
            if key in flag_dict:
                flag = flag_dict[key]
                print '   flag = ',flag
                print '   type(flag) = ',type(flag)
            else:
                flag = 0
                print '   not flagged'

            # Read the star data.  From both findstars and the PSFEx used file.
            fs_data = read_findstars(exp_dir, root)
            n_tot = len(fs_data)
            n_fs = fs_data['star_flag'].sum()
            print '   n_tot = ',n_tot
            print '   n_fs = ',n_fs
            mask = fs_data['star_flag'] == 1

            used_data = read_used(exp_dir, root)
            n_used = len(used_data)
            print '   n_used = ',n_used
            if len(used_data) == 0:
                print '   No stars were used.'
                continue

            # Figure out which fs objects go with which used objects.
            index = find_fs_index(used_data, fs_data)
            print '   index = ',index

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

            # Check: This should be the same as the used bounds
            alt_used_xmin = fs_data['x'][index].min()
            alt_used_xmax = fs_data['x'][index].max()
            alt_used_ymin = fs_data['y'][index].min()
            alt_used_ymax = fs_data['y'][index].max()
            print '   bounds from findstars[index] = ',
            print alt_used_xmin,alt_used_xmax,alt_used_ymin,alt_used_ymax
 
            # These are useful to calculate as "special" positions for testing.
            corners = [ wcs.toWorld(galsim.PositionD(i,j)) 
                        for i,j in [ (0,0), (2048,0), (0,4096), (2048,4096) ] ]
            #print '   corners = ',corners

            # Also figure out the location of each tape bump.  (More "special" positions)
            print '   nbumps = ',len(tbdata[ccdnum])
            assert len(tbdata[ccdnum]) == 6
            bumps = [ wcs.toWorld(bump_center(bump)) for bump in tbdata[ccdnum] ]
            #print '   bumps = ',bumps

            # Measure the shpes and sizes of the stars used by PSFEx.
            x = used_data['X_IMAGE']
            y = used_data['Y_IMAGE']
            e1, e2, size, mask = measure_shapes(x, y, file_name, wcs)

            # Measure the model shapes, sizes.
            psfex_file_name = os.path.join(exp_dir, root + '_psfcat.psf')
            m_e1, m_e2, m_size, m_mask = measure_psfex_shapes(x, y, psfex_file_name, file_name)

            # Compute the correlation functions
            #rho1, rho2, rho3 = measure_rho(e1,e2,size,m_e1,m_e2,m_size,mask&m_mask)
            
            # Write out the interesting stats for this ccd into a file, which we can 
            # then pull all together into a single FITS catalog later.
            stat_file = os.path.join(exp_dir, root + ".json")
            stats = [ 
                expnum, ccdnum, run, exp, root, 
                date, time, year, filter, detpos,
                telra, teldec, ha, airmass,
                flag,
                tot_xmin, tot_xmax, tot_ymin, tot_ymax,
                fs_xmin, fs_xmax, fs_ymin, fs_ymax,
                used_xmin, used_xmax, used_ymin, used_ymax,
                tot_area, fs_area, used_area,
                n_tot, n_fs, n_used,
                corners[0].ra / galsim.degrees, corners[0].dec / galsim.degrees,
                corners[1].ra / galsim.degrees, corners[1].dec / galsim.degrees,
                corners[2].ra / galsim.degrees, corners[2].dec / galsim.degrees,
                corners[3].ra / galsim.degrees, corners[3].dec / galsim.degrees,
                bumps[0].ra / galsim.degrees, bumps[0].dec / galsim.degrees,
                bumps[1].ra / galsim.degrees, bumps[1].dec / galsim.degrees,
                bumps[2].ra / galsim.degrees, bumps[2].dec / galsim.degrees,
                bumps[3].ra / galsim.degrees, bumps[3].dec / galsim.degrees,
                bumps[4].ra / galsim.degrees, bumps[4].dec / galsim.degrees,
                bumps[5].ra / galsim.degrees, bumps[5].dec / galsim.degrees,
                ]
            with open(stat_file,'w') as f:
                json.dump(stats, f)

            if args.single_ccd:
                break

    print '\nFinished processing all exposures'


if __name__ == "__main__":
    main()
