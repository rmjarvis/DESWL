#! /usr/bin/env python
# Build a catalog with information about each exposure/ccdnum.
# This includes information from the image header, the blacklists, and the results
# of the rho statistics.

import astropy.io.fits as pyfits
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')
    parser.add_argument('--output', default='exposure_info.fits',
                        help='The name of the output fits file')

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
    psfex_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psfex-y1'
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

    Returns date, time, filter, ccdnum, detpos, telra, teldec, ha, airmass, sky, sigsky, fwhm, wcs
    """
    print 'Start read_image_header'
    print img_file
    import galsim
    #import fitsio
    import pyfits

    if img_file.endswith('fz'):
        hdu = 1
    else:
        hdu = 0

    # fitsio is a bit faster here.  11 sec/exp rather than 12, so not a huge difference, but still.
    with pyfits.open(img_file) as pyf:
    #with fitsio.FITS(img_file) as pyf:
        #print pyf
        #print pyf[hdu]
        h = pyf[hdu].header
        #h = pyf[hdu].read_header()
        #print 'opened'
        # DATE-OBS looks like '2012-12-03T07:38:54.174780', so split on T.
        date = h['DATE-OBS']
        date, time = date.strip().split('T',1)

        # FILTER looks like 'z DECam SDSS c0004 9260.0 1520.0', so split on white space
        filter = h['FILTER']
        filter = filter.split()[0]

        # CCDNUM is 1-62.  DETPOS is a string such as 'S29     '.  Strip off the whitespace.
        ccdnum = h['CCDNUM']
        ccdnum = int(ccdnum)
        detpos = h['DETPOS']
        detpos = detpos.strip()
        #print 'detpos = ',detpos

        # TELRA, TELDEC look like '-62:30:22.320'.  Use GalSim to convert to decimal degrees.
        telra = h['TELRA']
        teldec = h['TELDEC']
        telra = galsim.HMS_Angle(telra) / galsim.degrees
        teldec = galsim.DMS_Angle(teldec) / galsim.degrees
        #print 'telra, teldec = ',telra,teldec

        # HA looks like '03:12:02.70''.  Use GalSim to convert to decimal degrees.
        ha = h['HA']
        ha = galsim.HMS_Angle(ha) / galsim.degrees

        # A few more items to grab from the header:
        airmass = float(h['AIRMASS'])
        sky = float(h['SKYBRITE'])
        sigsky = float(h['SKYSIGMA'])
        fwhm = float(h['FWHM'])

        # Use Galsim to read WCS
        wcs = galsim.FitsWCS(header=h)
        #print 'wcs = ',wcs

    return date, time, filter, ccdnum, detpos, telra, teldec, ha, airmass, sky, sigsky, fwhm, wcs
 
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
    y = toYearFraction(d)
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


def main():
    import os
    import glob
    import galsim

    args = parse_args()

    tbdata = read_tapebump_file('/astro/u/mjarvis/rmjarvis/DESWL/psfex/mask_ccdnum.txt')

    flag_dict = read_blacklists(args.tag)

    work = os.path.expanduser(args.work)
    print 'work dir = ',work

    datadir = '/astro/u/astrodat/data/DES'

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    # Setup the columns we will put in the output catalog:
    expnum_col = []  # The exposure number
    ccdnum_col = []  # The ccd number
    run_col = []     # In which run did DESDM process this?
    exp_col = []     # What is the full text of the exposure name
    root_col = []    # Just the root name
    date_col = []    # The date as a string
    time_col = []    # The time as a string
    year_col = []    # The date as a decimal year
    filter_col = []  # Which filter is this exposure
    detpos_col = []  # The code for this CCD (e.g. S29)
    telra_col = []   # The ra of the telescope pointing (degrees)
    teldec_col = []  # The dec of the telescopt pointing (degrees)
    ha_col = []      # The hour angle (degrees)
    airmass_col = [] # The airmass
    sky_col = []     # The median sky level
    sigsky_col = []  # The mean noise level from the sky
    fwhm_col = []    # An estimate of the seeing
    flag_col = []    # A bitmask flag for the ccd (or possibly the whole exposure)
    corner0_ra_col = []  # The ra,dec of the 4 corners of the chip (degrees)
    corner0_dec_col = []
    corner1_ra_col = []
    corner1_dec_col = []
    corner2_ra_col = []
    corner2_dec_col = []
    corner3_ra_col = []
    corner3_dec_col = []
    bump0_ra_col = []  # The ra,dec of the centers of the 6 tape bumps (degrees)
    bump0_dec_col = []
    bump1_ra_col = []
    bump1_dec_col = []
    bump2_ra_col = []
    bump2_dec_col = []
    bump3_ra_col = []
    bump3_dec_col = []
    bump4_ra_col = []
    bump4_dec_col = []
    bump5_ra_col = []
    bump5_dec_col = []

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
            except Exception as e:
                print '   Caught ',e
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue
            print '   root, ccdnum = ',root,ccdnum

            try:
                (date, time, filter, ccdnum2, detpos, telra, teldec, ha, 
                    airmass, sky, sigsky, fwhm, wcs) = read_image_header(file_name)
                print '   date, time = ',date,time
                print '   filter, ccdnum, detpos = ', filter,ccdnum,detpos
                print '   telra, teldec, ha = ', telra, teldec, ha
                print '   airmass, sky, sigsky, fwhm = ', airmass, sky, sigsky, fwhm
                if ccdnum != ccdnum2:
                    raise ValueError("CCDNUM from FITS header doesn't match ccdnum from file name.")
            except Exception as e:
                print '   Caught ',e
                print '   Error reading fits header.  Skipping this file:',file_name
                continue

            year = convert_to_year(date, time)
            print '   year = ',year

            # Check if we have a blacklist flag for this chip
            key = (expnum, ccdnum)
            if key in flag_dict:
                flag = flag_dict[key]
                print '   flag = ',flag
            else:
                flag = 0
                print '   not flagged'

            # These are useful to calculate as "special" positions for testing.
            corners = [ wcs.toWorld(galsim.PositionD(i,j)) 
                        for i,j in [ (0,0), (2048,0), (0,4096), (2048,4096) ] ]
            #print '   corners = ',corners

            # Also figure out the location of each tape bump.  (More "special" positions)
            print '   nbumps = ',len(tbdata[ccdnum])
            assert len(tbdata[ccdnum]) == 6
            bumps = [ wcs.toWorld(bump_center(bump)) for bump in tbdata[ccdnum] ]
            #print '   bumps = ',bumps

            expnum_col.append(expnum)
            ccdnum_col.append(ccdnum)
            run_col.append(run)
            exp_col.append(exp)
            root_col.append(root)
            date_col.append(date)
            time_col.append(time)
            year_col.append(year)
            filter_col.append(filter)
            detpos_col.append(detpos)
            telra_col.append(telra)
            teldec_col.append(teldec)
            ha_col.append(ha)
            airmass_col.append(airmass)
            sky_col.append(sky)
            sigsky_col.append(sigsky)
            fwhm_col.append(fwhm)
            flag_col.append(flag)
            corner0_ra_col.append(corners[0].ra / galsim.degrees)
            corner0_dec_col.append(corners[0].dec / galsim.degrees)
            corner1_ra_col.append(corners[1].ra / galsim.degrees)
            corner1_dec_col.append(corners[1].dec / galsim.degrees)
            corner2_ra_col.append(corners[2].ra / galsim.degrees)
            corner2_dec_col.append(corners[2].dec / galsim.degrees)
            corner3_ra_col.append(corners[3].ra / galsim.degrees)
            corner3_dec_col.append(corners[3].dec / galsim.degrees)
            bump0_ra_col.append(bumps[0].ra / galsim.degrees)
            bump0_dec_col.append(bumps[0].dec / galsim.degrees)
            bump1_ra_col.append(bumps[1].ra / galsim.degrees)
            bump1_dec_col.append(bumps[1].dec / galsim.degrees)
            bump2_ra_col.append(bumps[2].ra / galsim.degrees)
            bump2_dec_col.append(bumps[2].dec / galsim.degrees)
            bump3_ra_col.append(bumps[3].ra / galsim.degrees)
            bump3_dec_col.append(bumps[3].dec / galsim.degrees)
            bump4_ra_col.append(bumps[4].ra / galsim.degrees)
            bump4_dec_col.append(bumps[4].dec / galsim.degrees)
            bump5_ra_col.append(bumps[5].ra / galsim.degrees)
            bump5_dec_col.append(bumps[5].dec / galsim.degrees)

            if args.single_ccd:
                break

    print '\nFinished processing all exposures'

    # Check the length required for string columns:
    run_len = max([ len(s) for s in run_col ])
    exp_len = max([ len(s) for s in exp_col ])
    root_len = max([ len(s) for s in root_col ])
    date_len = max([ len(s) for s in date_col ])
    time_len = max([ len(s) for s in time_col ])
    filter_len = max([ len(s) for s in filter_col ])
    detpos_len = max([ len(s) for s in detpos_col ])

    cols = pyfits.ColDefs([
        pyfits.Column(name='expnum', format='J', array=expnum_col),
        pyfits.Column(name='ccdnum', format='I', array=ccdnum_col),
        pyfits.Column(name='run', format='A%d'%run_len, array=run_col),
        pyfits.Column(name='exp', format='A%d'%exp_len, array=exp_col),
        pyfits.Column(name='root', format='A%d'%root_len, array=root_col),
        pyfits.Column(name='date', format='A%d'%date_len, array=date_col),
        pyfits.Column(name='time', format='A%d'%time_len, array=time_col),
        pyfits.Column(name='year', format='E', array=year_col),
        pyfits.Column(name='filter', format='A%d'%filter_len, array=filter_col),
        pyfits.Column(name='detpos', format='A%d'%detpos_len, array=detpos_col),
        pyfits.Column(name='telra', format='E', unit='deg', array=telra_col),
        pyfits.Column(name='teldec', format='E', unit='deg', array=teldec_col),
        pyfits.Column(name='ha', format='E', unit='deg', array=ha_col),
        pyfits.Column(name='airmass', format='E', array=airmass_col),
        pyfits.Column(name='sky', format='E', array=sky_col),
        pyfits.Column(name='sigsky', format='E', array=sigsky_col),
        pyfits.Column(name='fwhm', format='E', array=fwhm_col),
        pyfits.Column(name='flag', format='J', array=flag_col),
        pyfits.Column(name='corner0_ra', format='E', unit='deg', array=corner0_ra_col),
        pyfits.Column(name='corner0_dec', format='E', unit='deg', array=corner0_dec_col),
        pyfits.Column(name='corner1_ra', format='E', unit='deg', array=corner1_ra_col),
        pyfits.Column(name='corner1_dec', format='E', unit='deg', array=corner1_dec_col),
        pyfits.Column(name='corner2_ra', format='E', unit='deg', array=corner2_ra_col),
        pyfits.Column(name='corner2_dec', format='E', unit='deg', array=corner2_dec_col),
        pyfits.Column(name='corner3_ra', format='E', unit='deg', array=corner3_ra_col),
        pyfits.Column(name='corner3_dec', format='E', unit='deg', array=corner3_dec_col),
        pyfits.Column(name='bump0_ra', format='E', unit='deg', array=bump0_ra_col),
        pyfits.Column(name='bump0_dec', format='E', unit='deg', array=bump0_dec_col),
        pyfits.Column(name='bump1_ra', format='E', unit='deg', array=bump1_ra_col),
        pyfits.Column(name='bump1_dec', format='E', unit='deg', array=bump1_dec_col),
        pyfits.Column(name='bump2_ra', format='E', unit='deg', array=bump2_ra_col),
        pyfits.Column(name='bump2_dec', format='E', unit='deg', array=bump2_dec_col),
        pyfits.Column(name='bump3_ra', format='E', unit='deg', array=bump3_ra_col),
        pyfits.Column(name='bump3_dec', format='E', unit='deg', array=bump3_dec_col),
        pyfits.Column(name='bump4_ra', format='E', unit='deg', array=bump4_ra_col),
        pyfits.Column(name='bump4_dec', format='E', unit='deg', array=bump4_dec_col),
        pyfits.Column(name='bump5_ra', format='E', unit='deg', array=bump5_ra_col),
        pyfits.Column(name='bump5_dec', format='E', unit='deg', array=bump5_dec_col),
        ])

    # Depending on the version of pyfits, one of these should work:
    try:
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
    except:
        tbhdu = pyfits.new_table(cols)
    tbhdu.writeto(args.output, clobber=True)


if __name__ == "__main__":
    main()
