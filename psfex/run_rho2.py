#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

 
import os
import numpy
from toFocal import toFocal

def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits.fz',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')
    parser.add_argument('--max_tiling', default=10,
                        help='maximum tiling to use')
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')

    # Options
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')
    parser.add_argument('--oldkeys', default=False, action='store_const', const=True,
                        help='Use old psfex_* keys in the cats files')

    args = parser.parse_args()
    return args


def read_data(args, work, limit_filters=None, subtract_mean=True, reserved=False):
    import astropy.io.fits as pyfits

    datadir = '/astro/u/astrodat/data/DES'

    RESERVED = 64
    BAD_CCDS = [2, 31, 61]

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    expinfo_file = '/astro/u/mjarvis/work/exposure_info_' + args.tag + '.fits'
    #expinfo_file = '/astro/u/mjarvis/work/exposure_info_y1a1-v01.fits'
    #expinfo_file = 'exposure_info_y1spte-v02.fits'
    print 'reading exposure_info file: ',expinfo_file
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    keys = ['ra', 'dec', 'x', 'y', 'e1', 'e2', 'size', 'psf_e1', 'psf_e2', 'psf_size']
    all_data = { key : [] for key in keys }
    all_keys = keys

    all_data['exp'] = []
    all_data['ccd'] = []
    all_keys = all_keys + ['exp', 'ccd' ]

    if subtract_mean:
        all_data['alt_e1'] = []
        all_data['alt_e2'] = []
        all_data['alt_size'] = []
        all_keys = all_keys + ['alt_e1', 'alt_e2', 'alt_size']
    if 'x' in keys:
        all_data['fov_x'] = []
        all_data['fov_y'] = []
        all_keys = all_keys + ['fov_x', 'fov_y']

    if args.oldkeys:
        inkeys =  [ k.replace('psf_','psfex_') for k in keys ]
    else:
        inkeys = keys

    all_filters = []  # This keeps track of the filter for each record
    all_tilings = []  # This keeps track of the tiling for each record
    filters = set()   # This is the set of all filters being used
    tilings = set()   # This is the set of all tilings being used

    ntot = 0
    nused = 0
    nreserved = 0
    ngood = 0

    cat_dir = os.path.join(work,'psf_cats')

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        if expnum not in expinfo['expnum']:
            print 'expnum is not in expinfo!'
            print 'expinfo[expnum] = ',expinfo['expnum']
            print 'Could not find information about this expnum.  Skipping ',run,exp
            continue
        k = numpy.nonzero(expinfo['expnum'] == expnum)[0][0]
        print 'k = ',k
        filter = expinfo['filter'][k]
        print 'filter[k] = ',filter
        if (limit_filters is not None) and (filter not in limit_filters):
            print 'Not doing this filter.'
            continue

        tiling = int(expinfo['tiling'][k])
        print 'tiling[k] = ',tiling

        if tiling == 0:
            # This shouldn't happen, but it did for a few exposures.  Just skip them, since this
            # might indicate some kind of problem.
            print 'tiling == 0.  Skip this exposure.'
            continue

        if tiling > args.max_tiling:
            print 'tiling is > %d.  Skip this exposure.'%args.max_tiling
            continue

        cat_file = os.path.join(cat_dir, exp + "_exppsf.fits")
        if not os.path.exists(cat_file):
            cat_file = os.path.join(cat_dir, exp + "_psf.fits")
        print 'cat_file = ',cat_file
        try:
            with pyfits.open(cat_file) as pyf:
                data = pyf[1].data
        except:
            print 'Unable to open cat_file %s.  Skipping this file.'%cat_file
            continue

        #ccdnums = numpy.unique(data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        print 'n = ',len(data)
        print 'nused = ',numpy.sum((data['flag'] & 1) != 0)
        print 'nreserved = ',numpy.sum((data['flag'] & 64) != 0)
        print 'ngood = ',numpy.sum(data['flag'] == 0)
        ntot += len(data)
        nused += numpy.sum((data['flag'] & 1) != 0)
        nreserved += numpy.sum((data['flag'] & 64) != 0)
        ngood += numpy.sum(data['flag'] == 0)

        bad_ccd = numpy.in1d(data['ccdnum'], BAD_CCDS)
        if args.use_reserved:
            mask = (data['flag'] == RESERVED) & ~bad_ccd
        else:
            mask = (data['flag'] == 0) & ~bad_ccd

        ngood = numpy.sum(mask)
        assert ngood == len(data[mask])
        if ngood == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        for key, inkey in zip(keys, inkeys):
            all_data[key].append(data[inkey][mask])

        all_data['exp'].append([expnum] * ngood)
        all_data['ccd'].append(data['ccdnum'][mask])

        if subtract_mean:
            e1 = data['e1'][mask]
            e2 = data['e2'][mask]
            s = data['size'][mask]
            p_e1 = data['psf_e1'][mask]
            p_e2 = data['psf_e2'][mask]
            p_s = data['psf_size'][mask]
            de1 = numpy.mean(e1-p_e1)
            de2 = numpy.mean(e2-p_e2)
            # Really want <(s^2 - p_s^2)/s^2> => 0  after subtracting ds
            # <1 - p_s^2/(s-ds)^2> = 0
            # <1 - p_s^2/s^2 (1 + 2ds/s + 3ds^2/s^2 + ...)  > = 0
            # 1 - <p_s^2/s^2> - 2ds<p_s^2/s^3> - 3ds^2<p_s^2/s^4> = 0
            a1 = numpy.mean(p_s**2/s**2)
            a2 = numpy.mean(p_s**2/s**3)
            a3 = numpy.mean(p_s**2/s**4)
            ds = (1. - a1) / (2.*a2)
            # Iterate once to refine
            ds = (1. - a1 - 3.*ds**2*a3) / (2.*a2)
            print 'de = ',de1,de2, 'mean e = ',numpy.mean(e1),numpy.mean(e2),
            print ' -> ',numpy.mean(e1-de1), numpy.mean(e2-de2)
            print 'ds = ',ds, 'mean s = ',numpy.mean(s),' -> ',numpy.mean(s-ds)
            print 'mean dt = ',numpy.mean( (s**2 - p_s**2) / s**2 ),
            print ' -> ',numpy.mean( ((s-ds)**2 - p_s**2) / (s-ds)**2 )
            all_data['alt_e1'].append(e1 - de1)
            all_data['alt_e2'].append(e2 - de2)
            all_data['alt_size'].append(s - ds)

        if 'x' in keys:
            # Convert to focal position.
            x,y = toFocal(data['ccdnum'][mask], data['x'][mask], data['y'][mask])
            # This comes back in units of mm.  Convert to arcsec.
            # 1 pixel = 15e-3 mm = 0.263 arcsec
            x *= 0.263/15e-3
            y *= 0.263/15e-3
            all_data['fov_x'].append(x)
            all_data['fov_y'].append(y)

        all_filters.extend( ([filter] * ngood) )
        all_tilings.extend( ([tiling] * ngood) )
        filters.add(filter)
        tilings.add(tiling)
    print '\nFinished processing all exposures'
    print 'filters = ',filters
    print 'tilings = ',tilings

    # Turn the data into a recarray
    print 'all_data.keys = ',all_data.keys()
    formats = ['f8'] * len(all_keys) + ['a1', 'i2']
    names = all_keys + ['filter', 'tiling']
    data = numpy.recarray(shape = (len(all_filters),),
                          formats = formats, names = names)
    print 'data.dtype = ',data.dtype
    for key in all_keys:
        data[key] = numpy.concatenate(all_data[key])
    data['filter'] = all_filters
    data['tiling'] = all_tilings
    print 'made recarray'

    print 'ntot = ',ntot
    print 'nused = ',nused
    print 'nreserved = ',nreserved
    print 'ngood = ',ngood

    return data, filters, tilings


def measure_rho(data, max_sep, tag=None, prefix='', use_xy=False, alt_tt=False):
    """Compute the rho statistics
    """
    import treecorr

    e1 = data[prefix+'e1']
    e2 = data[prefix+'e2']
    s = data[prefix+'size']
    p_e1 = data['psf_e1']
    p_e2 = data['psf_e2']
    p_s = data['psf_size']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (s**2-p_s**2)/s**2
    print 'mean de = ',numpy.mean(de1),numpy.mean(de2)
    print 'mean dt = ',numpy.mean(dt)

    if use_xy:
        x = data['fov_x']
        y = data['fov_y']
        print 'x = ',x
        print 'y = ',y

        ecat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=e1, g2=e2)
        decat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec',
                                 k=dt, g1=dt*e1, g2=dt*e2)
    else:
        ra = data['ra']
        dec = data['dec']
        print 'ra = ',ra
        print 'dec = ',dec

        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                                 k=dt, g1=dt*e1, g2=dt*e2)
    ecat.name = 'ecat'
    decat.name = 'decat'
    dtcat.name = 'dtcat'
    if tag is not None:
        for cat in [ ecat, decat, dtcat ]:
            cat.name = tag + ":"  + cat.name

    min_sep = 0.5
    bin_size = 0.2
    bin_slop = 0.1

    results = []
    for (cat1, cat2) in [ (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:
        print 'Doing correlation of %s vs %s'%(cat1.name, cat2.name)

        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        results.append(rho)

    if alt_tt:
        print 'Doing alt correlation of %s vs %s'%(dtcat.name, dtcat.name)

        rho = treecorr.KKCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
        rho.process(dtcat)
        results.append(rho)

    return results


def measure_cross_rho(tile_data, max_sep, tags=None, prefix=''):
    """Compute the rho statistics
    """
    import treecorr

    ntilings = len(tile_data)
    print 'len(tile_data) = ',ntilings

    de1 = [ d[prefix+'e1']-d['psf_e1'] for d in tile_data ]
    de2 = [ d[prefix+'e2']-d['psf_e2'] for d in tile_data ]
    dt = [ (d[prefix+'size']**2-d['psf_size']**2)/d['size']**2 for d in tile_data ]
    for k in range(len(tile_data)):
        print 'k = ',k
        print 'mean de = ',numpy.mean(de1[k]),numpy.mean(de2[k])
        print 'mean dt = ',numpy.mean(dt[k])

    ecats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg', 
                               g1=d['e1'], g2=d['e2']) 
              for d in tile_data ]
    for cat in ecats: cat.name = 'ecat'
    decats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg', 
                                g1=de1[k], g2=de2[k]) for k,d in enumerate(tile_data) ]
    for cat in decats: cat.name = 'decat'
    dtcats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg', 
                                k=dt[k], g1=d['e1']*dt[k], g2=d['e2']*dt[k])
               for k,d in enumerate(tile_data) ]
    for cat in dtcats: cat.name = 'dtcat'
    if tags is not None:
        for catlist in [ ecats, decats, dtcats]:
            for cat, tag in zip(catlist, tags):
                cat.name = tag + ":"  + cat.name

    min_sep = 0.5
    bin_size = 0.2
    bin_slop = 0.1

    results = []
    for (catlist1, catlist2) in [ (decats, decats),
                                  (ecats, decats),
                                  (dtcats, dtcats),
                                  (decats, dtcats),
                                  (ecats, dtcats) ]:

        catnames1 = [ cat.name for cat in catlist1 ]
        catnames2 = [ cat.name for cat in catlist2 ]
        print 'Doing correlation of %s vs %s'%(catnames1, catnames2)
        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
        # Avoid all auto correlations:
        for i in range(ntilings):
            for j in range(ntilings):
                if i == j: continue
                if catlist1 is catlist2 and i > j: continue
                print 'names: ',catlist1[i].name,catlist2[j].name
                rho.process_cross(catlist1[i], catlist2[j])
        varg1 = treecorr.calculateVarG(catlist1)
        varg2 = treecorr.calculateVarG(catlist2)
        rho.finalize(varg1, varg2)
        results.append(rho)

    return results

def write_stats(stat_file, rho1, rho2, rho3, rho4, rho5, corr_tt=None):
    import json

    stats = [
        rho1.meanlogr.tolist(),
        rho1.xip.tolist(),
        rho1.xip_im.tolist(),
        rho1.xim.tolist(),
        rho1.xim_im.tolist(),
        rho1.varxi.tolist(),
        rho2.xip.tolist(),
        rho2.xip_im.tolist(),
        rho2.xim.tolist(),
        rho2.xim_im.tolist(),
        rho2.varxi.tolist(),
        rho3.xip.tolist(),
        rho3.xip_im.tolist(),
        rho3.xim.tolist(),
        rho3.xim_im.tolist(),
        rho3.varxi.tolist(),
        rho4.xip.tolist(),
        rho4.xip_im.tolist(),
        rho4.xim.tolist(),
        rho4.xim_im.tolist(),
        rho4.varxi.tolist(),
        rho5.xip.tolist(),
        rho5.xip_im.tolist(),
        rho5.xim.tolist(),
        rho5.xim_im.tolist(),
        rho5.varxi.tolist(),
    ]
    if corr_tt is not None:
        stats.extend([
            corr_tt.xi.tolist(),
            corr_tt.varxi.tolist()
        ])
    #print 'stats = ',stats
    print 'stat_file = ',stat_file
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print 'Done writing ',stat_file


def filter_combinations(filters, single=True, combo=True):

    if single:
        use_filters = [ [f] for f in filters ]
    else:
        use_filters = []

    if combo:
        if 'r' in filters and 'i' in filters:
            use_filters.append(['r', 'i'])
        if 'r' in filters and 'i' in filters and 'z' in filters:
            use_filters.append(['r', 'i', 'z'])
        if 'g' in filters and 'r' in filters and 'i' in filters and 'z' in filters:
            use_filters.append(['g', 'r', 'i', 'z'])

    print 'use_filters = ',use_filters
    print 'tags = ',[ ''.join(filt) for filt in use_filters ]
    return use_filters


def do_canonical_stats(data, filters, tilings, work, prefix='', name='all', alt_tt=False):
    print 'Start CANONICAL: ',prefix,name
    # Measure the canonical rho stats using all pairs:
    use_filters = filter_combinations(filters)
    for filt in use_filters:
        print 'filter ',filt
        mask = numpy.in1d(data['filter'],filt)
        print 'sum(mask) = ',numpy.sum(mask)
        print 'len(data[mask]) = ',len(data[mask])
        tag = ''.join(filt)
        stats = measure_rho(data[mask], max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)

def do_cross_tiling_stats(data, filters, tilings, work, prefix='', name='cross'):
    print 'Start CROSS_TILING: ',prefix,name
    # Measure the rho stats using only cross-correlations between tiles.
    use_filters = filter_combinations(filters)
    for filt in use_filters:
        tile_data = []
        for til in tilings:
            print 'til = ',til
            mask = numpy.in1d(data['filter'],filt) & (data['tiling'] == til)
            print 'sum(mask) = ',numpy.sum(mask)
            print 'len(data[mask]) = ',len(data[mask])
            tile_data.append(data[mask])
        tag = ''.join(filt)
        tags = [ tag + ":" + str(til) for til in tilings ]
        stats = measure_cross_rho(tile_data, max_sep=300, tags=tags, prefix=prefix)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_cross_band_stats(data, filters, tilings, work, prefix='', name='crossband'):
    print 'Start CROSS_BAND: ',prefix,name
    # Measure the rho stats cross-correlating the different bands.
    use_filters = filter_combinations(filters, single=False)

    for filt in use_filters:
        print 'cross filters ',filt
        filt_data = []
        for f in filt:
            mask = data['filter'] == f
            filt_data.append(data[mask])
        stats = measure_cross_rho(filt_data, max_sep=300, tags=filt, prefix=prefix)
        tag = ''.join(filt)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_odd_even_stats(data, filters, tilings, work, prefix='', name='oddeven'):
    print 'Start ODD_EVEN: ',prefix,name
    # Measure the rho stats using only cross-correlations between odd vs even tilings.
    use_filters = filter_combinations(filters)

    for filt in use_filters:
        print 'odd/even ',filt
        odd = numpy.in1d(data['filter'], filt) & (data['tiling'] % 2 == 1)
        even = numpy.in1d(data['filter'], filt) & (data['tiling'] % 2 == 0)
        cats = [ data[odd], data[even] ]
        tag = ''.join(filt)
        tags = [ tag + ":odd", tag + ":even" ]
        stats = measure_cross_rho(cats, max_sep=300, tags=tags, prefix=prefix)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_fov_stats(data, filters, tilings, work, prefix='', name='fov'):
    print 'Start FOV: ',prefix,name
    # Measure the rho stats using the field-of-view positions.
    use_filters = filter_combinations(filters)
    for filt in use_filters:
        print 'filter ',filt
        mask = numpy.in1d(data['filter'],filt)
        print 'sum(mask) = ',numpy.sum(mask)
        print 'len(data[mask]) = ',len(data[mask])
        tag = ''.join(filt)
        stats = measure_rho(data[mask], max_sep=300, tag=tag,
                                                   prefix=prefix, use_xy=True)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def set_args(**kwargs):
    # Used when running from a shell python to setup args.
    class Namespace: 
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            for key in kwargs:
                setattr(self, key, kwargs[key])
        def __repr__(self):
            s = 'Namespace('
            lst = [ '%s=%r'%(k,v) for k,v in self.kwargs.items() ]
            s += ','.join(lst)
            s += ')'
            return s

    args =  Namespace(exp_match='*_[0-9][0-9].fits*',
                      exps='',
                      file='y1spte',
                      max_tiling=10,
                      runs='',
                      single_ccd=False,
                      oldkeys=False,
                      tag='v1',
                      work='~/work/psfex_rerun/v1')
    for key in kwargs:
        setattr(args, key, kwargs[key])
    return args

def write_data(data, file_name):
    import fitsio
    print "Writing data to ",file_name
    fitsio.write(file_name, data, clobber=True)

def main():

    args = parse_args()
    print 'args = ',args

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

    filters = None
    #filters = ['r']
    #filters = ['r', 'i']
    #filters = ['r', 'i', 'z']
    data, filters, tilings = read_data(args, work, limit_filters=filters, subtract_mean=False)

    print 'all filters = ',filters
    print 'all tilings = ',tilings

    out_file_name = os.path.join(work, "psf_%s.fits"%args.tag)
    write_data(data, out_file_name)

    for filt in filters:
        print 'n for filter %s = '%filt, numpy.sum(data['filter'] == filt)
    for til in tilings:
        print 'n for tiling %d = '%til, numpy.sum(data['tiling'] == til)

    gdata = numpy.where(data['filter'] == 'g')[0]
    rdata = numpy.where(data['filter'] == 'r')[0]
    idata = numpy.where(data['filter'] == 'i')[0]
    zdata = numpy.where(data['filter'] == 'z')[0]
    odddata = numpy.where(data['tiling']%2 == 1)[0]
    evendata = numpy.where(data['tiling']%2 == 0)[0]

    print 'len(gdata) = ',len(gdata)
    print 'len(rdata) = ',len(rdata)
    print 'len(idata) = ',len(idata)
    print 'len(zdata) = ',len(zdata)
    print 'len(odddata) = ',len(odddata)
    print 'len(evendata) = ',len(evendata)

    #filters = ['r', 'i']

    do_canonical_stats(data, filters, tilings, work, alt_tt=True)

    #do_cross_tiling_stats(data, filters, tilings, work)

    #do_cross_band_stats(data, filters, tilings, work)

    #do_odd_even_stats(data, filters, tilings, work)

    #do_fov_stats(data, filters, tilings, work)

    # Use subtract_mean=True to do these:
    #do_canonical_stats(data, filters, tilings, work, prefix='alt_', name='alt')

    #do_odd_even_stats(data, filters, tilings, work, prefix='alt_', name='altoddeven')


if __name__ == "__main__":
    main()
