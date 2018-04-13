#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import os
import numpy as np
from read_psf_cats import read_data

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description='Run rho stats on a set of exposures')

    # Drectory arguments
    parser.add_argument('--work', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--file', default='',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Options
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--use_psfex', default=False, action='store_const', const=True,
                        help='Use PSFEx rather than Piff model')
    parser.add_argument('--bands', default='grizY', type=str,
                        help='Limit to the given bands')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')

    args = parser.parse_args()
    return args


def measure_rho(data, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff'):
    """Compute the rho statistics
    """
    import treecorr

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    T = data['obs_T']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    print('mean e = ',np.mean(e1),np.mean(e2))
    print('mean T = ',np.mean(T))
    print('mean de = ',np.mean(de1),np.mean(de2))
    print('mean dT = ',np.mean(T-p_T))
    print('mean dT/T = ',np.mean(dt))

    if use_xy:
        x = data['fov_x']
        y = data['fov_y']
        print('x = ',x)
        print('y = ',y)

        ecat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=e1, g2=e2)
        decat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec',
                                 k=dt, g1=dt*e1, g2=dt*e2)
    else:
        ra = data['ra']
        dec = data['dec']
        print('ra = ',ra)
        print('dec = ',dec)

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

    bin_config = dict(
        sep_units = 'arcmin',
        bin_slop = 0.1,

        min_sep = 0.5,
        max_sep = max_sep,
        bin_size = 0.2,

        #min_sep = 2.5,
        #max_sep = 250,
        #nbins = 20,
    )

    results = []
    for (cat1, cat2) in [ (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(dtcat.name, dtcat.name))

        rho = treecorr.KKCorrelation(bin_config, verbose=2)
        rho.process(dtcat)
        results.append(rho)

    return results


def measure_cross_rho(tile_data, max_sep, tags=None, prefix='piff'):
    """Compute the rho statistics
    """
    import treecorr

    ntilings = len(tile_data)
    print('len(tile_data) = ',ntilings)

    de1 = [ d['obs_e1']-d[prefix+'_e1'] for d in tile_data ]
    de2 = [ d['obs_e2']-d[prefix+'_e2'] for d in tile_data ]
    dt = [ (d['obs_T']-d[prefix+'_T'])/d['obs_T'] for d in tile_data ]
    for k in range(len(tile_data)):
        print('k = ',k)
        print('mean de = ',np.mean(de1[k]),np.mean(de2[k]))
        print('mean dt = ',np.mean(dt[k]))

    ecats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg',
                               g1=d['obs_e1'], g2=d['obs_e2'])
              for d in tile_data ]
    for cat in ecats: cat.name = 'ecat'
    decats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg',
                                g1=de1[k], g2=de2[k]) for k,d in enumerate(tile_data) ]
    for cat in decats: cat.name = 'decat'
    dtcats = [ treecorr.Catalog(ra=d['ra'], dec=d['dec'], ra_units='deg', dec_units='deg',
                                k=dt[k], g1=d['obs_e1']*dt[k], g2=d['obs_e2']*dt[k])
               for k,d in enumerate(tile_data) ]
    for cat in dtcats: cat.name = 'dtcat'
    if tags is not None:
        for catlist in [ ecats, decats, dtcats]:
            for cat, tag in zip(catlist, tags):
                cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = 'arcmin',
        bin_slop = 0.1,

        min_sep = 0.5,
        max_sep = max_sep,
        bin_size = 0.2,

        #min_sep = 2.5,
        #max_sep = 250,
        #nbins = 20,
    )

    results = []
    for (catlist1, catlist2) in [ (decats, decats),
                                  (ecats, decats),
                                  (dtcats, dtcats),
                                  (decats, dtcats),
                                  (ecats, dtcats) ]:

        catnames1 = [ cat.name for cat in catlist1 ]
        catnames2 = [ cat.name for cat in catlist2 ]
        print('Doing correlation of %s vs %s'%(catnames1, catnames2))
        rho = treecorr.GGCorrelation(bin_config, verbose=2)
        # Avoid all auto correlations:
        for i in range(ntilings):
            for j in range(ntilings):
                if i == j: continue
                if catlist1 is catlist2 and i > j: continue
                print('names: ',catlist1[i].name,catlist2[j].name)
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
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)


def band_combinations(bands, single=True, combo=True):

    if single:
        use_bands = [ [b] for b in bands ]
    else:
        use_bands = []

    if combo:
        if 'r' in bands and 'i' in bands:
            use_bands.append(['r', 'i'])
        if 'r' in bands and 'i' in bands and 'z' in bands:
            use_bands.append(['r', 'i', 'z'])
        if 'g' in bands and 'r' in bands and 'i' in bands and 'z' in bands:
            use_bands.append(['g', 'r', 'i', 'z'])

    print('use_bands = ',use_bands)
    print('tags = ',[ ''.join(band) for band in use_bands ])
    return use_bands


def do_canonical_stats(data, bands, tilings, work, prefix='piff', name='all', alt_tt=False):
    print('Start CANONICAL: ',prefix,name)
    # Measure the canonical rho stats using all pairs:
    use_bands = band_combinations(bands)
    for band in use_bands:
        print('band ',band)
        mask = np.in1d(data['band'],band)
        print('sum(mask) = ',np.sum(mask))
        print('len(data[mask]) = ',len(data[mask]))
        tag = ''.join(band)
        stats = measure_rho(data[mask], max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)

def do_cross_tiling_stats(data, bands, tilings, work, prefix='piff', name='cross'):
    print('Start CROSS_TILING: ',prefix,name)
    # Measure the rho stats using only cross-correlations between tiles.
    use_bands = band_combinations(bands)
    for band in use_bands:
        tile_data = []
        for til in tilings:
            print('til = ',til)
            mask = np.in1d(data['band'],band) & (data['tiling'] == til)
            print('sum(mask) = ',np.sum(mask))
            print('len(data[mask]) = ',len(data[mask]))
            tile_data.append(data[mask])
        tag = ''.join(band)
        tags = [ tag + ":" + str(til) for til in tilings ]
        stats = measure_cross_rho(tile_data, max_sep=300, tags=tags, prefix=prefix)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_cross_band_stats(data, bands, tilings, work, prefix='piff', name='crossband'):
    print('Start CROSS_BAND: ',prefix,name)
    # Measure the rho stats cross-correlating the different bands.
    use_bands = band_combinations(bands, single=False)

    for band in use_bands:
        print('cross bands ',band)
        band_data = []
        for b in band:
            mask = data['band'] == b
            band_data.append(data[mask])
        stats = measure_cross_rho(band_data, max_sep=300, tags=band, prefix=prefix)
        tag = ''.join(band)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_odd_even_stats(data, bands, tilings, work, prefix='piff', name='oddeven'):
    print('Start ODD_EVEN: ',prefix,name)
    # Measure the rho stats using only cross-correlations between odd vs even tilings.
    use_bands = band_combinations(bands)

    for band in use_bands:
        print('odd/even ',band)
        odd = np.in1d(data['band'], band) & (data['tiling'] % 2 == 1)
        even = np.in1d(data['band'], band) & (data['tiling'] % 2 == 0)
        cats = [ data[odd], data[even] ]
        tag = ''.join(band)
        tags = [ tag + ":odd", tag + ":even" ]
        stats = measure_cross_rho(cats, max_sep=300, tags=tags, prefix=prefix)
        stat_file = os.path.join(work, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)


def do_fov_stats(data, bands, tilings, work, prefix='piff', name='fov'):
    print('Start FOV: ',prefix,name)
    # Measure the rho stats using the field-of-view positions.
    use_bands = band_combinations(bands)
    for band in use_bands:
        print('band ',band)
        mask = np.in1d(data['band'],band)
        print('sum(mask) = ',np.sum(mask))
        print('len(data[mask]) = ',len(data[mask]))
        tag = ''.join(band)
        stats = measure_rho(data[mask], max_sep=300, tag=tag, prefix=prefix, use_xy=True)
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
    print("Writing data to ",file_name)
    fitsio.write(file_name, data, clobber=True)

def main():

    args = parse_args()
    print('args = ',args)

    # Make the work directory if it does not exist yet.
    work = os.path.expanduser(args.work)
    print('work dir = ',work)
    try:
        if not os.path.isdir(work):
            os.makedirs(work)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    if args.use_psfex:
        prefix='psfex'
    else:
        prefix='piff'

    if args.file != '':
        print('Read file ',args.file)
        with open(args.file) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        print('File included %d exposures'%len(exps))
    else:
        exps = args.exps
        print('Explicit listing of %d exposures'%len(exps))
    exps = sorted(exps)

    keys = ['ra', 'dec', 'x', 'y', 'obs_e1', 'obs_e2', 'obs_T',
            prefix+'_e1', prefix+'_e2', prefix+'_T']

    data, bands, tilings = read_data(exps, work, keys, limit_bands=args.bands, prefix=prefix,
                                     use_reserved=args.use_reserved, frac=args.frac)

    print('all bands = ',bands)
    #print('all tilings = ',tilings)

    out_file_name = os.path.join(work, "psf_%s.fits"%args.tag)
    write_data(data, out_file_name)

    for band in bands:
        print('n for band %s = '%band, np.sum(data['band'] == band))
    #for til in tilings:
        #print('n for tiling %d = '%til, np.sum(data['tiling'] == til))

    gdata = np.where(data['band'] == 'g')[0]
    rdata = np.where(data['band'] == 'r')[0]
    idata = np.where(data['band'] == 'i')[0]
    zdata = np.where(data['band'] == 'z')[0]
    #odddata = np.where(data['tiling']%2 == 1)[0]
    #evendata = np.where(data['tiling']%2 == 0)[0]

    print('len(gdata) = ',len(gdata))
    print('len(rdata) = ',len(rdata))
    print('len(idata) = ',len(idata))
    print('len(zdata) = ',len(zdata))
    #print('len(odddata) = ',len(odddata))
    #print('len(evendata) = ',len(evendata))

    #bands = ['r', 'i']

    do_canonical_stats(data, bands, tilings, work, alt_tt=False, prefix=prefix)

    #do_cross_tiling_stats(data, bands, tilings, work, prefix=prefix)

    #do_cross_band_stats(data, bands, tilings, work, prefix=prefix)

    #do_odd_even_stats(data, bands, tilings, work, prefix=prefix)

    #do_fov_stats(data, bands, tilings, work, prefix=prefix)


if __name__ == "__main__":
    main()
