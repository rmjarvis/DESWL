#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import os
import numpy as np
from toFocal import toFocal

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
                        help='list of exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Options
    parser.add_argument('--max_tiling', default=10,
                        help='maximum tiling to use')
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--bands', default='grizY', type=str,
                        help='Limit to the given bands')
    parser.add_argument('--use_psfex', default=False, action='store_const', const=True,
                        help='Use PSFEx rather than Piff model')

    args = parser.parse_args()
    return args


def read_data(args, work, limit_bands=None, prefix='piff'):
    import fitsio

    RESERVED = 64
    BAD_CCDS = [2, 31, 61]

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
    all_data = { key : [] for key in keys }
    all_keys = keys

    all_data['exp'] = []
    all_data['ccd'] = []
    all_keys = all_keys + ['exp', 'ccd' ]

    if 'x' in keys:
        all_data['fov_x'] = []
        all_data['fov_y'] = []
        all_keys = all_keys + ['fov_x', 'fov_y']

    inkeys = keys

    all_bands = []  # This keeps track of the band for each record
    #all_tilings = []  # This keeps track of the tiling for each record
    bands = set()   # This is the set of all bands being used
    #tilings = set()   # This is the set of all tilings being used

    n_reject_mean_dt = 0
    n_reject_mean_de1 = 0
    n_reject_mean_de2 = 0
    n_reject_std_dt = 0
    n_reject_std_de1 = 0
    n_reject_std_de2 = 0
    n_reject_rho2 = 0

    for exp in exps:

        expnum = int(exp)
        try:
            expinfo = fitsio.read(os.path.join(work, exp, 'exp_info_%d.fits'%expnum))
        except Exception as e:
            #print('Caught: ',e)
            #print('Skip this exposure')
            continue

        if expnum not in expinfo['expnum']:
            print('expnum is not in expinfo!')
            print('expinfo[expnum] = ',expinfo['expnum'])
            print('Could not find information about this expnum.  Skipping ',run,exp)
            continue
        i = np.nonzero(expinfo['expnum'] == expnum)[0][0]
        #print('i = ',i)
        band = expinfo['band'][i]
        #print('band[k] = ',band)
        if (limit_bands is not None) and (band not in limit_bands):
            #print('Not doing band = %s.'%band)
            continue

        print('Start work on exp = ',exp)
        #tiling = int(expinfo['tiling'][k])
        #print('tiling[k] = ',tiling)

        #if tiling == 0:
            # This shouldn't happen, but it did for a few exposures.  Just skip them, since this
            # might indicate some kind of problem.
            #print('tiling == 0.  Skip this exposure.')
            #continue

        #if tiling > args.max_tiling:
            #print('tiling is > %d.  Skip this exposure.'%args.max_tiling)
            #continue

        for k in range(len(expinfo)):
            ccdnum = expinfo[k]['ccdnum']
            if expinfo[k]['flag'] != 0:
                #print('Skipping ccd %d because it is blacklisted: '%ccdnum, expinfo[k]['flag'])
                continue
            if ccdnum in BAD_CCDS:
                #print('Skipping ccd %d because it is BAD'%ccdnum)
                continue

            cat_file = os.path.join(work, exp, "psf_cat_%d_%d.fits"%(expnum,ccdnum))
            #print('cat_file = ',cat_file)
            try:
                data = fitsio.read(cat_file)
                flag = data[prefix+'_flag']
            except (OSError, IOError):
                #print('Unable to open cat_file %s.  Skipping this file.'%cat_file)
                continue

            ntot = len(data)
            nused = np.sum((flag & 1) != 0)
            nreserved = np.sum((flag & RESERVED) != 0)
            ngood = np.sum(flag == 0)
            #print('nused = ',nused)
            #print('nreserved = ',nreserved)
            #print('ngood = ',ngood)

            if args.use_reserved:
                mask = (flag == RESERVED) | (flag == RESERVED+1)
            else:
                mask = (flag == 0)
            #print('mask = ',mask)

            T = data['obs_T']
            e1 = data['obs_e1']
            e2 = data['obs_e2']
            dT = (data[prefix + '_T'] - data['obs_T'])
            de1 = (data[prefix + '_e1'] - data['obs_e1'])
            de2 = (data[prefix + '_e2'] - data['obs_e2'])
            used = (flag == 0)
            print(expnum, ccdnum, len(dT), band)
            print('T = ',np.mean(T[used]),np.std(T[used]))
            print('e1 = ',np.mean(e1[used]),np.std(e1[used]))
            print('e2 = ',np.mean(e2[used]),np.std(e2[used]))
            print('dT/T = ',np.mean(dT[used]/T[used]),np.std(dT[used]/T[used]))
            print('de1 = ',np.mean(de1[used]),np.std(de1[used]))
            print('de2 = ',np.mean(de2[used]),np.std(de2[used]))
            rho2 = (e1 - 1j*e2) * (de1 + 1j*de2)
            print('mean rho2 = ',np.mean(rho2))
            if abs(np.mean(dT[used]/T[used])) > 0.01:
                print('mean dT/T = %f on ccd %d.'%(np.mean(dT[used]/T[used]),ccdnum))
                #n_reject_mean_dt += 1
                #continue
            if abs(np.mean(de1[used])) > 0.01:
                print('mean de1 = %f on ccd %d.'%(np.mean(de1[used]),ccdnum))
                #n_reject_mean_de1 += 1
                #continue
            if abs(np.mean(de2[used])) > 0.01:
                print('mean de2 = %f on ccd %d.'%(np.mean(de2[used]),ccdnum))
                #n_reject_mean_de2 += 1
                #continue
            if abs(np.std(dT[used]/T[used])) > 0.1:
                print('std dT/T = %f on ccd %d.'%(np.std(dT[used]/T[used]),ccdnum))
                #n_reject_std_dt += 1
                #continue
            if abs(np.std(de1[used])) > 0.1:
                print('std de1 = %f on ccd %d.'%(np.std(de1[used]),ccdnum))
                #n_reject_std_de1 += 1
                #continue
            if abs(np.std(de2[used])) > 0.1:
                print('std de2 = %f on ccd %d.'%(np.std(de2[used]),ccdnum))
                #n_reject_std_de2 += 1
                #continue
            if abs(np.mean(rho2)) > 5.e-4:
                print('mean rho2 = %f on ccd %d.'%(np.mean(rho2),ccdnum))
                #n_reject_rho2 += 1
                #continue

            good = (abs(dT/T) < 0.1) & (abs(de1) < 0.1) & (abs(de2) < 0.1)
            mask = mask & good

            ngood = np.sum(mask)
            #print('ngood = ',ngood,'/',len(data))
            assert ngood == len(data[mask])
            if ngood == 0:
                print('All objects in ccd %d are flagged.'%ccdnum)
                print('Probably due to astrometry flags. Skip this exposure.')
                continue

            for key, inkey in zip(keys, inkeys):
                all_data[key].append(data[inkey][mask])

            all_data['exp'].append([expnum] * ngood)
            all_data['ccd'].append([ccdnum] * ngood)

            if 'x' in keys:
                # Convert to focal position.
                x,y = toFocal(ccdnum, data['x'][mask], data['y'][mask])
                # This comes back in units of mm.  Convert to arcsec.
                # 1 pixel = 15e-3 mm = 0.263 arcsec
                x *= 0.263/15e-3
                y *= 0.263/15e-3
                all_data['fov_x'].append(x)
                all_data['fov_y'].append(y)

            all_bands.extend( ([band] * ngood) )
            #all_tilings.extend( ([tiling] * ngood) )
            bands.add(band)
            #tilings.add(tiling)

    print('\nFinished processing all exposures')
    print('bands = ',bands)
    #print('tilings = ',tilings)

    print('n_reject_mean_dt = ',n_reject_mean_dt)
    print('n_reject_mean_de1 = ',n_reject_mean_de1)
    print('n_reject_mean_de2 = ',n_reject_mean_de2)
    print('n_reject_std_dt = ',n_reject_std_dt)
    print('n_reject_std_de1 = ',n_reject_std_de1)
    print('n_reject_std_de2 = ',n_reject_std_de2)
    print('n_reject_rho2 = ',n_reject_rho2)

    # Turn the data into a recarray
    print('all_data.keys = ',all_data.keys())
    formats = ['f8'] * len(all_keys) + ['a1', 'i2']
    #names = all_keys + ['band', 'tiling']
    names = all_keys + ['band']
    data = np.recarray(shape = (len(all_bands),),
                          formats = formats, names = names)
    print('data.dtype = ',data.dtype)
    for key in all_keys:
        data[key] = np.concatenate(all_data[key])
    data['band'] = all_bands
    #data['tiling'] = all_tilings
    print('made recarray')

    tilings = None
    return data, bands, tilings


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

    min_sep = 0.5
    bin_size = 0.2
    bin_slop = 0.1

    results = []
    for (cat1, cat2) in [ (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(dtcat.name, dtcat.name))

        rho = treecorr.KKCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
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
        print('Doing correlation of %s vs %s'%(catnames1, catnames2))
        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
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

    data, bands, tilings = read_data(args, work, limit_bands=args.bands, prefix=prefix)

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
