#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
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


def measure_rho(data,max_sep):
    """Compute the rho statistics
    """
    import numpy
    import treecorr

    ra = data['ra']
    dec = data['dec']
    e1 = data['e1']
    e2 = data['e2']
    s = data['size']
    p_e1 = data['psfex_e1']
    p_e2 = data['psfex_e2']
    p_s = data['psfex_size']

    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                            g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             g1=(e1-p_e1), g2=(e2-p_e2))
    dt = (s**2-p_s**2)/s**2
    print 'mean dt = ',numpy.mean(dt)
    dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             k=dt, g1=dt*e1, g2=dt*e2)

    min_sep = 0.3
    bin_size = 0.2
    bin_slop = 0.3

    results = []
    for (cat1, cat2) in [ (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:

        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        results.append(rho)

    return results


def measure_cross_rho(data,max_sep):
    """Compute the rho statistics
    """
    import numpy
    import treecorr

    ntilings = len(data)

    ra = [ d['ra'] for d in data ]
    dec = [ d['dec'] for d in data ]
    e1 = [ d['e1'] for d in data ]
    e2 = [ d['e2'] for d in data ]
    s = [ d['size'] for d in data ]
    p_e1 = [ d['psfex_e1'] for d in data ]
    p_e2 = [ d['psfex_e2'] for d in data ]
    p_s = [ d['psfex_size'] for d in data ]


    ecats = [ treecorr.Catalog(ra=ra[i], dec=dec[i], ra_units='deg', dec_units='deg', 
                               g1=e1[i], g2=e2[i]) for i in range(ntilings) ]
    decats = [ treecorr.Catalog(ra=ra[i], dec=dec[i], ra_units='deg', dec_units='deg', 
                                g1=(e1[i]-p_e1[i]), g2=(e2[i]-p_e2[i])) for i in range(ntilings) ]
    dt = [ (s[i]**2-p_s[i]**2)/s[i]**2 for i in range(ntilings) ]
    dtcats = [ treecorr.Catalog(ra=ra[i], dec=dec[i], ra_units='deg', dec_units='deg', 
                                k=dt[i], g1=dt[i]*e1[i], g2=dt[i]*e2[i]) for i in range(ntilings) ]

    min_sep = 0.3
    bin_size = 0.2
    bin_slop = 0.3

    results = []
    for (catlist1, catlist2) in [ (decats, decats),
                                  (ecats, decats),
                                  (dtcats, dtcats),
                                  (decats, dtcats),
                                  (ecats, dtcats) ]:
        rho = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='arcmin',
                                     bin_size=bin_size, bin_slop=bin_slop, verbose=2)
        # Avoid all auto correlations:
        for i in range(ntilings):
            for j in range(ntilings):
                if i == j: continue
                if catlist1 is catlist2 and i > j: continue
                rho.process_cross(catlist1[i], catlist2[j])
        varg1 = treecorr.calculateVarG(catlist1)
        varg2 = treecorr.calculateVarG(catlist2)
        rho.finalize(varg1, varg2)
        results.append(rho)

    return results


def main():
    import os
    import json
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

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    #expinfo_file = 'exposure_info_' + args.tag + '.fits'
    expinfo_file = 'exposure_info_y1spte_r_v01.fits'
    print 'reading exposure_info file: ',expinfo_file
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    all_lists = {}
    keys = [ 'ra', 'dec', 'e1', 'e2', 'size', 'psfex_e1', 'psfex_e2', 'psfex_size' ]

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        if expnum not in expinfo['expnum']:
            print 'expnum is not in expinfo!'
            print 'expinfo[expnum] = ',expinfo['expnum']
            raise RuntimeError("Could not find information about this expnum")
        k = numpy.nonzero(expinfo['expnum'] == expnum)[0][0]
        print 'k = ',k
        filter = expinfo['filter'][k]
        print 'filter[k] = ',filter

        tiling = int(expinfo['tiling'][k])
        print 'tiling[k] = ',tiling

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        cat_file = os.path.join(exp_dir, exp + "_psf.fits")
        print 'cat_file = ',cat_file
        try:
            with pyfits.open(cat_file) as pyf:
                data = pyf[1].data
        except:
            print 'Unable to open cat_file %s.  Skipping this file.'%cat_file
            continue

        ccdnums = numpy.unique(data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        stats = []

        mask = (data['flag'] == 0)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        # which filters to add this to
        filters = [ 'griz', filter ]
        if 'g' not in filter: filters = ['riz'] + filters

        # which tilings to add this to
        tilings = [ 0, tiling ]

        print 'filters = ',filters
        print 'tilings = ',tilings
        #print 'all_lists.keys = ',all_lists.keys()
        for filt in filters:
            if filt not in all_lists:
                all_lists[filt] = []
            #print 'all_lists[%s] = '%filt, all_lists[filt]
            for til in tilings:
                while til >= len(all_lists[filt]):
                    all_lists[filt].append({})
                #print 'all_lists[%s][%d] = '%(filt,til), all_lists[filt][til]
                for key in keys:
                    if key not in all_lists[filt][til]:
                        all_lists[filt][til][key] = []
                    #print 'all_lists[%s][%d][%s] = '%(filt,til,key), all_lists[filt][til][key]
                    all_lists[filt][til][key].append(data[key][mask])

    print '\nFinished processing all exposures'

    print 'all filters = ',all_lists.keys()
    for filt in all_lists:
        print 'filter ',filt
        for til in range(len(all_lists[filt])):
            print 'tiling ',til
            if len(all_lists[filt][til]) == 0:
                print 'No files to concatenate for filter %s, tiling %s.'%(filt,til)
                continue
            cols = []
            for key in keys:
                all_lists[filt][til][key] = numpy.concatenate(all_lists[filt][til][key])
                cols.append(pyfits.Column(name=key, format='E', array=all_lists[filt][til][key]))

            cols = pyfits.ColDefs(cols)

            # Depending on the version of pyfits, one of these should work:
            try:
                tbhdu = pyfits.BinTableHDU.from_columns(cols)
            except:
                tbhdu = pyfits.new_table(cols)
            cat_file = "psf_%s_%s.fits"%(filt,til)
            tbhdu.writeto(cat_file, clobber=True)

        # Measure the canonical rho stats using all pairs: (til=0)
        rho1, rho2, rho3, rho4, rho5 = measure_rho(all_lists[filt][0], max_sep=300)

        stat_file = os.path.join(work, "rho_%s.json"%filt)
        stats.append([
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
        ])
        #print 'stats = ',stats
        print 'stat_file = ',stat_file
        with open(stat_file,'w') as fp:
            json.dump(stats, fp)
        print 'Done writing ',stat_file

        # Measure the cross-only rho stats:
        rho1, rho2, rho3, rho4, rho5 = measure_cross_rho(all_lists[filt][1:], max_sep=300)

        stat_file = os.path.join(work, "rho_cross_%s.json"%filt)
        stats.append([
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
        ])
        #print 'stats = ',stats
        print 'stat_file = ',stat_file
        with open(stat_file,'w') as fp:
            json.dump(stats, fp)
        print 'Done writing ',stat_file

    print '\nFinished writing json files'

if __name__ == "__main__":
    main()
