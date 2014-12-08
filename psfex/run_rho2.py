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


def read_blacklists(tag):
    """Read the psfex blacklist file and the other blacklists.

    Returns a dict indexed by the tuple (expnum, ccdnum) with the bitmask value.
    """
    import astropy.io.fits as pyfits
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

def measure_rho(ra,dec,e1,e2,s,m_e1,m_e2,m_s,max_sep):
    """Compute the rho statistics
    """
    import numpy
    import treecorr
    import galsim

    #print 'e1 = ',e1
    #print 'm_e1 = ',m_e1
    #print 'de1 = ',e1-m_e1
    #print 'mean e1 = ',numpy.mean(e1)
    #print 'mean de1 = ',numpy.mean(e1-m_e1)
    #print 'e2 = ',e2
    #print 'm_e2 = ',m_e2
    #print 'de2 = ',e2-m_e2
    #print 'mean e2 = ',numpy.mean(e2)
    #print 'mean de2 = ',numpy.mean(e2-m_e2)
    #print 's = ',s
    #print 'm_s = ',m_s
    #print 'ds = ',s-m_s
    #print 'mean s = ',numpy.mean(s)
    #print 'mean ds = ',numpy.mean(s-m_s)

    # From Barney's paper: http://arxiv.org/pdf/0904.3056v2.pdf
    # rho1 = < (e-em)* (e-em) >     Barney originally called this D1.
    # rho2 = Re < e* (e-em) >       Barney's D2 is actually 2x this.
    # rho3 = < (s-sm) (s-sm) >      Not in Barney's paper, but an obvious extension.
    # rho4 = < s (s-sm) >           Ditto.

    #print 'ra = ',ra
    #print 'dec = ',dec
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                            g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             g1=(e1-m_e1), g2=(e2-m_e2))
    scat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                            k=s)
    dscat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             k=(s-m_s))

    rho1 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho1.process(decat)
    #print 'rho1 = ',rho1.xip
    #print 'rho1.xip_im = ',rho1.xip_im
    #print 'rho1.xim = ',rho1.xim
    #print 'rho1.xim_im = ',rho1.xim_im
    #print 'rho1.sigma = ',numpy.sqrt(rho1.varxi)

    rho2 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho2.process(ecat, decat)
    #print 'rho2 = ',rho2.xip
    #print 'rho2.xip_im = ',rho2.xip_im
    #print 'rho2.xim = ',rho2.xim
    #print 'rho2.xim_im = ',rho2.xim_im
    #print 'rho2.sigma = ',numpy.sqrt(rho2.varxi)

    rho3 = treecorr.KKCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho3.process(dscat)
    #print 'rho3 = ',rho3.xi
    #print 'rho3.sigma = ',numpy.sqrt(rho3.varxi)

    rho4 = treecorr.KKCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho4.process(scat, dscat)
    #print 'rho4 = ',rho4.xi
    #print 'rho4.sigma = ',numpy.sqrt(rho4.varxi)

    return rho1,rho2,rho3,rho4


def main():
    import os
    import glob
    import galsim
    import json
    import numpy
    import astropy.io.fits as pyfits

    args = parse_args()

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

    expinfo_file = 'exposure_info.fits'
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    ra_list = { 'griz' : [], 'riz' : [] }
    dec_list = { 'griz' : [], 'riz' : [] }
    e1_list = { 'griz' : [], 'riz' : [] }
    e2_list = { 'griz' : [], 'riz' : [] }
    s_list = { 'griz' : [], 'riz' : [] }
    pe1_list = { 'griz' : [], 'riz' : [] }
    pe2_list = { 'griz' : [], 'riz' : [] }
    ps_list = { 'griz' : [], 'riz' : [] }

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

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        cat_file = os.path.join(exp_dir, exp + "_psf.fits")
        with pyfits.open(cat_file) as pyf:
            data = pyf[1].data

        ccdnums = numpy.unique(data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        stats = []

        mask = (data['flag'] == 0)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        ra_list['griz'].append(data['ra'][mask])
        dec_list['griz'].append(data['dec'][mask])
        e1_list['griz'].append(data['e1'][mask])
        e2_list['griz'].append(data['e2'][mask])
        s_list['griz'].append(data['size'][mask])

        pe1_list['griz'].append(data['psfex_e1'][mask])
        pe2_list['griz'].append(data['psfex_e2'][mask])
        ps_list['griz'].append(data['psfex_size'][mask])

        if filter not in ra_list.keys():
            ra_list[filter] = []
            dec_list[filter] = []
            e1_list[filter] = []
            e2_list[filter] = []
            s_list[filter] = []
            pe1_list[filter] = []
            pe2_list[filter] = []
            ps_list[filter] = []

        if 'g' not in filter:
            ra_list['riz'].append(data['ra'][mask])
            dec_list['riz'].append(data['dec'][mask])
            e1_list['riz'].append(data['e1'][mask])
            e2_list['riz'].append(data['e2'][mask])
            s_list['riz'].append(data['size'][mask])

            pe1_list['riz'].append(data['psfex_e1'][mask])
            pe2_list['riz'].append(data['psfex_e2'][mask])
            ps_list['riz'].append(data['psfex_size'][mask])

        ra_list[filter].append(data['ra'][mask])
        dec_list[filter].append(data['dec'][mask])
        e1_list[filter].append(data['e1'][mask])
        e2_list[filter].append(data['e2'][mask])
        s_list[filter].append(data['size'][mask])

        pe1_list[filter].append(data['psfex_e1'][mask])
        pe2_list[filter].append(data['psfex_e2'][mask])
        ps_list[filter].append(data['psfex_size'][mask])

    print '\nFinished processing all exposures'

    for key in ra_list.keys():
        ra = numpy.concatenate(ra_list[key])
        dec = numpy.concatenate(dec_list[key])
        e1 = numpy.concatenate(e1_list[key])
        e2 = numpy.concatenate(e2_list[key])
        s = numpy.concatenate(s_list[key])
        pe1 = numpy.concatenate(pe1_list[key])
        pe2 = numpy.concatenate(pe2_list[key])
        ps = numpy.concatenate(ps_list[key])

        rho1, rho2, rho3, rho4 = measure_rho(ra,dec,e1,e2,s,pe1,pe2,ps,
                                             max_sep=300)

        stat_file = os.path.join(work, "rho_"+key+".json")
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
            rho3.xi.tolist(),
            rho3.varxi.tolist(),
            rho4.xi.tolist(),
            rho4.varxi.tolist()
        ])
        #print 'stats = ',stats
        with open(stat_file,'w') as f:
            json.dump(stats, f)

    print '\nFinished writing json files'

if __name__ == "__main__":
    main()
