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
    # rho3 = < (s^2-sm^2)/s^2 (s^2-sm^2)/s^2 >  Not in Barney's paper

    #print 'ra = ',ra
    #print 'dec = ',dec
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                            g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             g1=(e1-m_e1), g2=(e2-m_e2))
    dt = (s**2-m_s**2)/s**2
    print 'mean dt = ',numpy.mean(dt)
    dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', 
                             k=dt, g1=dt*e1, g2=dt*e2)

    rho1 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho1.process(decat)

    rho2 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho2.process(ecat, decat)

    rho3 = treecorr.KKCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho3.process(dtcat)

    rho4 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho4.process(dtcat, decat)

    rho5 = treecorr.GGCorrelation(min_sep=0.5, max_sep=max_sep, sep_units='arcmin',
                                  bin_size=0.1, verbose=1)
    rho5.process(dtcat)

    return rho1,rho2,rho3,rho4,rho5


def main():
    import os
    import glob
    import galsim
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

    expinfo_file = 'exposure_info_' + args.tag + '.fits'
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
        print key
        ra = numpy.concatenate(ra_list[key])
        dec = numpy.concatenate(dec_list[key])
        e1 = numpy.concatenate(e1_list[key])
        e2 = numpy.concatenate(e2_list[key])
        s = numpy.concatenate(s_list[key])
        pe1 = numpy.concatenate(pe1_list[key])
        pe2 = numpy.concatenate(pe2_list[key])
        ps = numpy.concatenate(ps_list[key])

        rho1, rho2, rho3, rho4, rho5 = measure_rho(ra,dec,e1,e2,s,pe1,pe2,ps, max_sep=300)

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
        with open(stat_file,'w') as f:
            json.dump(stats, f)

    print '\nFinished writing json files'

if __name__ == "__main__":
    main()
