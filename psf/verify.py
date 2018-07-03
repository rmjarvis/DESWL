#! /usr/bin/env python
# Verify that all PSF files were either made or blacklisted
# Typical usage:
#    python verify.py --tag=y1a1-v02 --file=y1all

import os
import sys
import glob
import logging
import fitsio
import pandas
import numpy as np
import numpy.lib.recfunctions as nprec
import galsim

# flag values for blacklist
NO_STARS_FLAG = 1
TOO_FEW_STARS_FLAG = 2
TOO_MANY_STARS_FLAG = 4
TOO_HIGH_FWHM_FLAG = 8
FINDSTARS_FAILURE = 16
PSF_FAILURE = 32
ERROR_FLAG = 64
LARGE_SIZE_SPREAD = 128

def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Verify PSF catalogs for a set of exposures')

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
    parser.add_argument('--blacklist', default=1, type=int,
                        help='add failed CCDs to the blacklist')

    args = parser.parse_args()
    return args


def main():

    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger('verify')

    args = parse_args()

    blacklist_file = '/astro/u/astrodat/data/DES/EXTRA/blacklists/psf'
    if args.tag:
        blacklist_file += '-' + args.tag
    blacklist_file += '.txt'

    work = os.path.expanduser(args.work)
    if work.endswith('y3_piff') and args.tag:
        work += '/' + args.tag

    if args.file != '':
        logger.info('Read file %s',args.file)
        with open(args.file) as fin:
            exps = { line.strip() for line in fin if line[0] != '#' }
        logger.info('File includes %d exposures',len(exps))
    elif args.exps is not None:
        exps = set(args.exps)
        logger.info('Explicit listing of %d exposures',len(exps))
    else:
        raise RuntimeError("Either file or exps is required")

    redo_exp = set()
    blacklist = set()

    all_expcat = []
 
    for exp in sorted(exps):
        expnum = int(exp)

        expname = os.path.join(work, exp, 'exp_psf_cat_%d.fits'%expnum)
        if not os.path.exists(expname):
            logger.warning('%s does not exist.  Add %s to the redo set',expname, exp)
            redo_exp.add(exp)
            # Make sure we actually redo everything in this exp by removing the psf_cat files.
            for psf_cat_file in glob.glob(os.path.join(work, exp, 'psf_cat_*.fits')):
                os.remove(psf_cat_file)
            continue

        logger.info("Reading %s",expname)
        try:
            expcat = fitsio.read(expname, ext='info')
            starcat = fitsio.read(expname, ext='stars')
        except Exception as e:
            logger.warning("Caught %s",e)
            logger.warning("Add %s to the redo set",exp)
            redo_exp.add(exp)
            continue
        deg = galsim.degrees
        ra = np.mean([ (ra * deg).wrap(0 * deg) for ra in starcat['ra'] ]) / deg
        dec = np.mean(starcat['dec'])
        expcat_with_radec = nprec.append_fields(expcat, ('ra', 'dec'), (ra, dec), usemask=False)
        all_expcat.append(expcat_with_radec)

        for row in expcat:
            ccdnum = row['ccdnum']
            flag = row['flag']
            psf_cat_file = os.path.join(work, exp, 'psf_cat_%d_%d.fits'%(expnum,ccdnum))

            if flag != 0:
                logger.info('  add (%s, %s) to blacklist, flag = %s',exp,ccdnum,flag)
                blacklist.add( (exp, ccdnum, flag) )
                if flag & ERROR_FLAG != 0:
                    logger.info('  ERROR_FLAG.  Redo this ccd (%s, %s)',exp,ccdnum)
                    redo_exp.add(exp)
                    if os.path.exists(psf_cat_file):
                        os.remove(psf_cat_file)
                continue

            if not os.path.exists(psf_cat_file):
                logger.info('  %d: %s does not exist.  Add to the redo set',ccdnum,psf_cat_file)
                blacklist.add( (exp, ccdnum, PSF_FAILURE) )
                redo_exp.add(exp)
                continue

            # Check that the Piff file exists
            piff_file_name = row['piff_file']
            if not os.path.exists(piff_file_name):
                logger.info('   %d: Piff file %s does not exist.',ccdnum,piff_file_name)
                blacklist.add( (exp, ccdnum, PSF_FAILURE) )
                redo_exp.add(exp)
                os.remove(psf_cat_file)
                continue

    logger.info('\nFinished verifying all exposures')
    if len(redo_exp) > 0:
        logger.info('%s exposures have missing Piff files',len(redo_exp))
        with open('redo_exp', 'w') as fout:
            for exp in sorted(redo_exp):
                fout.write('%s\n'%exp)
        logger.info('Wrote list of %d exposures to redo to redo_exp.',len(redo_exp))
    else:
        logger.info('All PSF files verified.')
        all_expcat = np.concatenate(all_expcat)
        all_expname = os.path.join(work, '%s_info.fits'%args.tag)
        fitsio.write(all_expname, all_expcat, clobber=True)
        logger.info('Wrote info file to %s',all_expname)

    if args.blacklist:
        logger.info('Logging blacklisted chips to %s',blacklist_file)
        if os.path.exists(blacklist_file):
            logger.info("Removing existing blacklist file")
            os.remove(blacklist_file)
        with open(blacklist_file,'w') as f:
            logger.info("Writing %d entries into blacklist file",len(blacklist))
            for exp, ccd, flag in sorted(blacklist):
                f.write("%s %s %s\n"%(exp,ccd,flag))


if __name__ == "__main__":
    main()
