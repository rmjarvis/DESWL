#! /usr/bin/env python
# Verify that all PSFEx files were either made or blacklisted
# Typical usage:
#    python verify.py --tag=y1a1-v02 --file=y1all

redo_exp = set()
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Build PSF catalogs for a set of runs/exposures')

    # Drectory arguments
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

    args = parser.parse_args()
    return args


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


def read_blacklists(tag):
    """Read the psfex blacklist file and the other blacklists.

    Returns a dict indexed by the tuple (expnum, ccdnum) with the bitmask value.
    """
    d = {}  # The dict will be indexed by (expnum, ccdnum)
    print 'reading blacklists'

    if False:
        import astropy.io.fits as pyfits
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

    args = parse_args()

    datadir = '/astro/u/astrodat/data/DES'

    print 'Reading blacklist files for tag = %s'%args.tag
    flag_dict = read_blacklists(args.tag)

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    for run,exp in zip(runs,exps):

        print 'Verifying exposure ',exp
        expnum = int(exp[6:])

        # The input directory from the main DESDM reduction location.
        input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))

        # This is where the PSFEx files should be.
        extra_dir = os.path.join(datadir,'EXTRA/red/%s/psfex-rerun/%s/%s/'%(run,args.tag,exp))

        # Get the file names in that directory.
        files = glob.glob('%s/%s'%(input_dir,args.exp_match))

        for file_name in files:
            #print '\nChecking ', file_name

            try:
                desdm_dir, root, ccdnum = parse_file_name(file_name)
            except:
                print '   Unable to parse file_name %s.  Skipping this file.'%file_name
                continue

            key = (expnum, ccdnum)
            if key in flag_dict:
                black_flag = flag_dict[key]
                #print '   blacklist flag = ',black_flag
                if black_flag & (113 << 15):
                    #print '   Catastrophic flag.  Skipping this file.'
                    continue
            else:
                black_flag = 0

            # Check that the .psf file exists
            psfex_file_name = os.path.join(extra_dir, root + '_psfcat.psf')

            if not os.path.exists(psfex_file_name):
                print '   PSFEx file %s does not exist.'%psfex_file_name
                if black_flag:
                    print '   In blacklist, but not catastrophic.'
                else:
                    redo_exp.add( (run,exp) )

    print '\nFinished verifying all exposures'
    if len(redo_exp) > 0:
        print '%s exposures have missing PSFEx files'
        with open('redo_exp', 'w') as fout:
            for run, exp in redo_exp:
                fout.write('%s  %s\n'%(run, exp))
        print 'Wrote list of %d exposures to redo to redo_exp.'%len(redo_exp)
    else:
        print 'All PSFEx files verified.'

if __name__ == "__main__":
    main()
