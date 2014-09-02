#! /usr/bin/env python
# Program to run single file in psfex

import argparse,os,glob,re,pyfits,random,copy
import numpy as np

parser = argparse.ArgumentParser(description='Run single file')

# Directory arguments
parser.add_argument('--cat_dir', default='/astro/u/rarmst/soft/bin/',
                    help='location of sextrator executable')
parser.add_argument('--psf_dir', default='/astro/u/rarmst/soft/bin/',
                    help='location of psfex executable')
parser.add_argument('--findstars_dir', default='/astro/u/mjarvis/bin',
                    help='location wl executables')
parser.add_argument('--output', default='./',
                    help='location of outputs')
parser.add_argument('--clear_output', default=False, action='store_const', const=True,
                    help='should the output directory be cleared before writing new files?')

# Exposure inputs
parser.add_argument('--file', default=None,
                    help='name of file')
parser.add_argument('--exp_match', default='',
                    help='regexp to search for files in exp_dir')

parser.add_argument('--exps', default='', nargs='+',
                    help='list of exposures to run')
parser.add_argument('--runs', default='', nargs='+',
                    help='list of runs')

# Configuration files
parser.add_argument('--config_cat',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/default.sex',
                    help='sextractor config file')
parser.add_argument('--config_psf',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/new.psfex',
                    help='psfex config file')
parser.add_argument('--config_findstars',
                    default='wl.config +wl_desdm.config +wl_firstcut.config',
                    help='wl config file')
parser.add_argument('--param_file',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.param_psfex',
                    help='sextractor param file')
parser.add_argument('--filt_file',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.conv',
                    help='name of sextractor filter file')
parser.add_argument('--star_file',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/sex.nnw',
                    help='name of sextractor star file')
parser.add_argument('--tapebump_file',
                    default='/astro/u/mjarvis/rmjarvis/DESWL/psfex/mask_ccdnum.txt',
                    help='name of tape bump file')
parser.add_argument('--make_symlinks', default=True, action='store_const', const=False,
                    help='make symlinks from $DESDATA/EXTRA/red/$run/psfex-rerun/$exp')

# options
parser.add_argument('--rm_files', default=True, action='store_const', const=False,
                    help='remove unpacked files after finished')
parser.add_argument('--run_psfex', default=True, action='store_const', const=False,
                    help='run psfex on files')
parser.add_argument('--use_findstars', default=False, action='store_const', const=True,
                    help='use findstars results in psfex')
parser.add_argument('--mag_cut', default=-1, type=float,
                    help='remove the top mags using mag_auto')
parser.add_argument('--nstars', default=10, type=int,
                    help='use median of brightest nstars for min mag')
parser.add_argument('--use_tapebumps', default=True, action='store_const', const=False,
                    help='avoid stars in or near tape bumps')
parser.add_argument('--tapebump_extra', default=2, type=float,
                    help='How much extra room around tape bumps to exclude stars in units of FWHM')
parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                    help='Only do 1 ccd per exposure (used for debugging)')


args = parser.parse_args()


# Read and parse the tapebump file if we are going to need it.
if args.use_tapebumps:
    raw_tbdata = np.genfromtxt(args.tapebump_file, delimiter=',')
    # repackage this as  dict (key = ccdnum) of lists of tuples (ymin, xmin, ymax, xmax)
    tbdata = {}
    for d in raw_tbdata:
        ccdnum = int(d[0])
        if ccdnum not in tbdata:
            tbdata[ccdnum] = []
        tbdata[ccdnum].append( (int(d[1]), int(d[2]), int(d[3]), int(d[4])) )
    print 'read in tape bump file.  %d bumps for %d chips'%(len(raw_tbdata),len(tbdata))

# This will be the function to process the tapebump data for a particular chip:
def exclude_tapebumps(tbd, data, extra):
    """
    tbd = tbdata for a particular ccdnum,
    data = the input data for the stars
    extra = how much extra distance around the tape bumps to exclude stars in pixels
    """
    # I'm sure it doesn't matter, but add an extra 0.5 pixel to the slop because the tape bump
    # values are in terms of which pixels to include as part of the tape bump.  So the edges
    # are an extra 0.5 pixel outside of that.
    extra += 0.5
    x = data['X_IMAGE']
    y = data['Y_IMAGE']
    masks = [(y>tb[0]-extra) & (x>tb[1]-extra) & (y<tb[2]+extra) & (x<tb[3]+extra) for tb in tbd]
    mask = np.any(masks, axis=0)
    if sum(mask) > 0:
        print '   masking %d stars for being in or near a bump'%sum(mask)
        print '       tapebumps = ',tbd
        print '       excluded x,y = ',zip(x[mask],y[mask])
    return data[np.logical_not(mask)]

# How many stars are too few or too many?
FEW_STARS = 50
MANY_STARS = 500
# How high is a high FWHM?  1.8 arcsec / 0.26 arcsec/pixel = 6.9 pixels
HIGH_FWHM = 6.9

# flag values
NO_STARS_FLAG = 1
TOO_FEW_STARS_FLAG = 2
TOO_MANY_STARS_FLAG = 4
TOO_HIGH_FWHM_FLAG = 8
ERROR_FLAG = 16

class NoStarsException(Exception):
    pass

blacklist_file = os.path.join(args.output,'psfex_blacklist')

def log_blacklist(run, exp, ccdnum, flag):
    try:
        with open(blacklist_file,'a') as f:
            f.write("%s %s %d %d\n"%(run,exp,ccdnum,flag))
    except OSError:
        print 'Error opening blacklist.  Wait and try again.'
        import time
        time.sleep(1)
        return log_blacklist(file,run,exp,ccdnum,flag)

# Make the output directory if it does not exist yet.
try:
    os.mkdir(args.output)
except OSError:
    pass

for run,exp in zip(args.runs,args.exps):

    print 'run, exp = ',run,exp

    odir = os.path.join(args.output,exp)
    logfile = os.path.join(odir, 'log.'+exp)

    if args.clear_output:
        import shutil
        try:
            shutil.rmtree(odir)
        except OSError:
            pass

    datadir='/astro/u/astrodat/data/DES'
    input_dir = os.path.join(datadir,'OPS/red/%s/red/%s/'%(run,exp))
    files=[]
    
    # We could build the file name but we might as just well search for 
    # what is there
    files = glob.glob('%s/%s'%(input_dir,args.exp_match))

    try:
        os.makedirs(odir)
    except:
        if not os.path.exists(odir): raise

    for file in files:
        print 'Processing ', file
        flag = 0

        # find out if the file is fpacked by the extension
        dir, base_file = os.path.split(file)

        # find the base filename
        if os.path.splitext(base_file)[1] == '.fz':
            do_unpack = True
            base_file=os.path.splitext(base_file)[0]
        else:
            do_unpack = False
        print '   dir, base_file = ',dir,base_file
    
        # remove the .fits extension for the root
        match=re.search('(.*)\.fits.*',base_file)

        if match:
            root=match.group(1)
        else:
            print "Cannot find base name for "+base_file+" skipping"
            continue

        # Get the chipnum.  We could do this by reading the file and looking for the 
        # CCDNUM keyword.  But it's simpler to just use the filename.  We should be
        # able to assume that the file names are not mangled.
        ccdnum = int(root.split('_')[-1])
        print '   ccdnum = ',ccdnum

        try:

            if do_unpack:
                funpack_file=odir+'/'+base_file
            
                # If an unpacked file does not exist in the output directory 
                if not os.path.exists(funpack_file):
                    print '   unpacking fz file'
                    cmd='funpack -O %s %s  >>%s 2>&1' % (funpack_file,file,logfile)
                    ok=os.system(cmd)
                img_file=funpack_file
                
            else:
                # If the file is not fpacked, make a symlink into the output directory
                if not os.path.exists(odir+'/'+base_file):
                    # make symlink to local directory
                    os.system('ln -s %s %s  >>%s 2>&1'%(file,odir,logfile))
                img_file=odir+'/'+base_file
            hdu=0

            # extract the saturation level, this is how desdm runs sextractor
            # we need the fwhm for class star
            # Also need the fwhm for doing the tape bumps.
            with pyfits.open(img_file) as pyfile:
                sat=-1
                fwhm=4.
                try:
                    sat=pyfile[hdu].header['SATURATE']
                    fwhm=pyfile[hdu].header['FWHM']
                except:
                    print 'Cannot process ',img_file
                    continue
                print '   fwhm = ',fwhm
            if fwhm > HIGH_FWHM:
                flag |= TOO_HIGH_FWHM_FLAG

            # This is the file that holds the vignettes 
            psfcat_file=odir+'/'+root+'_psfcat.fits'

            # if the sextractor catalog does not exist 
            if not os.path.exists(psfcat_file):
                print '   running sextractor'
                cat_cmd="{cat_dir}/sex {img_file}[0] -c {cat_config} -CATALOG_NAME {cat_file} -CATALOG_TYPE FITS_LDAC -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {img_file}[2] -PARAMETERS_NAME {param_file} -FILTER_NAME {filter_file}  -STARNNW_NAME {star_file} -DETECT_MINAREA 3 -SEEING_FWHM {fwhm} -SATUR_LEVEL {sat}  >>{logfile} 2>&1".format(
                    cat_dir=args.cat_dir, img_file=img_file, cat_config=args.config_cat,
                    cat_file=psfcat_file,param_file=args.param_file,filter_file=args.filt_file,
                    star_file=args.star_file,fwhm=fwhm,logfile=logfile,sat=sat)
                os.system(cat_cmd)

            # if we want to use only the stars selected by findstars
            if args.use_findstars:
        
                star_file=odir+'/'+root+'_findstars.fits'    

                # run find stars
                print '   running findstars'
                findstars_cmd='%s/findstars %s root=%s cat_ext=_psfcat.fits stars_file=%s input_prefix=%s/ >>%s 2>&1'%(
                        args.findstars_dir,args.config_findstars,root,star_file,odir,logfile)
                os.system(findstars_cmd)
            
                wlcat_file=psfcat_file.replace('psfcat','psfcat_findstars')
                initfile=pyfits.open(psfcat_file)
            
                findstars_file=pyfits.open(star_file)
                mask=findstars_file[1].data['star_flag']==1
                print '   found %d stars'%np.count_nonzero(mask)
                # create new sextractor file with only these entries
                data=initfile[2].data[mask]
                if len(data) == 0:
                    # Can't really do any of the rest of this, so skip out to the end.
                    raise NoStarsException()
                if len(data) < FEW_STARS:
                    flag |= TOO_FEW_STARS_FLAG
                if len(data) > MANY_STARS:
                    flag |= TOO_MANY_STARS_FLAG
            
                # Need to make different copy of these to not fail
                hdu1=copy.copy(initfile[0])
                hdu2=copy.copy(initfile[1])
            
                hdu = pyfits.BinTableHDU(data)
                hdu.name='LDAC_OBJECTS'
                list=pyfits.HDUList([hdu1,hdu2, hdu])
                list.writeto(wlcat_file,clobber=True)
        
                if args.rm_files:
                    os.system('rm %s >>%s 2>&1'%(psfcat_file,logfile))
                
                    # assign the psfcat_file to the new file
                    psfcat_file=wlcat_file

            # If we want to cut the brighest magnitudes
            if args.mag_cut>0 or args.use_tapebumps:

                # get the brightest 10 stars that have flags=0 and take the median just in case some
                # strange magnitudes were selected
                hdu=2
                pyfile=pyfits.open(psfcat_file)
                flags_mask=pyfile[hdu].data['FLAGS']==0
                data=pyfile[hdu].data
                print '   len(data) = ',len(data)

                if args.mag_cut > 0:
                    magcut_file=psfcat_file.replace('psfcat','psfcat_magcut_%0.1f'%args.mag_cut)
                    mags=data['MAG_AUTO'][flags_mask]
                    mags.sort()
                    min_star=np.median(mags[0:args.nstars])
                    print '   min mag = ',mags[0]
                    print '   median of brightest %d is '%args.nstars,min_star
            
                    mag_mask=data['MAG_AUTO']>min_star+args.mag_cut
                    print '   select stars dimmer than ',min_star+args.mag_cut
                    print '   which includes %d stars'%np.count_nonzero(mag_mask)
            
                    data=data[mag_mask]
                    print '   after exclude bright: len(data) = ',len(data)

                if args.use_tapebumps:
                    data = exclude_tapebumps(tbdata[ccdnum], data, args.tapebump_extra * fwhm)
                    print '   after exclude tapebumps: len(data) = ',len(data)

                # Recheck this.
                if len(data) < FEW_STARS:
                    flag |= TOO_FEW_STARS_FLAG
        
                # Need to make different copy of these to not fail
                hdu1=copy.copy(pyfile[0])
                hdu2=copy.copy(pyfile[1])
            
                hdu = pyfits.BinTableHDU(data)
                hdu.name='LDAC_OBJECTS'
                list=pyfits.HDUList([hdu1,hdu2, hdu])
                list.writeto(magcut_file,clobber=True)
            
                if args.rm_files:
                    os.system('rm %s >>%s 2>&1'%(psfcat_file,logfile))
                
                    # assign the psfcat_file to the new file
                    psfcat_file=magcut_file

    
            psf_file=psfcat_file.replace('fits','psf')
            catout_file=psfcat_file.replace('fits','used.fits')
    
            if args.run_psfex:
                print '   running psfex'
                psf_cmd=('%s/psfex %s -c %s -OUTCAT_TYPE FITS_LDAC -OUTCAT_NAME %s >>%s 2>&1' % (
                            args.psf_dir,psfcat_file,args.config_psf,catout_file,logfile))
                os.system(psf_cmd)
        
            if do_unpack and args.rm_files:
                rm_cmd='rm %s >>%s 2>&1'%(img_file,logfile)
                os.system(rm_cmd)
        
            if args.make_symlinks:
                print '   making symlink'
                link_dir = os.path.join(datadir,'EXTRA/red/%s/psfex-rerun/%s/'%(run,exp))
                print '   link_dir = ',link_dir
                try: 
                    os.makedirs(link_dir, mode=0775)
                    print '   made directory ',link_dir
                except OSError:
                    if not os.path.isdir(link_dir): raise
                psf_link = os.path.join(link_dir,os.path.basename(psf_file))
                try:
                    os.remove(psf_link)
                    print '   removed existing link: ',psf_link
                except OSError:
                    pass
                os.symlink(psf_file,psf_link)
                print '   made link: ',psf_link
                catout_link = os.path.join(link_dir,os.path.basename(catout_file))
                try:
                    os.remove(catout_link)
                    print '   removed existing link: ',catout_link
                except OSError:
                    pass
                os.symlink(catout_file,catout_link)
                print '   made link: ',catout_link

        except NoStarsException:
            print 'No stars.  Log this in the blacklist and continue.'
            flag |= NO_STARS_FLAG
        except Exception as e:
            print 'Caught exception: ',e
            print 'Log this in the blacklist and continue.'
            flag |= ERROR_FLAG

        if flag:
            log_blacklist(run,exp,ccdnum,flag)

        if args.single_ccd:
            break


