import galsim
from galsim import pyfits
import os
import sys
import math
import numpy
import shutil
import subprocess

# Setup various file and directory names:
work_dir = '/direct/astro+astronfs03/workarea/mjarvis'
tile_name = 'DES0436-5748'
out_dir = os.path.join(work_dir,tile_name)
print 'tile = ',tile_name
print 'out_dir = ',out_dir

meds_dir = os.path.join('/astro/u/astrodat/data/DES/meds/011/20130820000021_DES0436-5748')
meds_file = tile_name + '-r-meds-011.fits.fz'
print 'meds_file = ',meds_file

# Open the current meds file:
meds = pyfits.open(os.path.join(meds_dir,meds_file))
print 'Opened meds file'

# First pull out some information from the meta catalog
meta = meds[3].data[0]
magzp_ref = meta['magzp_ref']           # 30
se_hdu = meta['se_hdu']                 # 2
se_wt_hdu = meta['se_wt_hdu']           # 4
se_badpix_hdu = meta['se_badpix_hdu']   # 3
sky_hdu = meta['sky_hdu']               # 2
seg_hdu = meta['seg_hdu']               # 2
coadd_hdu = meta['coadd_hdu']           # 2
coadd_wt_hdu = meta['coadd_wt_hdu']     # 3
coadd_seg_hdu = meta['coadd_seg_hdu']   # 0  NB This is a bug.
fake_coadd_seg = meta['fake_coadd_seg'] # 0
des_root = meta['DESDATA']              # /astro/u/astrodat/data/DES
cat_file = meta['cat_file']             # ...-r-meds-input-011.dat
coadd_file = meta['coadd_file']         # ...-r.fits.fz
coadd_srclist = meta['coadd_srclist']   # ...-r-meds-srclist-011.dat
coaddcat_file = meta['coaddcat_file']   # ...-i_cat.fits
coaddseg_file = meta['coaddseg_file']   # ..._r_seg.fits.fz
print 'Read meta information'

# Note: all the above hdu numbers use 1-based numbering (ala cfitsio).  Galsim and pyfits both
# use 0-based, so need to subtract 1 from each.
se_hdu -= 1
se_wt_hdu -= 1
se_badpix_hdu -= 1
sky_hdu -= 1
seg_hdu -= 1
coadd_hdu -= 1
coadd_wt_hdu -= 1

# Read in the coadd files
coadd_image = galsim.fits.read(coadd_file, hdu=coadd_hdu)
coadd_weight = galsim.fits.read(coadd_file, hdu=coadd_wt_hdu)
coadd_segmap = galsim.fits.read(coaddseg_file, hdu=1)  # NB coadd_seg_hdu is wrong!
coadd_cat = pyfits.open(coaddcat_file)[1].data
print 'Opened all coadd files'
print 'nobjects = ',len(coadd_cat)

# Make arrays of the single-epoch images
image_cat = meds[2].data
image_path = image_cat['image_path']
sky_path = image_cat['sky_path']
seg_path = image_cat['seg_path']
magzp = image_cat['magzp']
scale = image_cat['scale']
nimages = len(image_cat)  # 138
print 'nimages = ',nimages

# Make the new srclist:
# Take the existing one and just update the image names.
with open(coadd_srclist,'r') as src:
    src_rows = [ line.split() for line in src ]
src_cols = zip(*src_rows)
out_path = [ os.path.join(out_dir,os.path.basename(file)) for file in image_path ]
# Copy over the sky file to our output directory to get everything in one place.
for file in sky_path[1:]:
    shutil.copy2(file, out_dir)
sky_path = [ os.path.join(out_dir,os.path.basename(file)) for file in sky_path ]
# We make new seg files now, so use out_dir for those paths as well.  But don't copy the originals.
seg_path = [ os.path.join(out_dir,os.path.basename(file)) for file in seg_path ]
src_cols[0] = list(out_path[1:])
src_cols[1] = list(sky_path[1:])
src_cols[2] = list(seg_path[1:])
src_rows = zip(*src_cols)
out_src_file = os.path.join(out_dir,os.path.basename(coadd_srclist))
with open(out_src_file,'w') as out_src:
    for row in src_rows:
        out_src.write('%s %s %s %s\n'%(row))

# Read the image files and make blank ones.
# In principal, we could read the files and get the bounds from that, but 
# that's slow.  So we just use the known bounds and use that for everything.
se_bounds = galsim.BoundsI(1,2048,1,4096)
image_wcs = [ galsim.FitsWCS(file_name) for file_name in image_path ]
images = [ galsim.Image(wcs=wcs, bounds=se_bounds) for wcs in image_wcs ]
# Except the first one (the coadd image) is different
coadd_bounds = galsim.BoundsI(1,10000,1,10000)
images[0] = galsim.Image(wcs=image_wcs[0], bounds=coadd_bounds)
print 'Initialized the blank images'

# Setup a PSF that varies according to a second order function across the image
# Each chip will have its own PSF functions for size, e1, e2.
# size = 0.9 arcsec + 0.2 arcsec * F_size(x,y)
# e1 = 0.05 * F_e1(x,y)
# e2 = 0.05 * F_e2(x,y)
# where each F has the form: a + b T1(x) + c T1(y) + d T2(x) + e T1(x) T1(y) + f T2(y)
# Similar for e1, e2.  So need 18 parameters.  Each parameter is chosen to 
# fall between -1 and 1, and the polynomial factors are Chebyshev polynomials
# using the appropriate bounds (different for x,y) also vary between -1 and 1.
ud = galsim.UniformDeviate(8675309)
psf_params = [ [ 2.*ud()-1. for i in range(18) ] for j in range(nimages) ]
print 'Initialized psf parameters'

# Finally use the coadd_cat and the one in the meds file to draw each galaxy on the appropriate
# new image file.
meds_cat = meds[1].data

nobj = len(coadd_cat)  # 60505
print 'nobj = ',nobj
assert nobj == len(meds_cat)
min_flux = 1.e100

input_flags = coadd_cat['FLAGS']
mag = coadd_cat['MAG_AUTO']
flux = coadd_cat['FLUX_AUTO']
ixx = coadd_cat['X2WIN_IMAGE']
ixy = coadd_cat['XYWIN_IMAGE']
iyy = coadd_cat['Y2WIN_IMAGE']
hlr = coadd_cat['FLUX_RADIUS'] * 1.18 # This is approximate, since FLUX_RADIUS is based on a 
                                      # Gaussian sigma, not half-light radius. But close enough.
e1 = (ixx-iyy)/(ixx+iyy)
e2 = 2.*ixy/(ixx+iyy)

print 'nobj with flags != 0 = ',(input_flags != 0).sum()
print 'nobj with mag > 24 = ',(mag > 24).sum()
print 'nobj with flux <= 0 = ',(flux <= 0).sum()
print 'nobj with ixx <= 0 = ',(ixx <= 0).sum()
print 'nobj with iyy <= 0 = ',(iyy <= 0).sum()
print 'nobj with hlr <= 0 = ',(hlr <= 0).sum()
print 'nobj with hlr > 8 = ',(hlr > 8).sum()
print 'nobj with e > 0.8 = ',(e1*e1 + e2*e2 > 0.8**2).sum()

# Set flags indicating why an object was masked:
flags = numpy.zeros(len(coadd_cat), dtype=int)
flags[ input_flags != 0 ] |= 1                      # 1 = input flag
flags[ (mag > 24) | (flux <= 0) ] |= 2              # 2 = too faint
flags[ (ixx <= 0) | (iyy <= 0) | (hlr <= 0)  ] |= 4 # 4 = bad size
flags[ hlr > 8 ] |= 8                               # 8 = too large
flags[ e1*e2 + e2*e2 > 0.8**2 ] |= 16               # 16 = too elliptical

print 'nobj passing all masks = ',(flags == 0).sum()

nf_pm1 = 0
nf_tot = 0

# Use spread model to decide which things we will call stars.
is_star = coadd_cat['SPREAD_MODEL'] < 0.003
print 'nstars = ',is_star.sum()
print 'nstars passing masks = ',is_star[flags==0].sum()

# Set up a place to store the truth values.  Initialize with -999
true_g1 = numpy.empty(len(coadd_cat))
true_g2 = numpy.empty(len(coadd_cat))
true_hlr = numpy.empty(len(coadd_cat))
true_g1[:] = -999.
true_g2[:] = -999.
true_hlr[:] = -999.
nexp = numpy.zeros(len(coadd_cat), dtype=int)
psf_fwhm = numpy.zeros(len(coadd_cat))
psf_g1 = numpy.zeros(len(coadd_cat))
psf_g2 = numpy.zeros(len(coadd_cat))
wcs_scale = numpy.zeros(len(coadd_cat))
wcs_g1 = numpy.zeros(len(coadd_cat))
wcs_g2 = numpy.zeros(len(coadd_cat))
wcs_theta = numpy.zeros(len(coadd_cat))

do_print = False    # Allows the ability to turn on debug printing dynamically.

for obj_num in range(len(coadd_cat)):
    coadd_info = coadd_cat[obj_num]
    meds_info = meds_cat[obj_num]
    number = coadd_info['NUMBER']

    if do_print: print '\ncoadd number ',number
    if flags[obj_num] != 0: continue

    print 'coadd number ',number
    assert coadd_info['NUMBER'] == meds_info['number']

    # Should already have skipped objects with flags
    assert coadd_info['FLAGS'] == 0

    # For the position, use the WIN value.  Supposedly, this is the most accurate one.
    ra = coadd_info['ALPHAWIN_J2000'] * galsim.degrees
    dec = coadd_info['DELTAWIN_J2000'] * galsim.degrees
    world_pos = galsim.CelestialCoord(ra,dec)
    if do_print: print 'world_pos = ',ra,dec,world_pos

    # Determine if this is a star
    spread = coadd_info['SPREAD_MODEL']
    if do_print: print 'spread = ',spread,' is_star? ',is_star[obj_num]

    if do_print: print 'mag = ',coadd_info['MAG_AUTO']
    assert coadd_info['MAG_AUTO'] <= 24

    flux = coadd_info['FLUX_AUTO']
    if do_print: print 'flux = ',flux,type(flux)
    assert flux > 0
    if flux < min_flux: 
        min_flux = flux

    if not is_star[obj_num]:
        # Get the parameters for building the galaxy
        ixx = coadd_info['X2WIN_IMAGE']
        ixy = coadd_info['XYWIN_IMAGE']
        iyy = coadd_info['Y2WIN_IMAGE']
        hlr = coadd_info['FLUX_RADIUS'] * 1.18 
        # This is approximate, since FLUX_RADIUS is based on a Gaussian sigma, not half-light 
        # radius.  But close enough.
        if do_print: print 'hlr = ',hlr,type(hlr)
        assert ixx > 0. and iyy > 0. and hlr > 0.
        assert hlr <= 8

        # These are in pixel units.  We want arcsec, so scale by the pixel size.
        # This ignores the distortion, but it's fine, since these are just reference
        # values -- we don't really care if they are a little bit wrong.  
        # Another way they are wrong is that they are the _observed_ size and shape, but 
        # we are going to use them as though they are the intrinsic galaxy size and shape.
        hlr *= 0.26

        # Build the galaxy based on the coadd values
        gal = galsim.Exponential(half_light_radius = float(hlr), flux = float(flux))
        e1 = (ixx-iyy)/(ixx+iyy)
        e2 = 2.*ixy/(ixx+iyy)
        if do_print: print 'e1,e2 = ',e1,e2
        assert e1*e1 + e2*e2 <= 0.8
        shear = galsim.Shear(e1=e1, e2=e2)
        gal.applyShear(shear)
        true_g1[obj_num] = shear.g1
        true_g2[obj_num] = shear.g2
        true_hlr[obj_num] = hlr

    # Figure out in which images this object was observed
    ncutout = meds_info['ncutout']
    file_id = meds_info['file_id']  # An array of indices

    if False:
        # Check that the MEDS treatment of the WCS is consistent with GalSim's.
        # MEDS uses X_IMAGE for the coadd image, and ALPHAMODEL, DELTAMODEL for the rest.
        # Also, the meds orig_row, orig_col use 0-based convention, the wcs uses 1-based.
        # The values should be accurate to around 1.e-3, so check to 1.e-2.
        x0 = coadd_info['X_IMAGE']
        y0 = coadd_info['Y_IMAGE']
        orig_row = meds_info['orig_row']
        orig_col = meds_info['orig_col']
        if do_print: print 'x0, y0 = ',x0,y0
        if do_print: print 'orig[0] = ',orig_col[0],orig_row[0]
        assert abs(x0 - (orig_col[0]+1)) < 1.e-2
        assert abs(y0 - (orig_row[0]+1)) < 1.e-2
        ra0 = coadd_info['ALPHAMODEL_J2000'] * galsim.degrees
        dec0 = coadd_info['DELTAMODEL_J2000'] * galsim.degrees
        world_pos0 = galsim.CelestialCoord(ra0,dec0)
        if do_print: print 'world_pos0 = ',ra0,dec0,world_pos0
        for k in range(1,ncutout):
            if do_print: print 'k = ',k
            id = file_id[k]
            if do_print: print 'id = ',id
            pos_k = image_wcs[id].toImage(world_pos0)
            if do_print: print 'pos_k = ',pos_k
            if do_print: print 'orign[',k,'] = ',orig_col[k],orig_row[k]
            assert abs(pos_k.x - (orig_col[k]+1)) < 1.e-2
            assert abs(pos_k.y - (orig_row[k]+1)) < 1.e-2
        # Despite all that, we will use the WIN values read in above from here on out.

    # Draw the object on each image
    for id in (id for id in file_id if id >= 0):
        assert id < nimages
        if do_print: print 'file_id = ',id,'  ',image_path[id]
        # The image onto which we will draw the object
        im = images[id]

        # Figure out the position to draw
        image_pos = im.wcs.toImage(world_pos)
        ix = int(math.floor(image_pos.x + 0.5))
        iy = int(math.floor(image_pos.y + 0.5))
        offset = image_pos - galsim.PositionD(ix,iy)
        if do_print: print 'image_pos = ',image_pos
        if do_print: print 'ix,iy,offset = ',ix,iy,offset

        # Build the PSF
        p = psf_params[id]
        if do_print: print 'p = ',p
        # Get the normalized position for the Chebyshev polynomials
        x = 2. * (image_pos.x - im.bounds.xmin) / (im.bounds.xmax - im.bounds.xmin) - 1.
        y = 2. * (image_pos.y - im.bounds.ymin) / (im.bounds.ymax - im.bounds.ymin) - 1.
        if do_print: print 'x,y = ',x,y
        f_size = p[0] + p[1]*x + p[2]*y + p[3]*(2.*x*x-1.) + p[4]*x*y + p[5]*(2.*y*y-1.)
        f_e1 = p[6] + p[7]*x + p[8]*y + p[9]*(2.*x*x-1.) + p[10]*x*y + p[11]*(2.*y*y-1.)
        f_e2 = p[12] + p[13]*x + p[14]*y + p[15]*(2.*x*x-1.) + p[16]*x*y + p[17]*(2.*y*y-1.)

        if do_print: print 'f = ',f_size,f_e1,f_e2
        fwhm = 0.9 + 0.1 * f_size  # ranges from 0.8 to 1.0
        e1 = 0.05 * f_e1           # ranges from -0.05 to 0.05
        e2 = 0.05 * f_e2           # ranges from -0.05 to 0.05
        # Well, the final range is not actuall restricted to these ranges.  Each term in the 
        # f's is limited to (-1,+1), but we add up a number of them, so the sum can be well
        # outside of this range.  But it's probably good enough.  The following calculation
        # shows empircally that about 2/3 of the f values fall within (-1,+1).
        if (-1 < f_size < 1): nf_pm1 += 1
        if (-1 < f_e1 < 1): nf_pm1 += 1
        if (-1 < f_e2 < 1): nf_pm1 += 1
        nf_tot += 3

        if do_print: print 'fwhm = ',fwhm
        if do_print: print 'e1,e2 = ',e1,e2
        psf = galsim.Gaussian(fwhm = fwhm)
        psf_shear = galsim.Shear(e1=e1, e2=e2)
        psf.applyShear(psf_shear)

        # Build the pixel
        pix = im.wcs.toWorld(galsim.Pixel(1.0), image_pos=image_pos)

        # Build the final object
        if is_star[obj_num]:
            psf.setFlux(float(flux))
            final = galsim.Convolve([psf, pix])
        else:
            final = galsim.Convolve([psf, pix, gal])

        # We want to adjust the flux for the relative zero point of the coadd catalog and
        # the single epoch catalog. The zeropoint is the magnitude that will produce 1 ADU.
        # So if the se zero point is larger than the coadd zero point, the object will appear
        # brighter (in ADU) on the single epoch image.  Thus, the flux rescaling is:
        final *= 10**( (magzp[id] - magzp_ref) / 2.5 )

        # Draw the object
        local_wcs = im.wcs.local(image_pos)
        if do_print: print 'local_wcs = ',local_wcs
        stamp = final.draw(wcs=local_wcs, use_true_center=False, offset=offset)
        stamp.setCenter(ix,iy)
        bounds = stamp.bounds & im.bounds
        if not bounds.isDefined(): 
            print 'Skipping object on coadd because stamp bounds did not overlap image bounds'
            print 'stamp.bounds = ',stamp.bounds
            print 'im.bounds = ',im.bounds
            print 'bounds = ',bounds
            continue
        if do_print:
            print 'drew stamp'
            print 'stamp.bounds = ',stamp.bounds
            print 'im.bounds = ',im.bounds
            print 'bounds = ',bounds
            print 'stamp.flux = ',stamp.array.sum()
            print 'stamp.max = ',stamp.array.max()
        im[bounds] += stamp[bounds]

        # Only add to the averages if we actually added the image.  So do this after the 
        # above if not bounds.isDefined() line.
        # Also, only do it for the single-eposures.  Not the coadd.
        if id > 0:
            psf_fwhm[obj_num] += fwhm
            psf_g1[obj_num] += psf_shear.g1
            psf_g2[obj_num] += psf_shear.g2
            wcs_decomp = local_wcs.getDecomposition()  # scale, shear, theta, flip
            if do_print: print 'wcs_decomp = ',wcs_decomp
            wcs_scale[obj_num] += wcs_decomp[0]
            wcs_g1[obj_num] += wcs_decomp[1].g1
            wcs_g2[obj_num] += wcs_decomp[1].g2
            wcs_theta[obj_num] += wcs_decomp[2].rad()
            # flip is always true
            assert wcs_decomp[3] == True
            nexp[obj_num] += 1

    if nexp[obj_num] > 0:
        psf_fwhm[obj_num] /= nexp[obj_num]
        psf_g1[obj_num] /= nexp[obj_num]
        psf_g2[obj_num] /= nexp[obj_num]
        wcs_scale[obj_num] /= nexp[obj_num]
        wcs_g1[obj_num] /= nexp[obj_num]
        wcs_g2[obj_num] /= nexp[obj_num]
        wcs_theta[obj_num] /= nexp[obj_num]
    
print 'fraction of f values between -1 and 1 = ',nf_pm1,'/',nf_tot,'=',float(nf_pm1)/nf_tot

# We will add a little bit of noise to the images.
print 'min_flux = ',min_flux
noise_sigma = min_flux/1000.
print 'add noise with sigma = ',noise_sigma
noise = galsim.GaussianNoise(ud, sigma = noise_sigma)

# For the coadd image, we just need to add noise and write the file to disk.
coadd_im = images[0]
coadd_im.addNoise(noise)
print 'Added noise to coadd image'
coadd_file = image_path[0]
print 'Original coadd file = ',coadd_file

# We will build a new hdulist for the new file and copy what we need from the old one.
# Also, we write this in uncompressed form and then fpack it to make sure that the 
# final result is funpack-able.
hdu_list = pyfits.open(coadd_file)
new_hdu_list = pyfits.HDUList()
# Copy the primary hdu
#new_hdu_list.append(hdu_list[0])

assert coadd_hdu == 1
coadd_im.write(hdu_list=new_hdu_list)
# copy over the header item SEXMGZPT
new_hdu_list[0].header['SEXMGZPT'] = hdu_list[coadd_hdu].header['SEXMGZPT']

# Next is the weight image
assert coadd_wt_hdu == 2
coadd_wt_im = galsim.fits.read(hdu_list=hdu_list[coadd_wt_hdu], compression='rice')
coadd_wt_im *= (1./noise_sigma**2) / coadd_wt_im.array.mean()
print 'coadd_wt_im.mean = ',coadd_wt_im.array.mean(),' should = ',1./noise_sigma**2
coadd_wt_im.write(hdu_list=new_hdu_list)

out_coadd_file = os.path.join(out_dir,os.path.basename(coadd_file))
print 'out_coadd_file = ',out_coadd_file
assert out_coadd_file.endswith('.fz')
if os.path.isfile(out_coadd_file): os.remove(out_coadd_file)
out_coadd_file = out_coadd_file[:-3]
if os.path.isfile(out_coadd_file): os.remove(out_coadd_file)
new_hdu_list.writeto(out_coadd_file, clobber=True)
print 'Wrote output coadd_file ',out_coadd_file

# Run fpack on the file
subprocess.Popen(['fpack','-D','-Y',out_coadd_file], close_fds=True).communicate()

# We will use the bounds of the coadd image below...
# GalSim doesn't currently convert easily between BoundsI and BoundsD...
coadd_bounds = galsim.BoundsD(coadd_im.bounds.xmin, coadd_im.bounds.xmax,
                              coadd_im.bounds.ymin, coadd_im.bounds.ymax)
# I add a slight border to make sure we don't draw things near the border twice. It seems 
# like that would be more problematic than missing some objects near the border.
coadd_bounds = coadd_bounds.addBorder(20.)
print 'using coadd_bounds = ',coadd_bounds

# Write the truth file
out_truth_file = os.path.join(out_dir,'end2end-truth.fits')
columns = [ pyfits.Column(name='id', format='J', array=coadd_cat['NUMBER']),
            pyfits.Column(name='flags', format='J', array=flags),
            pyfits.Column(name='is_star', format='I', array=is_star),
            pyfits.Column(name='true_g1', format='D', array=true_g1),
            pyfits.Column(name='true_g2', format='D', array=true_g2),
            pyfits.Column(name='true_hlr', format='D', array=true_hlr),
            pyfits.Column(name='ra', format='D', array=coadd_cat['ALPHAWIN_J2000']),
            pyfits.Column(name='dec', format='D', array=coadd_cat['DELTAWIN_J2000']),
            pyfits.Column(name='mag', format='D', array=coadd_cat['MAG_AUTO']),
            pyfits.Column(name='flux', format='D', array=coadd_cat['FLUX_AUTO']),
            pyfits.Column(name='nexp', format='I', array=nexp),
            pyfits.Column(name='mean_psf_g1', format='D', array=psf_g1),
            pyfits.Column(name='mean_psf_g2', format='D', array=psf_g2),
            pyfits.Column(name='mean_psf_fwhm', format='D', array=psf_fwhm),
            pyfits.Column(name='mean_wcs_g1', format='D', array=wcs_g1),
            pyfits.Column(name='mean_wcs_g2', format='D', array=wcs_g2),
            pyfits.Column(name='mean_wcs_scale', format='D', array=wcs_scale),
            pyfits.Column(name='mean_wcs_theta', format='D', array=wcs_theta),
            ]
coldefs = pyfits.ColDefs(columns)
table = pyfits.new_table(coldefs)
table.writeto(out_truth_file, clobber=True)

# Also, it turns out to be useful to have the coadd catalog in this directory as well:
# Update: now we make our own, which means the truth catalog isn't correct anymore (need to 
# match up objects by ra,dec).
#shutil.copy2(coaddcat_file, out_dir)

# Do the final processing on each single epoch image and write them to disk.
for image_num in range(1,nimages):
    im = images[image_num]
    file = image_path[image_num]
    print 'Finalize ',file

    # Get the catalog name.
    cat_file = file.replace('.fits.fz','_cat.fits')
    cat = pyfits.open(cat_file)[2].data  # 2 not 1!  hdu 1 has a bunch of meta data.
    if do_print: print 'catalog = ',cat_file
    print '   nobj in cat = ',len(cat)

    # Draw the objects in each image that weren't part of the coadd.
    # PSFEx will need these objects to get a good model of the PSF.
    nadded = 0
    for obj in cat:
        flags = obj['FLAGS']
        mag = obj['MAG_AUTO']
        flux = obj['FLUX_AUTO']
        x = obj['XWIN_IMAGE']
        y = obj['YWIN_IMAGE']
        spread = obj['SPREAD_MODEL']
        is_star = spread < 0.003

        if flags != 0: continue
        if mag > 24 or flux < 0: continue
        if ixx <= 0 or iyy <= 0 or hlr <= 0: continue

        image_pos = galsim.PositionD(x,y)
        world_pos = im.wcs.toWorld(image_pos)
        coadd_pos = coadd_im.wcs.toImage(world_pos)
        if do_print: print 'positions = ',image_pos,world_pos,coadd_pos
        if coadd_bounds.includes(coadd_pos):
            # Then we should have already drawn this object.
            if do_print: print 'coadd_pos is in coadd_bounds'
            continue

        if not is_star:
            ixx = obj['X2WIN_IMAGE']
            ixy = obj['XYWIN_IMAGE']
            iyy = obj['Y2WIN_IMAGE']
            hlr = obj['FLUX_RADIUS'] * 1.18 * 0.26
            gal = galsim.Exponential(half_light_radius = float(hlr), flux = float(flux))
            e1 = (ixx-iyy)/(ixx+iyy)
            e2 = 2.*ixy/(ixx+iyy)
            if (e1*e1 + e2*e2 > 0.8**2): continue
            gal.applyShear(e1=e1, e2=e2)

        ix = int(math.floor(x + 0.5))
        iy = int(math.floor(y + 0.5))
        offset = image_pos - galsim.PositionD(ix,iy)

        # Build the PSF the same way we did above
        p = psf_params[id]
        x = 2. * (image_pos.x - im.bounds.xmin) / (im.bounds.xmax - im.bounds.xmin) - 1.
        y = 2. * (image_pos.y - im.bounds.ymin) / (im.bounds.ymax - im.bounds.ymin) - 1.
        f_size = p[0] + p[1]*x + p[2]*y + p[3]*(2.*x*x-1.) + p[4]*x*y + p[5]*(2.*y*y-1.)
        f_e1 = p[6] + p[7]*x + p[8]*y + p[9]*(2.*x*x-1.) + p[10]*x*y + p[11]*(2.*y*y-1.)
        f_e2 = p[12] + p[13]*x + p[14]*y + p[15]*(2.*x*x-1.) + p[16]*x*y + p[17]*(2.*y*y-1.)

        fwhm = 0.9 + 0.1 * f_size
        e1 = 0.05 * f_e1
        e2 = 0.05 * f_e2
        psf = galsim.Gaussian(fwhm = fwhm)
        psf.applyShear(e1=e1,e2=e2)

        # Build the pixel
        pix = im.wcs.toWorld(galsim.Pixel(1.0), image_pos=image_pos)

        # Build the final object
        if is_star:
            final = galsim.Convolve([psf, pix])
        else:
            final = galsim.Convolve([psf, pix, gal])

        local_wcs = im.wcs.local(galsim.PositionD(x,y))
        stamp = final.draw(wcs=local_wcs, use_true_center=False, offset=offset)
        stamp.setCenter(ix,iy)

        bounds = stamp.bounds & im.bounds
        if not bounds.isDefined(): 
            continue
        im[bounds] += stamp[bounds]

        nadded += 1

    print '   Added ',nadded,' more objects from the single epoch catalog.'

    # Add the noise
    im.addNoise(noise)
    if do_print: print 'Added noise'

    # Add in the sky level, which is estimated in the background map by sextractor
    sky_im = galsim.fits.read(sky_path[image_num])
    im += sky_im

    # We will build a new hdulist for the new file and copy what we need from the old one.
    # Also, we write this in uncompressed form and then fpack it to make sure that the 
    # final result is funpack-able.
    hdu_list = pyfits.open(file)
    if do_print: 
        print 'hdu_list = ',hdu_list
        for h in hdu_list:
            print 'hdu = ',h.data
    new_hdu_list = pyfits.HDUList()

    # First is the image hdu.  We use the new image that we built.
    assert se_hdu==1
    if do_print: print 'im = ',im
    if do_print: print im.array
    im.write(hdu_list=new_hdu_list)
    if do_print: print 'after write im: ',new_hdu_list
    if do_print: print new_hdu_list[0].header
    if do_print: print new_hdu_list[0].data

    # Leave the badpix image the same.
    # TODO: It might be nice to add in artifacts in the image based on the bad pixel map.
    assert se_badpix_hdu==2
    badpix_im = galsim.fits.read(hdu_list=hdu_list[se_badpix_hdu], compression='rice')
    badpix_im = galsim.ImageS(badpix_im)
    if do_print: print 'badpix_im = ',badpix_im
    if do_print: print badpix_im.array
    if do_print: print badpix_im.array.astype(numpy.int32)
    if do_print: print 'max = ',badpix_im.array.max()
    badpix_im.write(hdu_list=new_hdu_list)
    if do_print: print 'after write badpix_im: ',new_hdu_list
    if do_print: print new_hdu_list[1].header
    if do_print: print new_hdu_list[1].data

    # Rescale the weight image to have the correct mean noise level.  We still let the
    # weight map be variable, but the nosie we add is constant.  (We can change this is 
    # we want using galsim.VariableGaussianNoise.)  However, it was requested that the 
    # mean level be accurate.  So we rescale the map to have the right mean.
    assert se_wt_hdu==3
    wt_im = galsim.fits.read(hdu_list=hdu_list[se_wt_hdu], compression='rice')
    wt_im *= (1./noise_sigma**2) / wt_im.array.mean()
    if do_print: print 'wt_im = ',wt_im
    if do_print: print wt_im.array
    wt_im.write(hdu_list=new_hdu_list)
    if do_print: print 'after write wt_im: ',new_hdu_list
    if do_print: print new_hdu_list[2].header
    if do_print: print new_hdu_list[2].data

    out_file = out_path[image_num]
    if do_print: print 'out_file = ',out_file
    if os.path.isfile(out_file): os.remove(out_file)
    assert out_file.endswith('.fz')
    out_file = out_file[:-3]
    if do_print: print 'out_file => ',out_file
    if os.path.isfile(out_file): os.remove(out_file)
    new_hdu_list.writeto(out_file, clobber=True)
    print '   Wrote file ',out_file

    # Run fpack on the file
    subprocess.Popen(['fpack','-D','-Y',out_file], close_fds=True).communicate()
    print '   Fpacked file ',out_file

    # Check that the file can be read correctly.
    f = pyfits.open(out_file + '.fz')
    if do_print: print 'f = ',f
    if do_print: print 'f[0] = ',f[0],f[0].data
    if do_print: print 'f[1] = ',f[1],f[1].data
    if do_print: print 'f[2] = ',f[2],f[2].data
    if do_print: print 'f[3] = ',f[3],f[3].data
    f[0].data  # Ignore the return value. Just check that it succeeds without raising an exception.
    f[1].data
    f[2].data
    f[3].data


print 'Done writing single-epoch files'
