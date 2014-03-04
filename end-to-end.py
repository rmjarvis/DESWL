import galsim
from galsim import pyfits
import os
import math
import numpy

# Setup various file and directory names:
work_dir = '/direct/astro+astronfs03/workarea/mjarvis'
tile_name = 'DES0436-5748'
out_dir = os.path.join(work_dir,tile_name)

meds_dir = os.path.join('/astro/u/astrodat/data/DES/meds/011/20130820000021_DES0436-5748')
meds_file = tile_name + '-r-meds-011.fits.fz'

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

flags = coadd_cat['FLAGS']
mag = coadd_cat['MAG_AUTO']
flux = coadd_cat['FLUX_AUTO']
ixx = coadd_cat['X2WIN_IMAGE']
ixy = coadd_cat['XYWIN_IMAGE']
iyy = coadd_cat['Y2WIN_IMAGE']
hlr = coadd_cat['FLUX_RADIUS'] * 1.18 # This is approximate, since FLUX_RADIUS is based on a 
                                      # Gaussian sigma, not half-light radius. But close enough.
e1 = (ixx-iyy)/(ixx+iyy)
e2 = 2.*ixy/(ixx+iyy)

print 'nobj with flags != 0 = ',(flags != 0).sum()
print 'nobj with mag > 24 = ',(mag > 24).sum()
print 'nobj with flux <= 0 = ',(flux <= 0).sum()
print 'nobj with ixx <= 0 = ',(ixx <= 0).sum()
print 'nobj with iyy <= 0 = ',(iyy <= 0).sum()
print 'nobj with hlr <= 0 = ',(hlr <= 0).sum()
print 'nobj with hlr > 8 = ',(hlr > 8).sum()
print 'nobj with e > 0.8 = ',(e1*e1 + e2*e2 > 0.8**2).sum()

mask = ( (flags != 0) |                         # Ignore anything with an input flag
         (mag > 24) | (flux <= 0) |             # Ignore faint things, or those with negative flux
         (ixx <= 0) | (iyy <= 0) | (hlr <= 0) | # Ignore negative sizes
         (hlr > 8) |                            # Ignore very large objects
         (e1*e1 + e2*e2 > 0.8**2) )             # Ignore very high ellipticities

print 'nobj passing all masks = ',(mask == 0).sum()

nf_pm1 = 0
nf_tot = 0

for k in range(len(coadd_cat)):
    if mask[k]: continue
    coadd_info = coadd_cat[k]
    meds_info = meds_cat[k]

    print 'coadd id ',coadd_info['NUMBER']
    assert coadd_info['NUMBER'] == meds_info['number']

    # Should already have skipped objects with flags
    assert coadd_info['FLAGS'] == 0

    # For the position, use the WIN value.  Supposedly, this is the most accurate one.
    ra = coadd_info['ALPHAWIN_J2000'] * galsim.degrees
    dec = coadd_info['DELTAWIN_J2000'] * galsim.degrees
    world_pos = galsim.CelestialCoord(ra,dec)
    print 'world_pos = ',ra,dec,world_pos

    # Determine if this is a star
    spread = coadd_info['SPREAD_MODEL']
    is_star = spread < 0.003
    print 'spread = ',spread,' is_star? ',is_star

    print 'mag = ',coadd_info['MAG_AUTO']
    assert coadd_info['MAG_AUTO'] <= 24

    flux = coadd_info['FLUX_AUTO']
    print 'flux = ',flux,type(flux)
    assert flux > 0
    if flux < min_flux: 
        min_flux = flux

    if not is_star:
        # Get the parameters for building the galaxy
        ixx = coadd_info['X2WIN_IMAGE']
        ixy = coadd_info['XYWIN_IMAGE']
        iyy = coadd_info['Y2WIN_IMAGE']
        hlr = coadd_info['FLUX_RADIUS'] * 1.18 
        # This is approximate, since FLUX_RADIUS is based on a Gaussian sigma, not half-light 
        # radius.  But close enough.
        print 'hlr = ',hlr,type(hlr)
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
        print 'e1,e2 = ',e1,e2
        assert e1*e1 + e2*e2 <= 0.8
        gal.applyShear(e1=e1, e2=e2)

    # Figure out in which images this object was observed
    ncutout = meds_info['ncutout']
    file_id = meds_info['file_id']  # An array of indices

    if True:
        # Check that the MEDS treatment of the WCS is consistent with GalSim's.
        # MEDS uses X_IMAGE for the coadd image, and ALPHAMODEL, DELTAMODEL for the rest.
        # Also, the meds orig_row, orig_col use 0-based convention, the wcs uses 1-based.
        # The values should be accurate to around 1.e-3, so check to 1.e-2.
        x0 = coadd_info['X_IMAGE']
        y0 = coadd_info['Y_IMAGE']
        orig_row = meds_info['orig_row']
        orig_col = meds_info['orig_col']
        #print 'x0, y0 = ',x0,y0
        #print 'orig[0] = ',orig_col[0],orig_row[0]
        assert abs(x0 - (orig_col[0]+1)) < 1.e-2
        assert abs(y0 - (orig_row[0]+1)) < 1.e-2
        ra0 = coadd_info['ALPHAMODEL_J2000'] * galsim.degrees
        dec0 = coadd_info['DELTAMODEL_J2000'] * galsim.degrees
        world_pos0 = galsim.CelestialCoord(ra0,dec0)
        #print 'world_pos0 = ',ra0,dec0,world_pos0
        for k in range(1,ncutout):
            #print 'k = ',k
            id = file_id[k]
            #print 'id = ',id
            pos_k = image_wcs[id].toImage(world_pos0)
            #print 'pos_k = ',pos_k
            #print 'orign[',k,'] = ',orig_col[k],orig_row[k]
            assert abs(pos_k.x - (orig_col[k]+1)) < 1.e-2
            assert abs(pos_k.y - (orig_row[k]+1)) < 1.e-2
        # Despite all that, we will use the WIN values read in above from here on out.

    # Draw the object on each image
    for id in (id for id in file_id if id >= 0):
        assert id < nimages
        print 'id = ',id,'  ',image_path[id]
        # The image onto which we will draw the object
        im = images[id]

        # Figure out the position to draw
        image_pos = im.wcs.toImage(world_pos)
        ix = int(math.floor(image_pos.x + 0.5))
        iy = int(math.floor(image_pos.y + 0.5))
        offset = image_pos - galsim.PositionD(ix,iy)
        print 'image_pos = ',image_pos
        #print 'ix,iy,offset = ',ix,iy,offset

        # Build the PSF
        p = psf_params[id]
        #print 'p = ',p
        # Get the normalized position for the Chebyshev polynomials
        x = 2. * (image_pos.x - im.bounds.xmin) / (im.bounds.xmax - im.bounds.xmin) - 1.
        y = 2. * (image_pos.y - im.bounds.ymin) / (im.bounds.ymax - im.bounds.ymin) - 1.
        #print 'x,y = ',x,y
        f_size = p[0] + p[1]*x + p[2]*y + p[3]*(2.*x*x-1.) + p[4]*x*y + p[5]*(2.*y*y-1.)
        f_e1 = p[6] + p[7]*x + p[8]*y + p[9]*(2.*x*x-1.) + p[10]*x*y + p[11]*(2.*y*y-1.)
        f_e2 = p[12] + p[13]*x + p[14]*y + p[15]*(2.*x*x-1.) + p[16]*x*y + p[17]*(2.*y*y-1.)

        #print 'f = ',f_size,f_e1,f_e2
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

        #print 'fwhm = ',fwhm
        #print 'e1,e2 = ',e1,e2
        psf = galsim.Gaussian(fwhm = fwhm)
        psf.applyShear(e1=e1,e2=e2)

        # Build the pixel
        pix = im.wcs.toWorld(galsim.Pixel(1.0), image_pos=image_pos)

        # Build the final object
        if is_star:
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
        #print 'local_wcs = ',local_wcs
        stamp = final.draw(wcs=local_wcs, use_true_center=False, offset=offset)
        stamp.setCenter(ix,iy)
        bounds = stamp.bounds & im.bounds
        if not bounds.isDefined(): continue
        im[bounds] += stamp[bounds]
    
print 'fraction of f values between -1 and 1 = ',nf_pm1,'/',nf_tot,'=',float(nf_pm1)/nf_tot

# We will add a little bit of noise to the images.
print 'min_flux = ',min_flux
print 'add noise with sigma = ',min_flux/1000.
noise = galsim.GaussianNoise(ud, sigma = min_flux/1000.)

# For the coadd image, we just need to add noise and write the file to disk.
coadd_im = images[0]
coadd_im.addNoise(noise)
print 'Added noise to coadd image'
coadd_file = image_path[0]
print 'Original coadd file = ',coadd_file
hdu_list = pyfits.open(coadd_file)
new_hdu = pyfits.HDUList()
coadd_im.write(hdu_list=new_hdu, compression='rice')
hdu_list[coadd_hdu] = new_hdu[1]
out_coadd_file = os.path.join(out_dir,os.path.basename(coadd_file))
coadd_im.write(out_coadd_file)
print 'Wrote output coadd_file ',out_coadd_file
# GalSim doesn't currently convert easily between BoundsI and BoundsD...
coadd_bounds = galsim.BoundsD(coadd_im.bounds.xmin, coadd_im.bounds.xmax,
                              coadd_im.bounds.ymin, coadd_im.bounds.ymax)
# I add a slight border to make sure we don't draw things near the border twice. It seems 
# like that would be more problematic than missing some objects near the border.
coadd_bounds = coadd_bounds.addBorder(20.)
print 'using coadd_bounds = ',coadd_bounds

# Do the final processing on each single epoch image and write them to disk.
for k in range(1,nimages):
    im = images[k]
    file = image_path[k]
    print 'Finalize ',file

    # Get the catalog name.
    cat_file = file.replace('.fits.fz','_cat.fits')
    cat = pyfits.open(cat_file)[2].data  # 2 not 1!  hdu 1 has a bunch of meta data.
    print 'catalog = ',cat_file
    print 'nobj in cat = ',len(cat)

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
        #print 'positions = ',image_pos,world_pos,coadd_pos
        if coadd_bounds.includes(coadd_pos):
            # Then we should have already drawn this object.
            #print 'coadd_pos is in coadd_bounds'
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
        if not bounds.isDefined(): continue
        im[bounds] += stamp[bounds]

        nadded += 1

    print 'Added ',nadded,' more objects from the single epoch catalog.'

    # Add the noise
    im.addNoise(noise)
    print 'Added noise'

    # Write the image to disk
    # We keep the original weight and badpix images, and just replace the actual image.
    print 'des_root = ',des_root
    print 'out_dir = ',out_dir
    hdu_list = pyfits.open(file)
    new_hdu = pyfits.HDUList()
    im.write(hdu_list=new_hdu, compression='rice')
    hdu_list[se_hdu] = new_hdu[1]
    out_file = os.path.join(out_dir,os.path.basename(file))
    im.write(out_file)
    print 'Wrote file ',out_file


