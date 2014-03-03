import galsim
from galsim import pyfits
import os
import math
import numpy

# Setup various file and directory names:
out_dir = '/direct/astro+astronfs03/workarea/mjarvis'

tile_name = 'DES0436-5748'
meds_dir = os.path.join('/astro/u/astrodat/data/DES/meds/011/20130820000021_DES0436-5748')
meds_file = tile_name + '-r-meds-011.fits.fz'

# Open the current meds file:
meds = pyfits.open(os.path.join(meds_dir,meds_file))

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

# Make arrays of the single-epoch images
image_cat = meds[2].data
image_path = image_cat['image_path']
sky_path = image_cat['sky_path']
seg_path = image_cat['seg_path']
magzp = image_cat['magzp']
scale = image_cat['scale']
nimages = len(image_cat)  # 138

# Read the image files and make blank ones.
# In principal, we could read the files and get the bounds from that, but 
# that's slow.  So we just use the known bounds and use that for everything.
se_bounds = galsim.BoundsI(1,2048,1,4096)
image_wcs = [ galsim.FitsWCS(file_name) for file_name in image_path ]
images = [ galsim.Image(wcs=wcs, bounds=se_bounds) for wcs in image_wcs ]
# Except the first one (the coadd image) is different
coadd_bounds = galsim.BoundsI(1,10000,1,10000)
images[0] = galsim.Image(wcs=image_wcs[0], bounds=coadd_bounds)

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

# Finally use the coadd_cat and the one in the meds file to draw each galaxy on the appropriate
# new image file.
meds_cat = meds[1].data

nobj = len(coadd_cat)  # 60505
assert nobj == len(meds_cat)
min_flux = 1.e100

for coadd_info, meds_info in zip(coadd_cat,meds_cat):
    print 'coadd id ',coadd_info['NUMBER']
    assert coadd_info['NUMBER'] == meds_info['number']

    # Skip anything with a flag.  This means either blending or saturation problems.
    flags = coadd_info['FLAGS']
    if flags != 0: 
        continue

    # For the position, use the WIN value.  Supposedly, this is the most accurate one.
    ra = coadd_info['ALPHAWIN_J2000'] * galsim.degrees
    dec = coadd_info['DELTAWIN_J2000'] * galsim.degrees
    world_pos = galsim.CelestialCoord(ra,dec)
    print 'world_pos = ',ra,dec,world_pos

    # Determine if this is a star
    spread = coadd_info['SPREAD_MODEL']
    is_star = spread < 0.003

    flux = coadd_info['FLUX_AUTO']
    if flux < min_flux: min_flux = flux

    if not is_star:
        # Get the parameters for building the galaxy
        ixx = coadd_info['X2WIN_IMAGE']
        ixy = coadd_info['XYWIN_IMAGE']
        iyy = coadd_info['Y2WIN_IMAGE']
        hlr = coadd_info['FLUX_RADIUS'] * 1.18 
        # This is approximate, since FLUX_RADIUS is based on a Gaussian sigma, not half-light 
        # radius.  But close enough.

        # Build the galaxy based on the coadd values
        print 'flux = ',flux,type(flux)
        print 'hlr = ',hlr,type(hlr)
        gal = galsim.Exponential(half_light_radius = float(hlr), flux = float(flux))
        e1 = (ixx-iyy)/(ixx+iyy)
        e2 = 2.*ixy/(ixx+iyy)
        print 'e1,e2 = ',e1,e2
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
        print 'x0, y0 = ',x0,y0
        print 'orig[0] = ',orig_col[0],orig_row[0]
        assert abs(x0 - (orig_col[0]+1)) < 1.e-2
        assert abs(y0 - (orig_row[0]+1)) < 1.e-2
        ra0 = coadd_info['ALPHAMODEL_J2000'] * galsim.degrees
        dec0 = coadd_info['DELTAMODEL_J2000'] * galsim.degrees
        world_pos0 = galsim.CelestialCoord(ra0,dec0)
        print 'world_pos0 = ',ra0,dec0,world_pos0
        for k in range(1,ncutout):
            print 'k = ',k
            id = file_id[k]
            print 'id = ',id
            pos_k = image_wcs[id].toImage(world_pos0)
            print 'pos_k = ',pos_k
            print 'orign[',k,'] = ',orig_col[k],orig_row[k]
            assert abs(pos_k.x - (orig_col[k]+1)) < 1.e-2
            assert abs(pos_k.y - (orig_row[k]+1)) < 1.e-2
        # Despite all that, we will use the WIN values read in above from here on out.

    # Draw the object on each image
    for id in (id for id in file_id if id >= 0):
        print 'id = ',id
        assert id < nimages
        # The image onto which we will draw the object
        im = images[id]

        # Figure out the position to draw
        image_pos = im.wcs.toImage(world_pos)
        ix = int(math.floor(image_pos.x + 0.5))
        iy = int(math.floor(image_pos.y + 0.5))
        offset = image_pos - galsim.PositionD(ix,iy)
        print 'image_pos = ',image_pos
        print 'ix,iy,offset = ',ix,iy,offset

        # Build the PSF
        p = psf_params[id]
        print 'p = ',p
        # Get the normalized position for the Chebyshev polynomials
        x = 2. * (image_pos.x - im.bounds.xmin) / (im.bounds.xmax - im.bounds.xmin) - 1.
        y = 2. * (image_pos.y - im.bounds.ymin) / (im.bounds.ymax - im.bounds.ymin) - 1.
        print 'x,y = ',x,y
        f_size = p[0] + p[1]*x + p[2]*y + p[3]*(2.*x*x-1.) + p[4]*x*y + p[5]*(2.*y*y-1.)
        f_e1 = p[6] + p[7]*x + p[8]*y + p[9]*(2.*x*x-1.) + p[10]*x*y + p[11]*(2.*y*y-1.)
        f_e2 = p[12] + p[13]*x + p[14]*y + p[15]*(2.*x*x-1.) + p[16]*x*y + p[17]*(2.*y*y-1.)
        print 'f = ',f_size,f_e1,f_e2
        fwhm = 0.9 + 0.1 * f_size  # ranges from 0.8 to 1.0
        e1 = 0.05 * f_e1           # ranges from -0.05 to 0.05
        e2 = 0.05 * f_e2           # ranges from -0.05 to 0.05
        # Well, the final range is not actuall restricted to these ranges.  Each term in the 
        # f's # is limited to (-1,+1), but we add up a number of them, so the sum can be well
        # outside of this range.  But it's probably good enough.
        print 'fwhm = ',fwhm
        print 'e1,e2 = ',e1,e2
        psf = galsim.Gaussian(fwhm = fwhm)
        psf.applyShear(e1=e1,e2=e2)

        # Build the pixel
        pix = im.wcs.toWorld(galsim.Pixel(1.0), image_pos=image_pos)

        # Build the final object
        if is_star:
            final = galsim.Convolve([psf, pix])
        else:
            final = galsim.Convolve([psf, pix, gal])

        # Draw the object
        local_wcs = im.wcs.local(image_pos)
        print 'local_wcs = ',local_wcs
        stamp = galsim.Image(256,256,wcs=local_wcs) # 256x256 is pretty big just to be safe...
        final.draw(stamp, use_true_center=False, offset=offset)
        stamp.setCenter(ix,iy)
        bounds = stamp.bounds & im.bounds
        if not bounds.isDefined(): continue
        im[bounds] += stamp[bounds]
    

print 'min_flux = ',min_flux
print 'add noise with sigma = ',min_flux/1000.
noise = galsim.GaussianNoise(ud, sigma = min_flux/1000.)

for id in range(nimages):
    im.addNoise(noise)
    print 'orig path = ',image_path[id]
    print 'des_root = ',des_root
    print 'out_dir = ',out_dir
    image_file = os.path.join(out_dir,os.path.basename(image_path[id]))
    print 'new path = ',image_file
    im.write(image_file)
