Items in this repository:

- end-to-end.py

  A script to generate low noise images corresponding to a single coadd tile,
  including all the corresponding single-epoch images.  It is intended to be 
  used to test that the weak lensing pipeline gets all the various sign 
  conventions correct for the wcs, the PSFEx interpolation, etc.

  - The starting point is an existing MEDS file from actual data.
  - It uses the existing data images for the WCS, although there are plans to 
    modify this to use either more extreme or simpler WCS functions to make the 
    tests either.
  - It uses an elliptical Gaussian PSF.  The variation across each chip uses 
    second order polynomials for the size, e1, and e2.
  - It adds a small amount of noise, set to be the minimum flux / 1000.
