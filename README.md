
This repository serves several purposes.

1. First, the issues page on GitHub can serve as a general clearing house for any
   issue related to the WL pipeline, tests, the catalogs, etc.  If you have any
   issue that you would like to bring to the attention of the WL pipeline testing
   subgroup, you may post it here.  

   We may repost the issue in a different location if there is a more appropriate
   place, or we may deal with it here.  But as a user, you do not need to worry about
   finding the right place to post any problem you are having.  Posting it here
   is perfectly fine.

   If there is any sensitive information or plots that you need to reference, we ask
   that you post them to redmine wiki

   https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Shear_Pipeline_Development_and_Testing

   and link to them from here.  But since none of us find the redmide issues page 
   terribly useful, we are encouraging the use of this GitHub page for all issues
   reporting related the WL pipeline.


2. Second, this repository does host some code that is related to the DES WL pipeline
   or testing thereof that does not have some other, more appropriate home.

    Items in this repository:

    - e2e

      This directory contains the code to make the so-called "end-to-end" test images.
      The script generates low noise images corresponding to a single coadd tile,
      including all the corresponding single-epoch images.  It is intended to be used
      to test that the weak lensing pipeline gets all the various sign conventions
      correct for the wcs, the PSFEx interpolation, etc.

        - The starting point is an existing MEDS file from actual data.
        - It uses the existing data images for the WCS, although there are plans to 
            modify this to use either more extreme or simpler WCS functions to make the 
            tests either.
        - It uses an elliptical Gaussian PSF.  The variation across each chip uses 
            second order polynomials for the size, e1, and e2.
        - It adds a small amount of noise, set to be the minimum flux / 1000.

    - psfex

      This directory contains the code to rerun PSFEx on the red images using custom
      parameters and star selection appropriate for the requirements of the weak
      lensing pipeline.  In particular:

        - It uses our own star-galaxy separation to select more reliable stars.
        - It removes bright stars that suffer from the bright-fatter problem.
        - It removes stars that fall in or next to the tape bumps to avoid these
            stars from contaminating the fits.*
        - It flags chips that are considered bad for any of the following reasons:*
            - Too few stars (<50)
            - Too many stars (>500)
            - Too large seeing (>1.8 arcsec FWHM)
            - Bad rho statistics (TBD)

      Note: Items above with the * are planned.  Not yet done.
          


Other repos of relevance:

- The im3shape shear estimation code:

   - https://bitbucket.org/joezuntz/im3shape

- The ngmix shear estimation code: 

   - https://github.com/esheldon/ngmix

   - https://github.com/esheldon/gmix_meds

- Code for making the meds files (also includes the shapelet PSF estimation code
  as well as a shapelet shear estimation code that we used to use but has too
  large noise biases):

   - https://github.com/rmjarvis/deswl_shapelets

- Code for computing shear correlation functions (i.e. the corr2 executable
  and also the treecorr python module):

   - https://github.com/rmjarvis/TreeCorr

- GalSim

   - https://github.com/GalSim-developers/GalSim

