This repository has some scripts that we use as part of to the DES WL pipeline
as well as scripts used for testing the pipeline.

In addition, the issues page may be used as a common location for any issues
people find while testing the DES shear catalogs.  If there is any sensitive
information or plots that you need to reference, we ask that you post them to
redmine wiki

https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Shear_Pipeline_Development_and_Testing

and link to them from here.  But since none of us find the redmide issues page 
terribly usable, we are encouraging the use of this GitHub page for all issues
reporting related the WL pipeline.



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
