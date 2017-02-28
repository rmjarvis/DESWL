from astropy.io import fits as pyfits
import numpy

imcat = pyfits.open('final_release/sva1_gold_r1.0_im3shape.fits')[1].data
ngcat = pyfits.open('final_release/sva1_gold_r1.0_ngmix.fits')[1].data
fcat = pyfits.open('final_release/sva1_gold_r1.0_wlinfo.fits')[1].data

sva1 = fcat['SVA1_FLAG'] == 0
ngmix = fcat['NGMIX_FLAG'] == 0
im3shape = fcat['IM3SHAPE_FLAG'] == 0

all_ng = sva1 & ngmix
all_im = sva1 & im3shape
ngmask = numpy.in1d(ngcat['COADD_OBJECTS_ID'], fcat['COADD_OBJECTS_ID'][all_ng])
immask = numpy.in1d(imcat['COADD_OBJECTS_ID'], fcat['COADD_OBJECTS_ID'][all_im])

ngcat = ngcat[ngmask]
imcat = imcat[immask]

nge1 = ngcat['E_1']
nge2 = ngcat['E_2']
ngvare = (ngcat['E_COV_1_1'] + ngcat['E_COV_2_2'])
ngw = ngcat['W']
ngs = ngcat['SENS_AVG']

sigsn_ng = numpy.sqrt(numpy.sum( (ngw**2 * (nge1**2 + nge2**2 - ngvare)) ) /
                      (2. * numpy.sum( (ngw**2 * ngs**2) )))
print
print 'ngmix has sigma_SN = ',sigsn_ng

ntot_ng = ((sigsn_ng**2 * numpy.sum(ngw * ngs)**2) /
           numpy.sum(ngw**2 * (ngs**2 * sigsn_ng**2 + ngvare/2.)))
neff_ng = ntot_ng / 139 / 60**2

print 'ntot = ',ntot_ng,' neff = ',neff_ng
nalt_ng = numpy.sum(ngw*ngs)**2 / numpy.sum(ngw**2 * ngs**2) * (1. -
            numpy.sum(ngw**2*ngvare/2.) / numpy.sum(ngw**2*(nge1**2+nge2**2)/2.))
print 'alt = ',nalt_ng,' neff = ',nalt_ng / 139 / 60**2
print

print 'mean shape = ',numpy.sum( (ngw * (nge1 + 1j*nge2)) ) / numpy.sum( ngw * ngs )
print 'sigma = ',numpy.sqrt(numpy.sum( (ngw**2 * (nge1**2 + nge2**2)) ) / numpy.sum( ngw * ngs )**2)
print

ime1 = imcat['E_1'] - imcat['NBC_C1']
ime2 = imcat['E_2'] - imcat['NBC_C2']
imw = imcat['w']
#imvare = 2. / imcat['snr']**2
# 1/w ~= 0.24**2 + vare
imvare = 1./imw - 0.24**2
imvare[imvare < 0.] = 0.
ims = imcat['NBC_M'] + 1.

print 'ngal with w > 0.24**-2 = ',numpy.sum(imw > 0.24**-2)
print 'mean w = ',numpy.mean(imw)
print 'median w = ',numpy.median(imw)
print 'max w = ',numpy.max(imw)
imw[imw > 0.24**-2] = 0.24**-2
print 'after set max to 0.24**-2:'
print 'mean w = ',numpy.mean(imw)
print 'median w = ',numpy.median(imw)
print 'max w = ',numpy.max(imw)

sigsn_im = numpy.sqrt(numpy.sum( (imw**2 * (ime1**2 + ime2**2 - imvare)) ) /
                      (2. * numpy.sum( (imw**2 * ims**2) )))
print
print 'im3shape has sigma_SN = ',sigsn_im

ntot_im = ((sigsn_im**2 * numpy.sum(imw * ims)**2) /
           numpy.sum(imw**2 * (ims**2 * sigsn_im**2 + imvare/2.)))
neff_im = ntot_im / 139 / 60**2

print 'ntot = ',ntot_im,' neff = ',neff_im
nalt_im = numpy.sum(imw*ims)**2 / numpy.sum(imw**2 * ims**2) * (1. -
            numpy.sum(imw**2*imvare/2.) / numpy.sum(imw**2*(ime1**2+ime2**2)/2.))
print 'alt = ',nalt_im,' neff = ',nalt_im / 139 / 60**2
print

print 'mean shape = ',numpy.sum( (imw * (ime1 + 1j*ime2)) ) / numpy.sum( imw )
print 'sigma = ',numpy.sqrt(numpy.sum( (imw**2 * (ime1**2 + ime2**2)) ) / numpy.sum( imw * ims )**2)

print
print 'CFHT definition:'
ntot_ng = numpy.sum(ngw)**2 / numpy.sum(ngw**2)
neff_ng = ntot_ng / 139 / 60**2
ntot_im = numpy.sum(imw)**2 / numpy.sum(imw**2)
neff_im = ntot_im / 139 / 60**2
print 'ngmix: ',ntot_ng, neff_ng
print 'im3shape: ',ntot_im, neff_im
sige_ng = numpy.sqrt(ntot_ng * numpy.sum(ngw**2*(nge1**2+nge2**2)/2)/numpy.sum(ngw*ngs)**2)
sige_im = numpy.sqrt(ntot_im * numpy.sum(imw**2*(ime1**2+ime2**2)/2)/numpy.sum(imw*ims)**2)
print 'ngmix sige = ',sige_ng
print 'im3shape sige = ',sige_im
