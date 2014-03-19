import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('e2e_im3shape.pdf')

truth = pyfits.open('end2end-truth.fits')[1]
im3shape = pyfits.open('end2end-im3shape-1.fits')[1]

# Need to make a matched set.  Truth has all number from 1 to N, im3shape is missing some.
# This next line is super slow.  So we write the answer to disk to avoid running this line 
# every time.
if False:
    mask1 = [ i for i in range(len(truth.data)) if i+1 in im3shape.data['identifier'] ]
    with open('mask1.pkl','wb') as fout:
        pickle.dump(mask1,fout)
else:
    with open('mask1.pkl','rb') as fin:
        mask1 = pickle.load(fin)

t1 = truth.data[mask1]
i1 = im3shape.data

# Some of the items are flagged.  Remove these.
# Also select out just the galaxy objects.
mask2 = (t1['flags'] == 0) & (t1['is_star'] == 0) & (i1['flag'] == 0)

# Plot true g vs measured g
tg1 = t1['true_g1'][mask2]
tg2 = t1['true_g2'][mask2]
thlr = t1['true_hlr'][mask2]
tid = t1['id'][mask2]
tra = t1['ra'][mask2]
tdec = t1['dec'][mask2]
tflux = t1['flux'][mask2]
tmag = t1['mag'][mask2]

ie1 = i1['e1'][mask2]
ie2 = i1['e2'][mask2]
ir = i1['radius'][mask2]
iid = i1['identifier'][mask2]
ira = i1['ra'][mask2]
idec = i1['dec'][mask2]
iflux = i1['bulge_flux'][mask2] + i1['disc_flux'][mask2]
isnr = i1['snr'][mask2]
imag = -2.5*np.log10(iflux[iflux > 0])
imag2 = -2.5*np.log10(isnr[isnr > 0])

plt.clf()
plt.scatter(tg1,ie1,rasterized=True)
plt.xlabel('True g1')
plt.ylabel('im3shape e1')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tg2,ie2,rasterized=True)
plt.xlabel('True g2')
plt.ylabel('im3shape e2')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tg1,ie2,rasterized=True)
plt.xlabel('True g1')
plt.ylabel('im3shape e2')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tg2,ie1,rasterized=True)
plt.xlabel('True g2')
plt.ylabel('im3shape e1')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(thlr,ir,rasterized=True)
plt.xlabel('True hlr')
plt.ylabel('im3shape radius')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(thlr[ir<10],ir[ir<10],rasterized=True)
plt.xlabel('True hlr')
plt.ylabel('im3shape radius')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tflux,iflux,rasterized=True)
plt.xlabel('True flux')
plt.ylabel('im3shape bulge_flux + disc_flux')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tflux,isnr,rasterized=True)
plt.xlabel('True flux')
plt.ylabel('im3shape snr')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tmag[iflux>0],imag,rasterized=True)
plt.xlabel('True mag')
plt.ylabel('-2.5log10(flux)')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tmag[isnr>0],imag2,rasterized=True)
plt.xlabel('True mag')
plt.ylabel('-2.5log10(snr)')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tra,ira,rasterized=True)
plt.xlabel('RA from truth')
plt.ylabel('RA from im3shape')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tdec,idec,rasterized=True)
plt.xlabel('Dec from truth')
plt.ylabel('Dec from im3shape')
plt.grid()
pp.savefig()

plt.clf()
plt.scatter(tid,iid,rasterized=True)
plt.xlabel('ID from truth')
plt.ylabel('ID from im3shape')
plt.grid()
pp.savefig()

pp.close()
