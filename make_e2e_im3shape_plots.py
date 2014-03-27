import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('e2e_im3shape-2.pdf')

truth = pyfits.open('end2end-truth.fits')[1]
im3shape = pyfits.open('end2end-im3shape-2.fits')[1]

tdata = truth.data
idata = im3shape.data
print 'tdata has %d columns, %d entries'%(len(tdata.columns), len(tdata))
print 'idata has %d columns, %d entries'%(len(idata.columns), len(idata))

# In v2, there are duplicates, so pull out just one copy of each entry.
im3shape_ids = im3shape.data['identifier']
if len(np.unique(im3shape_ids)) != len(im3shape_ids):
    unique_mask = np.unique(im3shape_ids, return_index=True)[1]
    idata = idata[unique_mask]
    print 'Removed duplicates from im3shape catalog.  Now has %d entries'%len(idata)

# Need to make a matched set.  Truth has all number from 1 to N, im3shape is missing some.
# This next line is super slow.  So we write the answer to disk to avoid running this line 
# every time.
if False:
    mask1 = [ i for i in range(len(tdata)) if i+1 in idata['identifier'] ]
    with open('mask1.pkl','wb') as fout:
        pickle.dump(mask1,fout)
    print 'Built mask1'
else:
    with open('mask1.pkl','rb') as fin:
        mask1 = pickle.load(fin)
    print 'Read mask1 from disk.'
tdata = tdata[mask1]

# Some of the items are flagged.  Remove these.
# Also select out just the galaxy objects.
mask2 = (tdata['flags'] == 0) & (tdata['is_star'] == 0) & (idata['flag'] == 0)

# Extract values that we want to plot
tg1 = tdata['true_g1'][mask2]
tg2 = tdata['true_g2'][mask2]
thlr = tdata['true_hlr'][mask2]
tid = tdata['id'][mask2]
tra = tdata['ra'][mask2]
tdec = tdata['dec'][mask2]
tflux = tdata['flux'][mask2]
tmag = tdata['mag'][mask2]

ie1 = idata['e1'][mask2]
ie2 = idata['e2'][mask2]
ir = idata['radius'][mask2]
iid = idata['identifier'][mask2]
ira = idata['ra'][mask2]
idec = idata['dec'][mask2]
iflux = idata['bulge_flux'][mask2] + idata['disc_flux'][mask2]
isnr = idata['snr'][mask2]
imag = -2.5*np.log10(iflux)
imag[iflux <= 0] = 99
imag2 = -2.5*np.log10(isnr)
imag2[isnr <= 0] = 99
irr = idata['radius_ratio'][mask2]
iba = idata['bulge_A'][mask2]
ida = idata['disc_A'][mask2]
ibi = idata['bulge_index'][mask2]
idi = idata['disc_index'][mask2]
ibf = idata['bulge_flux'][mask2]
idf = idata['disc_flux'][mask2]
ifr = idata['flux_ratio'][mask2]
ilike = idata['likelihood'][mask2]
ilnlike = np.log(np.abs(ilike))
immin = idata['model_min'][mask2]
immax = idata['model_max'][mask2]
ideb = idata['delta_e_bulge'][mask2]
idtb = idata['delta_theta_bulge'][mask2]
iraas = idata['ra_as'][mask2]
idecas = idata['dec_as'][mask2]
print 'Extracted all data fields'

m = 2.3

# For v2, this is a first-pass effort to color code "good" shear values.
good1 = (abs(m*tg1 - ie1) < 0.05) & (abs(m*tg2 - ie2) < 0.05)
good2 = (abs(m*tg1 - ie1) < 0.1) & (abs(m*tg2 - ie2) < 0.1) & ~good1
bad = (~good1) & (~good2)

def plt_scatter(xval, yval, xlabel, ylabel, mask=None, m=None):
    """Make a scatter plot of two values with labels.
    If mask is not None, it is applied to xval, yval, and the color-coding good1, good2, bad.
    If m is not None, the line y=mx is drawn.
    """
    plt.clf()
    if mask is not None:
        xval = xval[mask]
        yval = yval[mask]
        good1m = good1[mask]
        good2m = good2[mask]
        badm = bad[mask]
    else:
        good1m = good1
        good2m = good2
        badm = bad
    plt.grid()
    plt.scatter(xval[badm],yval[badm],s=0.4,rasterized=True,color='red')
    plt.scatter(xval[good2m],yval[good2m],s=0.4,rasterized=True,color='blue')
    plt.scatter(xval[good1m],yval[good1m],s=0.4,rasterized=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if m is not None:
        if m > 1.:
            mline, = plt.plot([-1./m,1./m],[-1.,1.],'c-')
        else:
            mline, = plt.plot([-1.,1.],[-m,m],'c-')
        plt.legend([mline], ['m = %f'%m])
    pp.savefig()
    print 'Plotted %s vs %s'%(ylabel, xlabel)

plt_scatter(tg1, ie1, 'True g1', 'im3shape e1', m=m)
plt_scatter(tg2, ie2, 'True g2', 'im3shape e2', m=m)
plt_scatter(tg2, ie1, 'True g2', 'im3shape e1')
plt_scatter(tg1, ie2, 'True g1', 'im3shape e2')
plt_scatter(thlr, ir, 'True hlr', 'im3shape radius')
plt_scatter(thlr, ir,'True hlr', 'im3shape radius', ir < 10)
plt_scatter(tflux, iflux, 'True flux', 'im3shape bulge_flux + disc_flux')
plt_scatter(tflux, iflux, 'True flux', 'im3shape bulge_flux + disc_flux', abs(iflux) < 1.e2)
plt_scatter(tflux, iflux, 'True flux', 'im3shape bulge_flux + disc_flux', abs(iflux) < 1.e1)
plt_scatter(tflux, isnr, 'True flux', 'im3shape snr')
plt_scatter(tflux, isnr, 'True flux', 'im3shape snr', abs(isnr) < 1.e5)
plt_scatter(tflux, isnr, 'True flux', 'im3shape snr', abs(isnr) < 1.e4)
plt_scatter(tmag, imag, 'True mag', '-2.5 log10(flux)', iflux > 0)
plt_scatter(tmag, imag2, 'True mag', '-2.5 log10(snr)', isnr > 0)

#plt_scatter(ir, irr, 'radius', 'radius_ratio')  # radius_ratio = 1
plt_scatter(iba, ida, 'bulge_A', 'disc_A')
plt_scatter(iba, ida, 'bulge_A', 'disc_A', (abs(iba) < 1.e3) & (abs(ida) < 1.e3))
plt_scatter(iba, ida, 'bulge_A', 'disc_A', (abs(iba) < 1.e2) & (abs(ida) < 1.e2))
#plt_scatter(ibi, idi, 'bulge_index', 'disc_index')  # bulge_index = 4, disc_index = 1
plt_scatter(ibf, idf, 'bulge_flux', 'disc_flux')
plt_scatter(ibf, idf, 'bulge_flux', 'disc_flux', (abs(ibf) < 1.e2) & (abs(idf) < 1.e2))
plt_scatter(ibf, idf, 'bulge_flux', 'disc_flux', (abs(ibf) < 1.e1) & (abs(idf) < 1.e1))
plt_scatter(iflux, ifr, 'total flux', 'flux_ratio')
plt_scatter(iflux, ifr, 'total flux', 'flux_ratio', (abs(ibf) < 1.e2) & (abs(idf) < 1.e2))
plt_scatter(iflux, ifr, 'total flux', 'flux_ratio', (abs(ibf) < 1.e1) & (abs(idf) < 1.e1))
plt_scatter(isnr, ilike, 'snr', 'likelihood')
plt_scatter(isnr, ilike, 'snr', 'likelihood', abs(ilike) < 1.e6)
plt_scatter(isnr, ilike, 'snr', 'likelihood', abs(ilike) < 1.e5)
plt_scatter(imag2, ilnlike, '-2.5 log10(snr)', 'ln(abs(likelihood))', isnr > 0)
plt_scatter(immin, immax, 'model_min', 'model_max')
#plt_scatter(ideb, idtb, 'delta_e_bulge', 'delta_theta_bulge')  # both = 0
plt_scatter(iraas, idecas, 'ra_as', 'dec_as')
plt_scatter(iraas, idecas, 'ra_as', 'dec_as', (abs(iraas) < 1) & (abs(idecas) < 1))

good1 = (iba > 10) & (ida > 0.05)
good2 = (iba > 0) & (ida > 0) & ~good1
bad = (~good1) & (~good2)

plt_scatter(iba, ida, 'bulge_A', 'disc_A', (abs(iba) < 1.e2) & (abs(ida) < 1.e2))
plt_scatter(iba, ida, 'bulge_A', 'disc_A', (abs(iba) < 1.e3) & (abs(ida) < 1.e3))
plt_scatter(tg1, ie1, 'True g1', 'im3shape e1')
plt_scatter(tg2, ie2, 'True g2', 'im3shape e2')
plt_scatter(thlr, ir,'True hlr', 'im3shape radius', ir < 10)
plt_scatter(tg1, ie1, 'True g1', 'im3shape e1',good1, m=m)
plt_scatter(tg2, ie2, 'True g2', 'im3shape e2',good1, m=m)
plt_scatter(thlr, ir,'True hlr', 'im3shape radius', good1 & (ir<10))
plt_scatter(ibf, idf, 'bulge_flux', 'disc_flux', (abs(ibf) < 1.e2) & (abs(idf) < 1.e2))
plt_scatter(ibf, idf, 'bulge_flux', 'disc_flux', (abs(ibf) < 1.e1) & (abs(idf) < 1.e1))
plt_scatter(isnr, ilike, 'snr', 'likelihood', abs(ilike) < 1.e5)
plt_scatter(imag2, ilnlike, '-2.5 log10(snr)', 'ln(abs(likelihood))', isnr > 0)

pp.close()
