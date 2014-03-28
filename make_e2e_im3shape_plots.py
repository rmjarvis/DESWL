import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('e2e_im3shape-2.pdf')

truth = pyfits.open('end2end-truth.fits')[1]
im3shape = pyfits.open('end2end-im3shape-2.fits')[1]

tdata = truth.data
xdata = im3shape.data
print 'tdata has %d columns, %d entries'%(len(tdata.columns), len(tdata))
print 'xdata has %d columns, %d entries'%(len(xdata.columns), len(xdata))

# In v2, there are duplicates, so pull out just one copy of each entry.
im3shape_ids = im3shape.data['identifier']
if len(np.unique(im3shape_ids)) != len(im3shape_ids):
    unique_mask = np.unique(im3shape_ids, return_index=True)[1]
    xdata = xdata[unique_mask]
    print 'Removed duplicates from im3shape catalog.  Now has %d entries'%len(xdata)

# Need to make a matched set.  Truth has all number from 1 to N, im3shape is missing some.
# This next line is super slow.  So we write the answer to disk to avoid running this line 
# every time.
mask_file = 'im3shape_mask1.pkl'
try:
    with open(mask_file,'rb') as fin:
        mask1 = pickle.load(fin)
    print 'Read mask1 from',mask_file
except:
    mask1 = [ i for i in range(len(tdata)) if i+1 in xdata['identifier'] ]
    print 'Built mask1'
    with open(mask_file,'wb') as fout:
        pickle.dump(mask1,fout)
    print 'Wrote mask to file:',mask_file

tdata = tdata[mask1]

# Some of the items are flagged.  Remove these.
# Also select out just the galaxy objects.
mask2 = (tdata['flags'] == 0) & (tdata['is_star'] == 0) & (xdata['flag'] == 0)

# Extract values that we want to plot
tid = tdata['id'][mask2]
tg1 = tdata['true_g1'][mask2]
tg2 = tdata['true_g2'][mask2]
thlr = tdata['true_hlr'][mask2]
tra = tdata['ra'][mask2]
tdec = tdata['dec'][mask2]
tflux = tdata['flux'][mask2]
tmag = tdata['mag'][mask2]

xid = xdata['identifier'][mask2]
xe1 = xdata['e1'][mask2]
xe2 = xdata['e2'][mask2]
xr = xdata['radius'][mask2]
xra = xdata['ra'][mask2]
xdec = xdata['dec'][mask2]
xflux = xdata['bulge_flux'][mask2] + xdata['disc_flux'][mask2]
xsnr = xdata['snr'][mask2]
xmag = -2.5*np.log10(xflux)
xmag[xflux <= 0] = 99
xmag2 = -2.5*np.log10(xsnr)
xmag2[xsnr <= 0] = 99
xrr = xdata['radius_ratio'][mask2]
xba = xdata['bulge_A'][mask2]
xda = xdata['disc_A'][mask2]
xbi = xdata['bulge_index'][mask2]
xdi = xdata['disc_index'][mask2]
xbf = xdata['bulge_flux'][mask2]
xdf = xdata['disc_flux'][mask2]
xfr = xdata['flux_ratio'][mask2]
xlike = xdata['likelihood'][mask2]
xlnlike = np.log(np.abs(xlike))
xmmin = xdata['model_min'][mask2]
xmmax = xdata['model_max'][mask2]
xdeb = xdata['delta_e_bulge'][mask2]
xdtb = xdata['delta_theta_bulge'][mask2]
xraas = xdata['ra_as'][mask2]
xdecas = xdata['dec_as'][mask2]
print 'Extracted all data fields'

m = 2.3

# For v2, this is a first-pass effort to color code "good" shear values.
good1 = (abs(m*tg1 - xe1) < 0.05) & (abs(m*tg2 - xe2) < 0.05)
good2 = (abs(m*tg1 - xe1) < 0.1) & (abs(m*tg2 - xe2) < 0.1) & ~good1
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

plt_scatter(tg1, xe1, 'True g1', 'im3shape e1', m=m)
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2', m=m)
plt_scatter(tg2, xe1, 'True g2', 'im3shape e1')
plt_scatter(tg1, xe2, 'True g1', 'im3shape e2')
plt_scatter(thlr, xr, 'True hlr', 'im3shape radius')
plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 10)
plt_scatter(tflux, xflux, 'True flux', 'im3shape bulge_flux + disc_flux')
plt_scatter(tflux, xflux, 'True flux', 'im3shape bulge_flux + disc_flux', abs(xflux) < 1.e2)
plt_scatter(tflux, xflux, 'True flux', 'im3shape bulge_flux + disc_flux', abs(xflux) < 1.e1)
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr')
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr', abs(xsnr) < 1.e5)
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr', abs(xsnr) < 1.e4)
plt_scatter(tmag, xmag, 'True mag', '-2.5 log10(flux)', xflux > 0)
plt_scatter(tmag, xmag2, 'True mag', '-2.5 log10(snr)', xsnr > 0)

#plt_scatter(xr, xrr, 'radius', 'radius_ratio')  # radius_ratio = 1
plt_scatter(xba, xda, 'bulge_A', 'disc_A')
plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e3) & (abs(xda) < 1.e3))
plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e2) & (abs(xda) < 1.e2))
#plt_scatter(xbi, xdi, 'bulge_index', 'disc_index')  # bulge_index = 4, disc_index = 1
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux')
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 1.e2) & (abs(xdf) < 1.e2))
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 1.e1) & (abs(xdf) < 1.e1))
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio')
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio', (abs(xbf) < 1.e2) & (abs(xdf) < 1.e2))
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio', (abs(xbf) < 1.e1) & (abs(xdf) < 1.e1))
plt_scatter(xsnr, xlike, 'snr', 'likelihood')
plt_scatter(xsnr, xlike, 'snr', 'likelihood', abs(xlike) < 1.e6)
plt_scatter(xsnr, xlike, 'snr', 'likelihood', abs(xlike) < 1.e5)
plt_scatter(xmag2, xlnlike, '-2.5 log10(snr)', 'ln(abs(likelihood))', xsnr > 0)
plt_scatter(xmmin, xmmax, 'model_min', 'model_max')
#plt_scatter(ideb, xdtb, 'delta_e_bulge', 'delta_theta_bulge')  # both = 0
plt_scatter(xraas, xdecas, 'ra_as', 'dec_as')
plt_scatter(xraas, xdecas, 'ra_as', 'dec_as', (abs(xraas) < 1) & (abs(xdecas) < 1))

good1 = (xba > 10) & (xda > 0.05)
good2 = (xba > 0) & (xda > 0) & ~good1
bad = (~good1) & (~good2)

plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e2) & (abs(xda) < 1.e2))
plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e3) & (abs(xda) < 1.e3))
plt_scatter(tg1, xe1, 'True g1', 'im3shape e1')
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2')
plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 10)
plt_scatter(tg1, xe1, 'True g1', 'im3shape e1',good1, m=m)
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2',good1, m=m)
plt_scatter(thlr, xr,'True hlr', 'im3shape radius', good1 & (xr<10))
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 1.e2) & (abs(xdf) < 1.e2))
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 1.e1) & (abs(xdf) < 1.e1))
plt_scatter(xsnr, xlike, 'snr', 'likelihood', abs(xlike) < 1.e5)
plt_scatter(xmag2, xlnlike, '-2.5 log10(snr)', 'ln(abs(likelihood))', xsnr > 0)

pp.close()
