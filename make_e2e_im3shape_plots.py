import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import sys

from matplotlib.backends.backend_pdf import PdfPages

dir = 'DES0436-5748'
orig_truth = pyfits.open(dir+'/end2end-truth.fits')[1].data
sex_cat = pyfits.open(dir+'/DES0436-5748_r_cat.fits')[1].data
match = pyfits.open(dir+'/match.fits')[1].data
xdata = pyfits.open('im3shape-end2end-v3.fits')[1].data

print 'orig_truth has %d columns, %d entries'%(len(orig_truth.columns), len(orig_truth))
print 'sex_cat has %d columns, %d entries'%(len(sex_cat.columns), len(sex_cat))
print 'match has %d columns, %d entries'%(len(match.columns), len(match))
print 'xdata has %d columns, %d entries'%(len(xdata.columns), len(xdata))

# The im3shape entries are not in order.  The identifier column - 1 gives the actual index to
# use for the truth entry.
truth = orig_truth[ match['index'][xdata['identifier']-1] ]
print 'truth has %d columns, %d entries'%(len(truth.columns), len(truth))

# Some of the items are flagged.  Remove these.
# Also select out just the galaxy objects.
mask = (match['ok'] == 1) & (truth['flags'] == 0) & (truth['is_star'] == 0) & (xdata['flag'] == 0)
print 'number of objects in original catalog = ',len(orig_truth)
print 'number of objects drawn = ',(orig_truth['flags'] == 0).sum()
print 'number of stars drawn = ',((orig_truth['flags'] == 0) &orig_truth['is_star']).sum()
print 'number detected by sextractor = ',len(sex_cat)
print 'number detected by sextractor with FLAGS==0: ',(sex_cat['FLAGS'] == 0).sum()
print 'number with good matches: ',match['ok'].sum()
print 'number of these that are stars = ',(match['ok'] & truth['is_star']).sum()
print 'number that were not actually drawn = ',(match['ok'] & (truth['flags'] != 0)).sum()
print 'number that im3shape marked as failure = ',(match['ok'] & (xdata['flag'] != 0)).sum()
print 'total passing all cuts = ',mask.sum()

# Extract values that we want to plot
tid = truth['id'][mask]
tg1 = truth['true_g1'][mask]
tg2 = truth['true_g2'][mask]
thlr = truth['true_hlr'][mask]
tra = truth['ra'][mask]
tdec = truth['dec'][mask]
tflux = truth['flux'][mask]
tmag = truth['mag'][mask]

xid = xdata['identifier'][mask]
xe1 = xdata['e1'][mask]
xe2 = xdata['e2'][mask]
xr = xdata['radius'][mask]
xra = xdata['ra'][mask]
xdec = xdata['dec'][mask]
xflux = xdata['mean_flux'][mask]
xsnr = xdata['snr'][mask]
xmag = -2.5*np.log10(xflux)
xmag[xflux <= 0] = 99
xmag2 = -2.5*np.log10(xsnr)
xmag2[xsnr <= 0] = 99
xrr = xdata['radius_ratio'][mask]
xba = xdata['bulge_A'][mask]
xda = xdata['disc_A'][mask]
xbi = xdata['bulge_index'][mask]
xdi = xdata['disc_index'][mask]
xbf = xdata['bulge_flux'][mask]
xdf = xdata['disc_flux'][mask]
xfr = xdata['flux_ratio'][mask]
xlike = xdata['likelihood'][mask]
xlnlike = np.log(np.abs(xlike))
xmmin = xdata['model_min'][mask]
xmmax = xdata['model_max'][mask]
xdeb = xdata['delta_e_bulge'][mask]
xdtb = xdata['delta_theta_bulge'][mask]
xraas = xdata['ra_as'][mask]
xdecas = xdata['dec_as'][mask]
print 'Extracted all data fields'

def simple_plots():
    """Make a few simple plots of truth vs meas
    """
    plt.clf()
    plt.axis([-0.3,0.3,-0.3,0.3])
    plt.grid()
    plt.xlabel('True g1')
    plt.ylabel('im3shape e1')
    plt.plot([-1.,-1.],[1.,1.],'c-')
    plt.scatter(truth['true_g1'],xdata['e1'],s=0.4,rasterized=True)
    plt.savefig('im3shape_e1.png')

    plt.clf()
    plt.axis([-0.3,0.3,-0.3,0.3])
    plt.grid()
    plt.xlabel('True g2')
    plt.ylabel('im3shape e2')
    plt.plot([-1.,-1.],[1.,1.],'c-')
    plt.scatter(truth['true_g2'],xdata['e2'],s=0.4,rasterized=True)
    plt.savefig('im3shape_e2.png')

simple_plots()
sys.exit()

pp = PdfPages('e2e_im3shape-3.pdf')
m = 1.0
tol1 = 0.01
tol2 = 0.05

# For v2, this is a first-pass effort to color code "good" shear values.
good1 = (abs(m*tg1 - xe1) < tol1) & (abs(m*tg2 - xe2) < tol1)
good2 = (abs(m*tg1 - xe1) < tol2) & (abs(m*tg2 - xe2) < tol2) & ~good1
bad = (~good1) & (~good2)


def plt_scatter(xval, yval, xlabel, ylabel, mask=None, m=None, title=None):
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
    plt.axis([min(xval), max(xval), min(yval), max(yval)])
    plt.grid()
    baddot = plt.scatter(xval[badm],yval[badm],s=0.4,rasterized=True,color='red')
    okdot = plt.scatter(xval[good2m],yval[good2m],s=0.4,rasterized=True,color='blue')
    gooddot = plt.scatter(xval[good1m],yval[good1m],s=0.4,rasterized=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if m is not None:
        if m > 1.:
            mline, = plt.plot([-1./m,1./m],[-1.,1.],'c-')
        else:
            mline, = plt.plot([-1.,1.],[-m,m],'c-')
        plt.legend([mline, gooddot, okdot, baddot], 
                   [
                        'm = %.2f'%m,
                        '|diff| < %.2f (N=%d)'%(tol1,good1.sum()),
                        '%.2f  < |diff| < %.2f (N=%d)'%(tol1,tol2,good2.sum()),
                        '|diff| > %.2f (N=%d)'%(tol2,bad.sum()),
                   ],
                   #loc = 2 # upper left
                   loc = 4 # lower right
                   )
    if title is not None:
        plt.title(title)
    pp.savefig()
    print 'Plotted %s vs %s'%(ylabel, xlabel)

plt_scatter(tg1, xe1, 'True g1', 'im3shape e1')
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2')
plt_scatter(tg1, xe1, 'True g1', 'im3shape e1', m=m)
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2', m=m)
plt_scatter(tg2, xe1, 'True g2', 'im3shape e1')
plt_scatter(tg1, xe2, 'True g1', 'im3shape e2')
plt_scatter(thlr, xr, 'True hlr', 'im3shape radius')
plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 10)
plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 2.5)
plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux')
plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux', (abs(xflux) < 1.e5) & (tflux < 1.e5))
plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux', (abs(xflux) < 2.e4) & (tflux < 2.e4))
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr')
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr',
            (abs(xsnr) < 1.e5) & (tflux < 1.e5))
plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr',
            (abs(xsnr) < 1.e4) & (tflux < 1.e4))
plt_scatter(tmag, xmag, 'True mag', '-2.5 log10(flux)', xflux > 0)
plt_scatter(tmag, xmag2, 'True mag', '-2.5 log10(snr)', xsnr > 0)

#plt_scatter(xr, xrr, 'radius', 'radius_ratio')  # radius_ratio = 1
plt_scatter(xba, xda, 'bulge_A', 'disc_A')
plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e3) & (abs(xda) < 1.e3))
plt_scatter(xba, xda, 'bulge_A', 'disc_A', (abs(xba) < 1.e2) & (abs(xda) < 1.e2))
#plt_scatter(xbi, xdi, 'bulge_index', 'disc_index')  # bulge_index = 4, disc_index = 1
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux')
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 1.e2) & (abs(xdf) < 1.e2))
plt_scatter(xbf, xdf, 'bulge_flux', 'disc_flux', (abs(xbf) < 2.) & (abs(xdf) < 4.))
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio')
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio', abs(xflux) < 1.e5)
plt_scatter(xflux, xfr, 'total flux', 'flux_ratio', abs(xflux) < 1.e4)
plt_scatter(xsnr, xlike, 'snr', 'likelihood')
plt_scatter(xsnr, xlike, 'snr', 'likelihood', abs(xlike) < 1.e6)
plt_scatter(xsnr, xlike, 'snr', 'likelihood', (abs(xlike) < 1.e5) & (xsnr < 1.e4))
plt_scatter(xmag2, xlnlike, '-2.5 log10(snr)', 'ln(abs(likelihood))', xsnr > 0)
plt_scatter(xmmin, xmmax, 'model_min', 'model_max')
plt_scatter(xmmin, xmmax, 'model_min', 'model_max', (abs(xmmin) < 0.001) & (abs(xmmax) < 0.15))
plt_scatter(xmmin, xmmax, 'model_min', 'model_max', (abs(xmmin) < 3.e-4) & (abs(xmmax) < 0.03))
plt_scatter(xmmin, xmmax, 'model_min', 'model_max', (abs(xmmin) < 3.e-6) & (abs(xmmax) < 0.03))
#plt_scatter(ideb, xdtb, 'delta_e_bulge', 'delta_theta_bulge')  # both = 0
plt_scatter(xraas, xdecas, 'ra_as', 'dec_as')
plt_scatter(xraas, xdecas, 'ra_as', 'dec_as', (abs(xraas) < 1) & (abs(xdecas) < 1))

mmin_limit = 1.e-6
plt_scatter(tg1, xe1, 'True g1', 'im3shape e1',
            (abs(xmmin) < mmin_limit), title='|model_min| < %f'%mmin_limit)
plt_scatter(tg2, xe2, 'True g2', 'im3shape e2',
            (abs(xmmin) < mmin_limit), title='|model_min| < %f'%mmin_limit)
plt_scatter(xmmin, xmmax, 'model_min', 'model_max',
            (abs(xmmin) < mmin_limit), title='|model_min| < mmin_limit, (N=%d)'%((abs(xmmin)<mmin_limit).sum()))

pp.close()
