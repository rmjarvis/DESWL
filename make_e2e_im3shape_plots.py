import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import sys
from matplotlib.backends.backend_pdf import PdfPages

dir = 'DES0436-5748'
truth_file = 'end2end-truth.fits'
sex_file = 'DES0436-5748_r_cat.fits'
match_file = 'match.fits'
im3shape_file = 'im3shape-end2end-v3.fits'

class Truth(object):
    def __init__(self, dir, file_name, sex_name, match_name):
        import astropy.io.fits as pyfits
        import os
        self.orig = pyfits.open(os.path.join(dir,file_name))[1].data
        self.sex = pyfits.open(os.path.join(dir,sex_name))[1].data
        self.match = pyfits.open(os.path.join(dir,match_name))[1].data

        self.name = 'truth'

        print 'orig_truth has %d columns, %d entries'%(len(self.orig.columns), len(self.orig))
        print 'sex_cat has %d columns, %d entries'%(len(self.sex.columns), len(self.sex))
        print 'match has %d columns, %d entries'%(len(self.match.columns), len(self.match))

        self.match_index = self.match['index']

        self.ok = ( (self.match['ok'] == 1) & 
                    (self.orig['flags'][self.match_index] == 0) &
                    (self.orig['is_star'][self.match_index] == 0) )
        print 'number of objects in original catalog = ',len(self.orig)
        print 'number of objects drawn = ',(self.orig['flags'] == 0).sum()
        print 'number of stars drawn = ',((self.orig['flags'] == 0) & self.orig['is_star']).sum()
        print 'number detected by sextractor = ',len(self.sex)
        print 'number detected by sextractor with FLAGS==0: ',(self.sex['FLAGS'] == 0).sum()
        print 'number with good matches: ',self.match['ok'].sum()
        print 'number of these that are stars = ',(self.match['ok'] & self['is_star']).sum()
        print 'number that were not drawn = ',(self.match['ok'] & (self['flags'] != 0)).sum()
        print 'truth has %d entries listed as ok'%(len(self.ok))

        # pull out the basic information into attributes, since we might want to access
        # these values multiple times.
        self.id = self.orig['id'][self.match_index]
        self.g1 = self.orig['true_g1'][self.match_index]
        self.g2 = self.orig['true_g2'][self.match_index]
        self.r = self.orig['true_hlr'][self.match_index]
        self.flux = self.orig['flux'][self.match_index]
        self.mag = self.orig['mag'][self.match_index]
        self.ra = self.orig['ra'][self.match_index]
        self.dec = self.orig['dec'][self.match_index]

    def __getitem__(self, key): 
        return self.orig[key][self.match_index]

class Im3Shape(object):
    def __init__(self, file_name):
        import astropy.io.fits as pyfits
        self.data = pyfits.open(im3shape_file)[1].data
        print 'im3shape has %d columns, %d entries'%(len(self.data.columns), len(self.data))

        self.name = 'im3shape'

        # The im3shape entries are not in order.  The identifier column - 1 gives the actual index
        # to use for the truth entry.
        self.index = self.data['identifier']-1

        self.ok = (self.data['flag'] == 0)

        # pull out the basic information that should be in all catalogs
        self.id = self.data['identifier']
        self.g1 = self.data['e1']
        self.g2 = self.data['e2']
        self.r = self.data['radius']
        self.flux = self.data['mean_flux']
        self.snr = self.data['snr']
        self.ra = self.data['ra']
        self.dec = self.data['dec']

    def __getitem__(self, key): 
        return self.data[key]

def simple_plots(truth, meas):
    """Make a few simple plots of truth vs meas
    """
    mask = truth.ok[meas.index] & meas.ok
    tg1 = truth.g1[meas.index][mask]
    tg2 = truth.g2[meas.index][mask]
    xg1 = meas.g1[mask]
    xg2 = meas.g2[mask]

    plt.clf()
    plt.axis([-0.3,0.3,-0.3,0.3])
    plt.grid()
    plt.xlabel('True g1')
    plt.ylabel('{0} e1'.format(meas.name))
    plt.plot([-1.,-1.],[1.,1.],'c-')
    plt.scatter(tg1,xg1,s=0.4,rasterized=True)
    plt.savefig('{0}_e1.png'.format(meas.name))

    plt.clf()
    plt.axis([-0.3,0.3,-0.3,0.3])
    plt.grid()
    plt.xlabel('True g2')
    plt.ylabel('{0} e2'.format(meas.name))
    plt.plot([-1.,-1.],[1.,1.],'c-')
    plt.scatter(tg2,xg2,s=0.4,rasterized=True)
    plt.savefig('{0}_e2.png'.format(meas.name))


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


def do_basic_plots(truth, meas):

    mask = truth.ok[meas.index] & meas.ok
    tg1 = truth.g1[meas.index][mask]
    tg2 = truth.g2[meas.index][mask]
    xg1 = meas.g1[mask]
    xg2 = meas.g2[mask]

    global tol1, tol2, good1, good2, bad

    tol1 = 0.01
    tol2 = 0.05
    good1 = (abs(tg1 - xg1) < tol1) & (abs(tg2 - xg2) < tol1)
    good2 = (abs(tg1 - xg1) < tol2) & (abs(tg2 - xg2) < tol2) & ~good1
    bad = (~good1) & (~good2)

    print 'number that %s marked as failure = %d'%(meas.name,(truth.ok[meas.index]&~meas.ok).sum())
    print 'total passing all cuts = ',mask.sum()
    print 'number within %f of correct shape = %d'%(tol1, good1.sum())

    # extract values that we want to plot
    tid = truth.id[meas.index][mask]
    tg1 = truth.g1[meas.index][mask]
    tg2 = truth.g2[meas.index][mask]
    thlr = truth.r[meas.index][mask]
    tflux = truth.flux[meas.index][mask]
    tmag = truth.mag[meas.index][mask]

    xid = meas.id[mask]
    xg1 = meas.g1[mask]
    xg2 = meas.g2[mask]
    xr = meas.r[mask]
    xflux = meas.flux[mask]
    xsnr = meas.snr[mask]
    xmag = -2.5*np.log10(xflux)
    xmag[xflux <= 0] = 99
    xmag2 = -2.5*np.log10(xsnr)
    xmag2[xsnr <= 0] = 99


    plt_scatter(tg1, xg1, 'True g1', 'im3shape e1')
    plt_scatter(tg2, xg2, 'True g2', 'im3shape e2')
    plt_scatter(tg1, xg1, 'True g1', 'im3shape e1', m=1)
    plt_scatter(tg2, xg2, 'True g2', 'im3shape e2', m=1)
    plt_scatter(tg2, xg1, 'True g2', 'im3shape e1')
    plt_scatter(tg1, xg2, 'True g1', 'im3shape e2')
    plt_scatter(thlr, xr, 'True hlr', 'im3shape radius')
    plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 10)
    plt_scatter(thlr, xr,'True hlr', 'im3shape radius', xr < 2.5)
    plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux')
    plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux',
                (abs(xflux) < 1.e5) & (tflux < 1.e5))
    plt_scatter(tflux, xflux, 'True flux', 'im3shape mean_flux',
                (abs(xflux) < 2.e4) & (tflux < 2.e4))
    plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr')
    plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr',
                (abs(xsnr) < 1.e5) & (tflux < 1.e5))
    plt_scatter(tflux, xsnr, 'True flux', 'im3shape snr',
                (abs(xsnr) < 1.e4) & (tflux < 1.e4))
    plt_scatter(tmag, xmag, 'True mag', '-2.5 log10(flux)', xflux > 0)
    plt_scatter(tmag, xmag2, 'True mag', '-2.5 log10(snr)', xsnr > 0)

    if hasattr(meas,'ra'):
        tra = truth.ra[meas.index][mask]
        tdec = truth.dec[meas.index][mask]
        xra = meas.ra[mask]
        xdec = meas.dec[mask]

        plt_scatter(tra, xra, 'True RA', 'im3shape RA')
        plt_scatter(tdec, xdec, 'True Dec', 'im3shape Dec')


def do_extra_im3shape_plots(truth, meas):

    mask = truth.ok[meas.index] & meas.ok
    tg1 = truth.g1[meas.index][mask]
    tg2 = truth.g2[meas.index][mask]
    xg1 = meas.g1[mask]
    xg2 = meas.g2[mask]

    xid = meas.id[mask]
    xr = meas.r[mask]
    xflux = meas.flux[mask]
    xsnr = meas.snr[mask]
    xmag = -2.5*np.log10(xflux)
    xmag[xflux <= 0] = 99
    xmag2 = -2.5*np.log10(xsnr)
    xmag2[xsnr <= 0] = 99

    xrr = im3shape['radius_ratio'][mask]
    xba = im3shape['bulge_A'][mask]
    xda = im3shape['disc_A'][mask]
    xbi = im3shape['bulge_index'][mask]
    xdi = im3shape['disc_index'][mask]
    xbf = im3shape['bulge_flux'][mask]
    xdf = im3shape['disc_flux'][mask]
    xfr = im3shape['flux_ratio'][mask]
    xlike = im3shape['likelihood'][mask]
    xlnlike = np.log(np.abs(xlike))
    xmmin = im3shape['model_min'][mask]
    xmmax = im3shape['model_max'][mask]
    xdeb = im3shape['delta_e_bulge'][mask]
    xdtb = im3shape['delta_theta_bulge'][mask]
    xraas = im3shape['ra_as'][mask]
    xdecas = im3shape['dec_as'][mask]
    print 'Extracted all data fields'

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
    plt_scatter(xraas, xdecas, 'ra_as', 'dec_as',
                (abs(xraas-0.26) < 0.1) & (abs(xdecas+0.26) < 0.1))

    mmin_limit = 1.e-6
    plt_scatter(tg1, xg1, 'True g1', 'im3shape e1',
                (abs(xmmin) < mmin_limit), title='|model_min| < %f'%mmin_limit)
    plt_scatter(tg2, xg2, 'True g2', 'im3shape e2',
                (abs(xmmin) < mmin_limit), title='|model_min| < %f'%mmin_limit)
    plt_scatter(xmmin, xmmax, 'model_min', 'model_max',
                (abs(xmmin) < mmin_limit),
                title='|model_min| < %f, (N=%d)'%(mmin_limit,(abs(xmmin)<mmin_limit).sum()))



truth = Truth(dir, truth_file, sex_file, match_file)
im3shape = Im3Shape(im3shape_file)

simple_plots(truth, im3shape)

pp = PdfPages('e2e_im3shape-3.pdf')

do_basic_plots(truth, im3shape)
do_extra_im3shape_plots(truth, im3shape)

pp.close()


