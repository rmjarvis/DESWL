import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('e2e_nfit-2.pdf')

truth = pyfits.open('end2end-truth.fits')[1]
nfit = pyfits.open('test-end-to-end-nfit-01.fits')[1]

tdata = truth.data
xdata = nfit.data
print 'tdata has %d columns, %d entries'%(len(tdata.columns), len(tdata))
print 'xdata has %d columns, %d entries'%(len(xdata.columns), len(xdata))

# Some of the items are flagged.  Remove these.
# Also select out just the galaxy objects.
mask = (tdata['flags'] == 0) & (tdata['is_star'] == 0) & (xdata['exp_flags'] == 0)

# Extract values that we want to plot
tid = tdata['id'][mask]
tg1 = tdata['true_g1'][mask]
tg2 = tdata['true_g2'][mask]
thlr = tdata['true_hlr'][mask]
tra = tdata['ra'][mask]
tdec = tdata['dec'][mask]
tflux = tdata['flux'][mask]
tmag = tdata['mag'][mask]

xid = xdata['number'][mask]
xg1 = xdata['exp_g'][mask][:,0]
xg2 = xdata['exp_g'][mask][:,1]
xcenx = xdata['exp_pars'][mask][:,0]
xceny = xdata['exp_pars'][mask][:,1]
xcen = np.sqrt(xcenx**2+xceny**2)
xe1 = xdata['exp_pars'][mask][:,2]
xe2 = xdata['exp_pars'][mask][:,3]
xr2 = xdata['exp_pars'][mask][:,4]
xr = np.sqrt(xr2)
xr[xr2<0] = -1
xflux = xdata['exp_pars'][mask][:,5]
xsnr = xdata['exp_s2n_w'][mask]
xmag = -2.5*np.log10(xflux)
xmag[xflux <= 0] = 99
xmag2 = -2.5*np.log10(xsnr)
xmag2[xsnr <= 0] = 99
xchi2 = xdata['exp_chi2per'][mask]
xdof = xdata['exp_dof'][mask]
xaic = xdata['exp_aic'][mask]
xbic = xdata['exp_bic'][mask]
xarate = xdata['exp_arate'][mask]
xtau = xdata['exp_tau'][mask]
xp = xdata['exp_P'][mask]
xq0 = xdata['exp_Q'][mask][:,0]
xq1 = xdata['exp_Q'][mask][:,1]
xq = np.sqrt(xq0**2+xq1**2)
xr00 = xdata['exp_R'][mask][:,0,0]
xr01 = xdata['exp_R'][mask][:,0,1]
xr10 = xdata['exp_R'][mask][:,1,0]
xr11 = xdata['exp_R'][mask][:,1,1]
xtrr = xr00+xr11
xdetr = xr00*xr11-xr01*xr10
print 'Extracted all data fields'

m1 = -1
m2 = 1

# For v2, this is a first-pass effort to color code "good" shear values.
good1 = (abs(m1*tg1 - xg1) < 0.05) & (abs(m2*tg2 - xg2) < 0.05)
good2 = (abs(m1*tg1 - xg1) < 0.1) & (abs(m2*tg2 - xg2) < 0.1) & ~good1
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
        if abs(m) > 1.:
            mline, = plt.plot([-1./m,1./m],[-1.,1.],'c-')
        else:
            mline, = plt.plot([-1.,1.],[-m,m],'c-')
        plt.legend([mline], ['m = %f'%m])
    pp.savefig()
    print 'Plotted %s vs %s'%(ylabel, xlabel)

plt_scatter(tg1, xg1, 'True g1', 'nfit g1', m=m1)
plt_scatter(tg2, xg2, 'True g2', 'nfit g2', m=m2)
plt_scatter(tg2, xg1, 'True g2', 'nfit g1')
plt_scatter(tg1, xg2, 'True g1', 'nfit g2')
plt_scatter(tg1, xe1, 'True g1', 'nfit raw e1', m=m1)
plt_scatter(tg2, xe2, 'True g2', 'nfit raw e2', m=m2)
plt_scatter(tg2, xe1, 'True g2', 'nfit raw e1')
plt_scatter(tg1, xe2, 'True g1', 'nfit raw e2')
plt_scatter(thlr, xr2, 'True hlr', 'nfit Irr')
plt_scatter(thlr, xr2, 'True hlr', 'nfit Irr', xr2<10)
plt_scatter(thlr, xr,'True hlr', 'sqrt(Irr)', xr2>0)
plt_scatter(thlr, xr,'True hlr', 'sqrt(Irr)', (xr2>0) & (xr2<10))
plt_scatter(thlr, xr,'True hlr', 'sqrt(Irr)', (xr2>0) & (xr2<1))
plt_scatter(tflux, xflux, 'True flux', 'nfit flux')
plt_scatter(tflux, xflux, 'True flux', 'nfit flux', xflux<1.e3)
plt_scatter(tflux, xsnr, 'True flux', 'nfit snr')
plt_scatter(tflux, xsnr, 'True flux', 'nfit snr', (xsnr<1.e3) & (tflux < 1.e4))
plt_scatter(tmag, xmag, 'True mag', '-2.5 log10(flux)', xflux > 0)
plt_scatter(tmag, xmag2, 'True mag', '-2.5 log10(snr)', xsnr > 0)

plt_scatter(xcenx, xceny, 'centroid_x', 'centroid_y')
plt_scatter(xcenx, xceny, 'centroid_x', 'centroid_y', xcen<1)
plt_scatter(xcenx, xceny, 'centroid_x', 'centroid_y', xcen<0.2)
plt_scatter(xcenx, xceny, 'centroid_x', 'centroid_y', xcen<0.05)
plt_scatter(xchi2, xsnr, 'chi2', 'snr')
plt_scatter(xchi2, xsnr, 'chi2', 'snr', xchi2<1.e4)
plt_scatter(xchi2, xsnr, 'chi2', 'snr', xchi2<1.e3)
plt_scatter(xchi2, xdof, 'chi2', 'dof')
plt_scatter(xchi2, xdof, 'chi2', 'dof', xchi2<1.e4)
plt_scatter(xchi2, xdof, 'chi2', 'dof', xchi2<1.e3)
plt_scatter(xaic, xbic, 'aic', 'bic')
plt_scatter(xaic, xbic, 'aic', 'bic', xaic<1.e6)
plt_scatter(xarate, xaic, 'arate', 'aic', xaic<1.e6)
plt_scatter(xarate, xtau, 'arate', 'tau')
plt_scatter(xarate, xr2, 'arate', 'Irr')
plt_scatter(xarate, xr2, 'arate', 'Irr',xr2<10)
plt_scatter(xp, xtau, 'P', 'tau')
plt_scatter(xp, xarate, 'P', 'arate')
plt_scatter(xp, xr2, 'P', 'Irr')
plt_scatter(xp, xr2, 'P', 'Irr',xr2<10)
plt_scatter(xq0, xq1, 'Q0', 'Q1')
plt_scatter(xq0, xq1, 'Q0', 'Q1', xq<1.e3)
plt_scatter(xq0, xq1, 'Q0', 'Q1', xq<300)
plt_scatter(xp, xq, 'P', '|Q|')
plt_scatter(xp, xq, 'P', '|Q|', xq<1.e3)
plt_scatter(xr00, xr11, 'R00', 'R11')
plt_scatter(xr00, xr11, 'R00', 'R11',xtrr<1.e6)
plt_scatter(xr00, xr11, 'R00', 'R11',xtrr<1.e5)
plt_scatter(xr00, xr11, 'R00', 'R11',(abs(xr00)<1.e4) & (abs(xr11)<1.e4))
plt_scatter(xr00, xr11, 'R00', 'R11',xr00**2+xr11**2<3.e3**2)
plt_scatter(xr00, xr01, 'R00', 'R01')
plt_scatter(xr00, xr01, 'R00', 'R01',xtrr<1.e6)
plt_scatter(xr00, xr01, 'R00', 'R01',(abs(xr00)<1.e5) & (abs(xr01)<1.e5))
plt_scatter(xr00, xr01, 'R00', 'R01',(abs(xr00)<1.e4) & (abs(xr01)<1.e4))
plt_scatter(xp, xtrr, 'P', 'Tr(R)')
plt_scatter(xp, xtrr, 'P', 'Tr(R)',xtrr<1.e6)
plt_scatter(xp, xtrr, 'P', 'Tr(R)',xtrr<1.e5)
plt_scatter(xp, xtrr, 'P', 'Tr(R)',(xtrr<1.e5) & (xp<500))
plt_scatter(xp, xdetr, 'P', 'Det(R)')
plt_scatter(xp, xdetr, 'P', 'Det(R)',xtrr<1.e5)
plt_scatter(xp, xdetr, 'P', 'Det(R)',(xtrr<1.e5) & (xp<1.e4) & (abs(xdetr)<1.e8))
plt_scatter(xp, xdetr, 'P', 'Det(R)',(xtrr<1.e5) & (xp<500) & (abs(xdetr)<1.e7))
plt_scatter(xq, xtrr, '|Q|', 'Tr(R)')
plt_scatter(xq, xtrr, '|Q|', 'Tr(R)',xtrr<1.e6)
plt_scatter(xq, xtrr, '|Q|', 'Tr(R)',xtrr<1.e5)
plt_scatter(xq, xtrr, '|Q|', 'Tr(R)',(xtrr<1.e5) & (xq<1000))
plt_scatter(xq, xdetr, '|Q|', 'Det(R)')
plt_scatter(xq, xdetr, '|Q|', 'Det(R)',xtrr<1.e5)
plt_scatter(xq, xdetr, '|Q|', 'Det(R)',(xtrr<1.e5) & (xq<1.e4) & (abs(xdetr)<1.e8))
plt_scatter(xq, xdetr, '|Q|', 'Det(R)',(xtrr<1.e5) & (xq<1000) & (abs(xdetr)<1.e7))

pp.close()
