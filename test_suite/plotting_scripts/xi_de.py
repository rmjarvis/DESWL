import pyfits
import numpy
import treecorr

def plot_dxi(gg, filename, sqrtn=1, label=None, ximinus=False):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mp
    plt.style.use('supermongo')
    plt.minorticks_on()

    meanr = numpy.exp(gg.meanlogr)
    xip = gg.xip
    xim = gg.xim
    sig = numpy.sqrt(gg.varxi)

    # Use Alexandre's xi file as a requirement.
    theta = [0.5, 1.49530470404, 2.91321328166, 5.67019971545, 11.0321639144,
                21.4548924111, 41.6931936543, 80.8508152859, 156.285886576,
                297.92139021, 300.]

    dxip = [1.4e-5, 3.82239654447e-06, 2.36185315415e-06, 1.26849547074e-06, 6.3282672138e-07,
            3.25623661098e-07, 1.747852053e-07, 8.75326181278e-08, 3.60247306537e-08,
            1.13521735321e-08, 1.125e-8]
    dxim = [5.1e-07, 5.1915039397e-07, 5.85603620143e-07, 5.69416647483e-07, 4.49138135875e-07,
            2.84585155218e-07, 1.5435276787e-07, 8.127342852e-08, 4.58443061967e-08,
            2.65823041365e-08, 2.6e-08]

    plt, ax = plt.subplots(1,1, sharex=True)

    pos_xip = xip.copy()
    pos_xip[xip<0] = 1.e-9
    neg_xip = -xip
    neg_xip[xip>0] = 1.e-9
    pos_xim = xim.copy()
    pos_xim[xim<0] = 1.e-9
    neg_xim = -xim
    neg_xim[xim>0] = 1.e-9

    #ax = axes[0]
    ax.fill( [theta[0]] + theta + [theta[-1]], [0.] + dxip + [0.], color = '#FFFF82')
    ax.plot(meanr, pos_xip, color='blue')
    ax.plot(meanr, neg_xip, color='blue', ls=':')
    ax.errorbar(meanr[xip>0], xip[xip>0], yerr=sig[xip>0]/sqrtn, color='blue', ls='')
    ax.errorbar(meanr[xip<0], -xip[xip<0], yerr=sig[xip<0]/sqrtn, color='blue', ls='')
    xip_line, = ax.plot(-meanr, xip, color='blue')
    xip_req = mp.Patch(color='#FFFF82')
    # Note: latex stuff doesn't work if you put it in the label='' form in each item.
    # But putting the strings explicitly in the legend works fine.
    if False:
        ax.legend([xip_line, xip_req],
                    [r'$\xi_+(\theta)$', r'$\xi_+$ requirement'],
                    loc='upper right', fontsize=24)
    ax.set_ylim( [5.e-8, 2.e-5] )
    if label is None:
        if ximinus:
            ax.set_ylabel(r'$\xi_{\Delta e}(\theta)$', fontsize=24)
        else:
            ax.set_ylabel(r'$\xi_{+\Delta e}(\theta)$', fontsize=24)
    else:
        ax.set_ylabel(label)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log')

    if ximinus:
        #ax = axes[1]
        #ax.fill( [theta[0]] + theta + [theta[-1]], [0.] + dxim + [0.], color = 'Lime')
        ax.plot(meanr, pos_xim, color='red')
        ax.plot(meanr, neg_xim, color='red', ls=':')
        ax.errorbar(meanr[xim>0], xim[xim>0], yerr=sig[xim>0]/sqrtn, color='red', ls='')
        ax.errorbar(meanr[xim<0], -xim[xim<0], yerr=sig[xim<0]/sqrtn, color='red', ls='')
        xim_line, = ax.plot(-meanr, xim, color='red')
        #xim_req = mp.Patch(color='#FFFF82')
        #xim_req = mp.Patch(color='Lime')
        #ax.legend([xim_line, xim_req],
                    #[r'$\xi_-(\theta)$', r'$\xi_-$ requirement'],
                    #loc='upper right')
        #ax.set_ylim( [1.e-8, 2.e-6] )
        #ax.set_ylim( [1.e-8, 2.e-5] )
        #ax.set_ylabel(r'$\xi_{+,-,\Delta e}(\theta)$')
        #ax.set_yscale('log', nonposy='clip')
        #ax.set_xscale('log')
        ax.legend([xip_line, xim_line],
                    [r'$\xi_+(\theta)$', r'$\xi_-(\theta)$'],
                    loc='upper right', fontsize=24)

    ax.set_xlim( [1.,100.] )
    ax.set_xlabel(r'$\theta$ (arcmin)', fontsize=24)

    plt.tight_layout()
    #plt.set_size_inches(8., 12.)
    plt.savefig(filename)


def plot_xi(gg, filename):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mp
    plt.style.use('supermongo')
    plt.minorticks_on()

    meanr = numpy.exp(gg.meanlogr)
    xip = gg.xip
    xim = gg.xim
    sig = numpy.sqrt(gg.varxi)

    plt, ax = plt.subplots(1,1, sharex=True)

    pos_xip = xip.copy()
    pos_xip[xip<0] = 1.e-9
    neg_xip = -xip
    neg_xip[xip>0] = 1.e-9
    pos_xim = xim.copy()
    pos_xim[xim<0] = 1.e-9
    neg_xim = -xim
    neg_xim[xim>0] = 1.e-9

    ax.plot(meanr, pos_xip, color='blue')
    ax.plot(meanr, neg_xip, color='blue', ls=':')
    ax.errorbar(meanr[xip>0], xip[xip>0], yerr=sig[xip>0], color='blue', ls='')
    ax.errorbar(meanr[xip<0], -xip[xip<0], yerr=sig[xip<0], color='blue', ls='')
    xip_line, = ax.plot(-meanr, xip, color='blue')
    #ax.set_ylim( [5.e-8, 2.e-5] )
    ax.set_ylabel(r'$\xi(\theta)$')
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log')

    if True:
        ax.plot(meanr, pos_xim, color='red')
        ax.plot(meanr, neg_xim, color='red', ls=':')
        ax.errorbar(meanr[xim>0], xim[xim>0], yerr=sig[xim>0], color='red', ls='')
        ax.errorbar(meanr[xim<0], -xim[xim<0], yerr=sig[xim<0], color='red', ls='')
        xim_line, = ax.plot(-meanr, xim, color='red')
        ax.legend([xip_line, xim_line],
                    [r'$\xi_+(\theta)$', r'$\xi_-(\theta)$'],
                    loc='upper right', fontsize=24)

    ax.set_xlim( [1.,500.] )
    ax.set_xlabel(r'$\theta$ (arcmin)')

    plt.tight_layout()
    #plt.set_size_inches(8., 12.)
    plt.savefig(filename)


def calculateEB(gg):
    # Make s a matrix, so we can eventually do the integral by doing a matrix product.
    r = numpy.exp(gg.logr)
    meanr = numpy.exp(gg.meanlogr) # Use the actual mean r for each bin
    s = numpy.outer(1./r, meanr)  
    ssq = s*s

    # xip_E(R) = 1/2 (xip(R) + xim(R) + int(dlogr xim(r) (4 - 12 R^2/r^2) H(r-R)
    # Let Gm be the factor we multiply by xim.  Note: s = r/R
    Gm = numpy.zeros_like(s)
    Gm[s>1.] = 4. - 12./ssq[s>1.]

    # xim_E(R) = 1/2 (xip(R) + xim(R) + int(dlogr xip(r) (4 - 12 r^2/R^2) r^2/R^2 H(R-r)
    # Let Gp be the factor we multiply by xip.
    Gp = numpy.zeros_like(s)
    Gp[s<1.] = (4. - 12.*ssq[s<1.]) * ssq[s<1.]

    # Now do the integral by taking the matrix products.
    # Note that dlogr = bin_size
    Gpxip = Gp.dot(gg.xip) * gg.bin_size
    Gmxim = Gm.dot(gg.xim) * gg.bin_size

    # The variance of each of these is the original varxi plus an extra term
    # dGpxip = int_r=0..2R [1/4 dlogr^2 (T+(s)^2 + T-(s)^2) Var(xi)]
    varGpxip = (Gp**2).dot(gg.varxi) * 0.25 * gg.bin_size**2
    varGmxim = (Gm**2).dot(gg.varxi) * 0.25 * gg.bin_size**2

    varxip = gg.varxi + varGpxip
    varxim = gg.varxi + varGmxim

    xipE = 0.5 * (gg.xip + gg.xim + Gmxim)
    ximE = 0.5 * (gg.xip + gg.xim + Gpxip)
    xipB = gg.xip - xipE
    ximB = gg.xim - ximE

    # In both cases, there is a part of the integral we cannot do.  We make the ansatz that
    # the B mode at the largest scales is consistent with zero.
    nmean = int(numpy.ceil(numpy.log(10) / gg.bin_size))  # Average over the last factor of 2 in scales.

    # For xip, using the Gmxim integral, the unknown part is of the form a + b R^2
    y = xipB[-nmean:]
    print 'y = ',y
    x = r[-nmean:]
    print 'x = ',x
    sw = 1./numpy.sqrt(varxip[-nmean:])
    A = numpy.vstack([ numpy.ones(len(x)) * sw, x**2 * sw ]).T
    print 'A = ',A
    a,b = numpy.linalg.lstsq(A,y*sw)[0]
    print 'a,b = ',a,b

    xipE += a + b*r**2
    xipB -= a + b*r**2
    print 'xipE -> ',xipE
    print 'xipB -> ',xipB

    # For xim, using the Gpxip integral, the unknown part is of the form a/R^2 + b/R^4
    y = numpy.concatenate([ximB[:nmean/2], ximB[-nmean:]])
    x = numpy.concatenate([r[:nmean/2], r[-nmean:]])
    sw = 1./numpy.sqrt(numpy.concatenate([varxim[:nmean/2], varxim[-nmean:]]))
    print 'y = ',y
    A = numpy.vstack([ x**-2 * sw, x**-4 * sw]).T
    print 'A = ',A
    a,b = numpy.linalg.lstsq(A,y*sw)[0]
    print 'a,b = ',a,b

    ximE += a*r**-2 + b*r**-4
    ximB -= a*r**-2 + b*r**-4

    return xipE, xipB, ximE, ximB, varxip, varxim

def smooth(xi, n):
    """Return a smoother version of xi where the values are the average of +-n points."""
    M = numpy.zeros( (len(xi), len(xi)) )
    x,y = numpy.indices( (len(xi), len(xi)) )
    M[ numpy.abs(x-y) < n ] = 1
    return M.dot(xi) / M.dot(numpy.ones_like(xi))


def plot_EB(gg, filename):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mp
    plt.style.use('supermongo')
    plt.minorticks_on()

    meanr = numpy.exp(gg.meanlogr)
    (xip_E, xip_B, xim_E, xim_B, varxip, varxim) = calculateEB(gg)

    # These are still the original density.  Sparsify by factor of 5.

    #xip_E = gg.xip
    #xip_B = gg.xim
    #varxip = gg.varxi

    n = int(numpy.ceil(numpy.log(5) / gg.bin_size))  # Average over the last factor of 2 in scales.
    xip_E = smooth(xip_E,n)
    xim_E = smooth(xim_E,n)
    xip_B = smooth(xip_B,n)
    xim_B = smooth(xim_B,n)
    varxip = smooth(varxip,n)/(2*n)
    varxim = smooth(varxim,n)/(2*n)

    plt, ax = plt.subplots(1,1, sharex=True)

    pos_xip_E = xip_E.copy()
    pos_xip_E[xip_E<0] = 1.e-14
    neg_xip_E = -xip_E
    neg_xip_E[xip_E>0] = 1.e-14
    pos_xip_B = xip_B.copy()
    pos_xip_B[xip_B<0] = 1.e-14
    neg_xip_B = -xip_B
    neg_xip_B[xip_B>0] = 1.e-14

    pos_xim_E = xim_E.copy()
    pos_xim_E[xim_E<0] = 1.e-14
    neg_xim_E = -xim_E
    neg_xim_E[xim_E>0] = 1.e-14
    pos_xim_B = xim_B.copy()
    pos_xim_B[xim_B<0] = 1.e-14
    neg_xim_B = -xim_B
    neg_xim_B[xim_B>0] = 1.e-14

    sigp = numpy.sqrt(varxip)
    sigm = numpy.sqrt(varxip)

    n2 = int(numpy.ceil(numpy.log(2) / gg.bin_size))
    e = slice(0,len(xip_E),n2)

    ax.plot(meanr, pos_xip_E, color='blue')
    ax.plot(meanr, neg_xip_E, color='blue', ls=':')
    ax.errorbar(meanr[e][xip_E[e]>0], xip_E[e][xip_E[e]>0], yerr=sigp[e][xip_E[e]>0],
                color='blue', ls='')
    ax.errorbar(meanr[e][xip_E[e]<0], -xip_E[e][xip_E[e]<0], yerr=sigp[e][xip_E[e]<0],
                color='blue', ls='')
    xip_E_line, = ax.plot(-meanr, xip_E, color='blue')

    ax.plot(meanr, pos_xim_E, color='cyan')
    ax.plot(meanr, neg_xim_E, color='cyan', ls=':')
    ax.errorbar(meanr[e][xim_E[e]>0], xim_E[e][xim_E[e]>0], yerr=sigm[e][xim_E[e]>0],
                color='cyan', ls='')
    ax.errorbar(meanr[e][xim_E[e]<0], -xim_E[e][xim_E[e]<0], yerr=sigm[e][xim_E[e]<0],
                color='cyan', ls='')
    xim_E_line, = ax.plot(-meanr, xim_E, color='cyan')

    ax.plot(meanr, pos_xip_B, color='red')
    ax.plot(meanr, neg_xip_B, color='red', ls=':')
    ax.errorbar(meanr[e][xip_B[e]>0], xip_B[e][xip_B[e]>0], yerr=sigp[e][xip_B[e]>0],
                color='red', ls='')
    ax.errorbar(meanr[e][xip_B[e]<0], -xip_B[e][xip_B[e]<0], yerr=sigp[e][xip_B[e]<0],
                color='red', ls='')
    xip_B_line, = ax.plot(-meanr, xip_B, color='red')

    ax.plot(meanr, pos_xim_B, color='magenta')
    ax.plot(meanr, neg_xim_B, color='magenta', ls=':')
    ax.errorbar(meanr[e][xim_B[e]>0], xim_B[e][xim_B[e]>0], yerr=sigm[e][xim_B[e]>0],
                color='magenta', ls='')
    ax.errorbar(meanr[e][xim_B[e]<0], -xim_B[e][xim_B[e]<0], yerr=sigm[e][xim_B[e]<0],
                color='magenta', ls='')
    xim_B_line, = ax.plot(-meanr, xim_B, color='magenta')

    ax.legend([xip_E_line, xim_E_line, xip_B_line, xim_B_line],
                [r'$\xi_{+\rm E}(\theta)$', r'$\xi_{-\rm E}(\theta)$',
                 r'$\xi_{+\rm B}(\theta)$', r'$\xi_{-\rm B}(\theta)$'],
                loc='upper right', fontsize=24)

    ax.set_ylim( [1.e-9, 6.e-4] )
    ax.set_ylabel(r'$\xi_{\rm E,B}(\theta)$')
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log')
    ax.set_xlim( [1.,500.] )
    ax.set_xlabel(r'$\theta$ (arcmin)')

    plt.tight_layout()
    #plt.set_size_inches(8., 12.)
    plt.savefig(filename)



def plot_map(gg, filename):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mp
    plt.style.use('supermongo')
    plt.minorticks_on()

    meanr = numpy.exp(gg.meanlogr)
    (mapsq, mapsq_im, mxsq, mxsq_im, varmapsq) = gg.calculateMapSq()

    # These are still the original density.  Sparsify by factor of 5.

    meanr = meanr[::5]
    mapsq = mapsq[::5]
    mxsq = mxsq[::5]
    varmapsq = varmapsq[::5]

    plt, ax = plt.subplots(1,1, sharex=True)

    pos_mapsq = mapsq.copy()
    pos_mapsq[mapsq<0] = 1.e-14
    neg_mapsq = -mapsq
    neg_mapsq[mapsq>0] = 1.e-14
    pos_mxsq = mxsq.copy()
    pos_mxsq[mxsq<0] = 1.e-14
    neg_mxsq = -mxsq
    neg_mxsq[mxsq>0] = 1.e-14

    sig = numpy.sqrt(varmapsq)

    ax.plot(meanr, pos_mapsq, color='blue')
    ax.plot(meanr, neg_mapsq, color='blue', ls=':')
    ax.errorbar(meanr[mapsq>0], mapsq[mapsq>0], yerr=sig[mapsq>0], color='blue', ls='')
    ax.errorbar(meanr[mapsq<0], -mapsq[mapsq<0], yerr=sig[mapsq<0], color='blue', ls='')
    mapsq_line, = ax.plot(-meanr, mapsq, color='blue')
    ax.set_ylim( [1.e-8, 6.e-6] )
    ax.set_ylabel(r'$M_{\rm ap}(\theta)$')
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log')

    ax.plot(meanr, pos_mxsq, color='red')
    ax.plot(meanr, neg_mxsq, color='red', ls=':')
    ax.errorbar(meanr[mxsq>0], mxsq[mxsq>0], yerr=sig[mxsq>0], color='red', ls='')
    ax.errorbar(meanr[mxsq<0], -mxsq[mxsq<0], yerr=sig[mxsq<0], color='red', ls='')
    mxsq_line, = ax.plot(-meanr, mxsq, color='red')
    ax.legend([mapsq_line, mxsq_line],
                [r'$M_{\rm ap}(\theta)$', r'$M_X(\theta)$'],
                loc='upper right', fontsize=24)

    ax.set_xlim( [1.,100.] )
    ax.set_xlabel(r'$\theta$ (arcmin)')

    plt.tight_layout()
    #plt.set_size_inches(8., 12.)
    plt.savefig(filename)


def load_data():
    imcat = pyfits.open('des_sv_wl_im3shape.fits')[1].data
    print 'imcat has %d objects'%len(imcat)
    ngcat = pyfits.open('des_sv_wl_ngmix.fits')[1].data
    print 'ngcat has %d objects'%len(ngcat)
    fcat = pyfits.open('des_sv_wl_info.fits')[1].data
    print 'gcat has %d objects'%len(fcat)

    gold = fcat['sva1_gold_flags'] == 0
    spte = fcat['sva1_spte_flags'] == 0
    mag = fcat['sva1_gold_mag_flags'] == 0
    ngmix = fcat['ngmix_flags'] == 0
    im3shape = fcat['im3shape_flags'] == 0

    all_ng = gold & spte & mag & ngmix
    all_im = gold & spte & mag & im3shape
    print 'ng selection includes %d galaxies'%numpy.sum(all_ng)
    print 'im selection includes %d galaxies'%numpy.sum(all_im)

    ngmask = numpy.in1d(ngcat['coadd_objects_id'], fcat['coadd_objects_id'][all_ng])
    immask = numpy.in1d(imcat['coadd_objects_id'], fcat['coadd_objects_id'][all_im])

    return fcat, ngcat[ngmask], imcat[immask]

def match_cats(raw_fcat, raw_ngcat, raw_imcat):
    print 'len raw cats = ',len(raw_fcat), len(raw_ngcat), len(raw_imcat)

    gold = raw_fcat['sva1_gold_flags'] == 0
    spte = raw_fcat['sva1_spte_flags'] == 0
    mag = raw_fcat['sva1_gold_mag_flags'] == 0
    ngmix = raw_fcat['ngmix_flags'] == 0
    im3shape = raw_fcat['im3shape_flags'] == 0
    all_ng = gold & spte & mag & ngmix
    all_im = gold & spte & mag & im3shape
    match = all_ng & all_im
    print 'match includes %d galaxies'%numpy.sum(match)

    ngmask = numpy.in1d(raw_ngcat['coadd_objects_id'], raw_fcat['coadd_objects_id'][match])
    immask = numpy.in1d(raw_imcat['coadd_objects_id'], raw_fcat['coadd_objects_id'][match])

    fcat = raw_fcat[match]
    ngcat = raw_ngcat[ngmask]
    imcat = raw_imcat[immask]
    print 'len matched cats = ',len(fcat), len(ngcat),len(imcat)

    return fcat, ngcat, imcat

def linear_fit(y, x, w=None, mask=None):
    import numpy
    # Follow along at http://mathworld.wolfram.com/LeastSquaresFitting.html
    n = len(x)
    if w is None:
        w = numpy.ones(n)
    if mask is not None:
        x = x[mask]
        y = y[mask]
        w = w[mask]
        n = len(x)
    xm = numpy.sum(w * x) / numpy.sum(w)
    ym = numpy.sum(w * y) / numpy.sum(w)
    ssxx = numpy.sum( w * numpy.abs(x-xm)**2 )
    ssyy = numpy.sum( w * numpy.abs(y-ym)**2 )
    ssxy = numpy.sum( w * numpy.conjugate(x-xm)*(y-ym) )
    m = ssxy / ssxx
    c = ym - m * xm
    s = numpy.sqrt( (ssyy - numpy.abs(m*ssxy)) / (n-2) )
    sigm = s / numpy.sqrt(ssxx)
    sigc = s * numpy.sqrt( 1./n + numpy.abs(xm)**2/ssxx )

    return m, c, sigm, sigc

def calculate_alpha(cat):
    if 'exp_e_1' in cat.columns.names:
        e = cat['exp_e_1'] + 1j*cat['exp_e_2']
        epsf = cat['psfrec_e_1'] + 1j*cat['psfrec_e_1']
        w = cat['exp_w']
    else:
        e = cat['e1'] - cat['nbc_c1'] + 1j*(cat['e2'] - cat['nbc_c2'])
        epsf = cat['mean_psf_e1_sky'] + 1j*cat['mean_psf_e2_sky']
        w = cat['w']

    alpha = linear_fit(e, epsf, w)
    print alpha[0:4:2]
    return alpha[0]


def calculate_mean(cat):
    if 'exp_e_1' in cat.columns.names:
        e = cat['exp_e_1'] + 1j*cat['exp_e_2']
        w = cat['exp_w']
        #s = cat['exp_e_sens_avg']
    else:
        e = cat['e1'] - cat['nbc_c1'] + 1j*(cat['e2'] - cat['nbc_c2'])
        w = cat['w']
        #s = cat['nbc_m'] + 1.

    #return numpy.sum(w*e) / numpy.sum(w*s)
    return numpy.sum(w*e) / numpy.sum(w)

def compute_xi(ra, dec, e1, e2, w, s):

    # Apply the mean sensitivity in each case
    #means = numpy.sum(w*s) / numpy.sum(w)
    #e1c = e1 / means
    #e2c = e2 / means

    # Build a TreeCorr catalog
    cat = treecorr.Catalog(ra=ra, dec=dec, g1=e1, g2=e2, k=s, w=w, ra_units='deg', dec_units='deg')

    # Compute the correlation function
    gg = treecorr.GGCorrelation(bin_size=0.05, min_sep=0.1, max_sep=500, sep_units='arcmin', output_dots=True, verbose=2)
    ss = treecorr.KKCorrelation(bin_size=0.05, min_sep=0.1, max_sep=500, sep_units='arcmin', output_dots=True, verbose=2)

    print 'Start corr2'
    gg.process(cat)
    ss.process(cat)

    gg.xip /= ss.xi
    gg.xim /= ss.xi
    gg.varxi /= ss.xi**2

    return gg


def compute_xi_ng(fcat, ngcat, mask=None):
    # Get the RA, Dec
    ra = fcat['ra']
    dec = fcat['dec']
    ng = numpy.in1d(fcat['coadd_objects_id'], ngcat['coadd_objects_id'])

    # Get the shear data
    nge1 = ngcat['exp_e_1']
    nge2 = ngcat['exp_e_2']
    ngw = ngcat['exp_w']
    ngs = ngcat['exp_e_sens_avg']

    if mask is None:
        return compute_xi(ra[ng], dec[ng], nge1, nge2, ngw, ngs)
    else:
        return compute_xi(ra[ng][mask], dec[ng][mask], nge1[mask], nge2[mask], ngw[mask], ngs[mask])


def compute_xi_im(fcat, imcat, mask=None):
    # Get the RA, Dec
    ra = fcat['ra']
    dec = fcat['dec']
    im = numpy.in1d(fcat['coadd_objects_id'], imcat['coadd_objects_id'])

    # Get the shear data
    ime1 = imcat['e1'] - imcat['nbc_c1']
    ime2 = imcat['e2'] - imcat['nbc_c2']
    imw = imcat['w']
    numpy.clip(imw, 0.0, 0.24**-2, imw)
    ims = imcat['nbc_m'] + 1.

    if mask is None:
        return compute_xi(ra[im], dec[im], ime1, ime2, imw, ims)
    else:
        return compute_xi(ra[im][mask], dec[im][mask], ime1[mask], ime2[mask], imw[mask], ims[mask])



def compute_dxi(fcat, ngcat, imcat, mask=None, single=False, raw_ngcat=None, raw_imcat=None):
    # Get the RA, Dec
    ra = fcat['ra']
    dec = fcat['dec']

    # Get the shear data
    nge1 = ngcat['exp_e_1']
    nge2 = ngcat['exp_e_2']
    ngw = ngcat['exp_w']
    ngs = ngcat['exp_e_sens_avg']

    ime1 = imcat['e1'] - imcat['nbc_c1']
    ime2 = imcat['e2'] - imcat['nbc_c2']
    imw = imcat['w']
    numpy.clip(imw, 0.0, 0.24**-2, imw)
    ims = imcat['nbc_m'] + 1.

    # Use the geometric mean for the weight.
    w = numpy.sqrt(imw * ngw)

    meanngs = numpy.sum(w*ngs) / numpy.sum(w)
    meanims = numpy.sum(w*ims) / numpy.sum(w)

    if raw_imcat is not None and raw_ngcat is not None:
        im_alpha1 = calculate_alpha(raw_imcat)
        im_alpha2 = calculate_alpha(imcat)
        ng_alpha1 = calculate_alpha(raw_ngcat)
        ng_alpha2 = calculate_alpha(ngcat)
        print 'dalpha_ng = ',ng_alpha2-ng_alpha1
        print 'dalpha_im = ',im_alpha2-im_alpha1
        nge1 -= (ng_alpha2-ng_alpha1).real * ngcat['psfrec_e_1']
        nge2 -= (ng_alpha2-ng_alpha1).real * ngcat['psfrec_e_2']
        ime1 -= (im_alpha2-im_alpha1).real * imcat['mean_psf_e1_sky']
        ime2 -= (im_alpha2-im_alpha1).real * imcat['mean_psf_e2_sky']

        im_mean1 = calculate_mean(raw_imcat)
        im_mean2 = calculate_mean(imcat)
        ng_mean1 = calculate_mean(raw_ngcat)
        ng_mean2 = calculate_mean(ngcat)
        print 'dmean_ng = ',ng_mean2-ng_mean1
        print 'dmean_im = ',im_mean2-im_mean1
        nge1 -= (ng_mean2-ng_mean1).real
        nge2 -= (ng_mean2-ng_mean1).real
        ime1 -= (im_mean2-im_mean1).real
        ime2 -= (im_mean2-im_mean1).real

    # Apply the mean sensitivity in each case
    nge1c = nge1 / meanngs
    nge2c = nge2 / meanngs
    ime1c = ime1 / meanims
    ime2c = ime2 / meanims

    print 'mean ngmix shear (matched) = ',numpy.mean(nge1c),numpy.mean(nge2c)
    print 'weighted = ',numpy.sum(nge1c*ngw)/numpy.sum(ngw),numpy.sum(nge2c*ngw)/numpy.sum(ngw)
    print 'mean im3shape shear (matched) = ',numpy.mean(ime1c),numpy.mean(ime2c)
    print 'weighted = ',numpy.sum(ime1c*imw)/numpy.sum(imw),numpy.sum(ime2c*imw)/numpy.sum(imw)

    de1 = nge1c - ime1c
    de2 = nge2c - ime2c
    print 'mean diff = ',numpy.sum(de1*w)/numpy.sum(w),numpy.sum(de2*w)/numpy.sum(w)

    if False:
        # This is completely invalid!  Just a test to see if it would work.
        # Which it did, but the selection-bias alpha subtraction is right thing to do.
        print 'mean diff = ',numpy.sum(de1*w)/numpy.sum(w),numpy.sum(de2*w)/numpy.sum(w)
        de1 -= numpy.sum(de1 * w) / numpy.sum(w)
        de2 -= numpy.sum(de2 * w) / numpy.sum(w)

    # Build a TreeCorr catalog with the difference.
    if mask is None:
        cat = treecorr.Catalog(ra=ra, dec=dec, g1=de1, g2=de2, w=w, ra_units='deg', dec_units='deg')
    else:
        cat = treecorr.Catalog(ra=ra[mask], dec=dec[mask], g1=de1[mask], g2=de2[mask], w=w[mask], ra_units='deg', dec_units='deg')

    # Compute the correlation function
    gg = treecorr.GGCorrelation(bin_size=0.2, min_sep=0.5, max_sep=300, sep_units='arcmin', bin_slop=0.5, output_dots=True, verbose=2)

    print 'Start corr2'
    gg.process(cat)

    if single:
        # Also compute the single-catalog correlations
        cat_ng = treecorr.Catalog(ra=ra, dec=dec, g1=nge1c, g2=nge2c, w=w, ra_units='deg', dec_units='deg')
        cat_im = treecorr.Catalog(ra=ra, dec=dec, g1=ime1c, g2=ime2c, w=w, ra_units='deg', dec_units='deg')

        gg_ng = treecorr.GGCorrelation(bin_size=0.2, min_sep=0.5, max_sep=300, sep_units='arcmin', bin_slop=0.5, output_dots=True, verbose=2)
        gg_ng.process(cat_ng)

        gg_im = treecorr.GGCorrelation(bin_size=0.2, min_sep=0.5, max_sep=300, sep_units='arcmin', bin_slop=0.5, output_dots=True, verbose=2)
        gg_im.process(cat_im)

        gg_x = treecorr.GGCorrelation(bin_size=0.2, min_sep=0.5, max_sep=300, sep_units='arcmin', bin_slop=0.5, output_dots=True, verbose=2)
        gg_x.process(cat_ng, cat_im)

        return gg, gg_ng, gg_im, gg_x

    else:
        return gg

def main():
    raw_fcat, raw_ngcat, raw_imcat = load_data()

    # dede plot
    if True:
        single = False
        match_fcat, match_ngcat, match_imcat = match_cats(raw_fcat, raw_ngcat, raw_imcat)
        gg = compute_dxi(match_fcat, match_ngcat, match_imcat, raw_ngcat=raw_ngcat, raw_imcat=raw_imcat, single=single)
        if single:
            gg, gg_ng, gg_im, gg_x = gg
        plot_dxi(gg, 'xi_de.pdf')
        plot_dxi(gg, 'xi_de.png')

        if single:
            plot_dxi(gg_ng, 'xi_ng.png')
            plot_dxi(gg_im, 'xi_im.png')
            plot_dxi(gg_x, 'xi_x.png')
            gg_ng.xip -= gg_im.xip
            gg_ng.xim -= gg_im.xim
            plot_dxi(gg_ng, 'xi_diff.png')

    # B-mode plot
    if False:
        gg_ng = compute_xi_ng(raw_fcat, raw_ngcat, raw_imcat)


if __name__ == '__main__':
    main()

