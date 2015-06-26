import pyfits
import numpy
import treecorr

def plot_xi(gg, filename, sqrtn=1, label=None):
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
                    loc='upper right')
    ax.set_ylim( [5.e-8, 2.e-5] )
    if label is None:
        ax.set_ylabel(r'$\xi_{+,\Delta e}(\theta)$')
    else:
        ax.set_ylabel(label)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log')

    if False:
        #ax = axes[1]
        ax.fill( [theta[0]] + theta + [theta[-1]], [0.] + dxim + [0.], color = 'Lime')
        ax.plot(meanr, pos_xim, color='blue')
        ax.plot(meanr, neg_xim, color='blue', ls=':')
        ax.errorbar(meanr[xim>0], xim[xim>0], yerr=sig[xim>0]/sqrtn, color='green', ls='')
        ax.errorbar(meanr[xim<0], -xim[xim<0], yerr=sig[xim<0]/sqrtn, color='green', ls='')
        xim_line, = ax.plot(-meanr, xim, color='green')
        #xim_req = mp.Patch(color='#FFFF82')
        xim_req = mp.Patch(color='Lime')
        #ax.legend([xim_line, xim_req],
                    #[r'$\xi_-(\theta)$', r'$\xi_-$ requirement'],
                    #loc='upper right')
        #ax.set_ylim( [1.e-8, 2.e-6] )
        ax.set_ylim( [1.e-8, 2.e-5] )
        ax.set_ylabel(r'$\xi_{+,-,\Delta e}(\theta)$')
        ax.set_yscale('log', nonposy='clip')
        ax.set_xscale('log')

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
    raw_ngcat = ngcat[ngmask]
    raw_imcat = imcat[immask]
    raw_fcat = fcat
    print 'len raw cats = ',len(raw_fcat), len(raw_ngcat), len(raw_imcat)

    match = all_ng & all_im
    print 'match includes %d galaxies'%numpy.sum(match)

    ngmask = numpy.in1d(ngcat['coadd_objects_id'], fcat['coadd_objects_id'][match])
    immask = numpy.in1d(imcat['coadd_objects_id'], fcat['coadd_objects_id'][match])

    fcat = fcat[match]
    ngcat = ngcat[ngmask]
    imcat = imcat[immask]
    print 'len matched cats = ',len(fcat), len(ngcat),len(imcat)

    return fcat, ngcat, imcat, raw_fcat, raw_ngcat, raw_imcat

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

    meanngs = numpy.sum(w*ngs) / numpy.sum(w)
    meanims = numpy.sum(w*ims) / numpy.sum(w)

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
    fcat, ngcat, imcat, raw_fcat, raw_ngcat, raw_imcat = load_data()
    gg = compute_dxi(fcat, ngcat, imcat, raw_ngcat=raw_ngcat, raw_imcat=raw_imcat)
    # gg, gg_ng, gg_im, gg_x = dxi.compute_dxi(fcat, ngcat, imcat, ngsnr > 20, single=True)
    plot_xi(gg, 'xip_de.pdf')


if __name__ == '__main__':
    main()

