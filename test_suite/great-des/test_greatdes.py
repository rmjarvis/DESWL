#import matplotlib
#matplotlib.use('Agg') # needs to be done before import pyplot

# Command from Troxel to look into later:
# plt.hist2d(svi3.rsnr,svi3.radius,bins=500,range=((0,100),(1,4)),norm=LogNorm())

def load_im_data():
    import pyfits
    import glob
    import numpy
    #files = sorted(glob.glob('results/im3shape/results-disc/main_cats/nbc.meds.*.fits'))
    files = sorted(glob.glob('results/im3shape/results_bord_rsnr/nbc.meds.*.rsnr.fits'))
    num = []
    rad = []
    rgp = []
    flux = []
    snr = []
    flag = []
    e = []
    round_snr = []
    is_disc = []
    niter = []
    chisq = []
    min_res = []
    max_res = []
    info = []
    psf_e = []
    psf_fwhm = []
    trcov = []
    detcov = []
    var_r = []
    var_e = []

    # Sort first by gnum then fnum
    fnum = [ int(f.split('.')[2]) for f in files ]
    gnum = [ int(f.split('.')[3][1:]) for f in files ]
    flist = sorted( zip(gnum, fnum, files) )

    for gnum, fnum, f in flist:
        #tokens = f.split('.')
        #fnum = int(tokens[2])
        #gnum = int(tokens[3][1:])
        with pyfits.open(f) as pyf:
            data = pyf[1].data
            print f, len(data)
            xnum = data['coadd_objects_id'].copy()
            xnum[:] += fnum * 10000
            xnum[:] += gnum * 10000000
            num.append(xnum)
            rad.append(data['radius'].copy())
            rgp.append(data['mean_rgpp_rp'].copy())
            flux.append(data['mean_flux'].copy())
            snr.append(data['snr'].copy())
            flag.append(data['error_flag'].copy())
            # Weird.  e1 seems to be flipped!
            e1 = -data['e1']
            e2 = data['e2']
            e.append( e1 + 1j * e2 )
            if 'round_snr' in data.columns.names:
                round_snr.append(data['round_snr'].copy())
            is_disc.append( (data['bulge_flux'] == 0.) )
            niter.append(data['levmar_iterations'].copy())
            chisq.append(data['chi2_pixel'].copy())
            min_res.append(data['min_residuals'].copy())
            max_res.append(data['max_residuals'].copy())
            info.append(data['info_flag'].copy())
            psf_e.append( -data['mean_psf_e1_sky'] + 1j * data['mean_psf_e2_sky'] )
            psf_fwhm.append(data['mean_psf_fwhm'].copy())
            c22 = data['covmat_2_2']
            c23 = -data['covmat_2_3']
            c33 = data['covmat_3_3']
            c44 = data['covmat_4_4']
            radius = data['radius']
            esq = e1**2 + e2**2
            trcov.append(c22 + c33)
            detcov.append(c22 * c33 - c23**2)
            var_r.append(c44 / radius**2)
            ve = (e1**2 * c22 + 2.*e1*e2 * c23 + e2**2 * c33) / esq
            ve[esq == 0.] = (c22 + c33)[esq==0.]/2.
            var_e.append(ve)

    num = numpy.concatenate(num)
    rad = numpy.concatenate(rad)
    rgp = numpy.concatenate(rgp)
    flux = numpy.concatenate(flux)
    snr = numpy.concatenate(snr)
    flag = numpy.concatenate(flag)
    e = numpy.concatenate(e)
    if len(round_snr) > 0:
        round_snr = numpy.concatenate(round_snr)
    is_disc = numpy.concatenate(is_disc)
    niter = numpy.concatenate(niter)
    chisq = numpy.concatenate(chisq)
    min_res = numpy.concatenate(min_res)
    max_res = numpy.concatenate(max_res)
    info = numpy.concatenate(info)
    psf_e = numpy.concatenate(psf_e)
    psf_fwhm = numpy.concatenate(psf_fwhm)
    trcov = numpy.concatenate(trcov)
    detcov = numpy.concatenate(detcov)
    var_r = numpy.concatenate(var_r)
    var_e = numpy.concatenate(var_e)

    return num, rad, rgp, flux, snr, flag, e, round_snr, is_disc, niter, chisq, min_res, max_res, info, psf_e, psf_fwhm, trcov, detcov, var_r, var_e



def load_im_sv():
    import pyfits
    import glob
    import numpy
    #files = sorted(glob.glob('results/im3shape/results-disc/main_cats/nbc.meds.*.fits'))
    files = sorted(glob.glob('/astro/u/astrodat/data/DES/wlpipe/im3shape_v9/full_cats/*.fits'))

    rad = []
    rgp = []
    snr = []
    flag = []
    e = []
    is_disc = []
    info = []

    for f in files:
        with pyfits.open(f) as pyf:
            data = pyf[1].data
            print f, len(data)
            rad.append(data['radius'].copy())
            rgp.append(data['mean_rgpp_rp'].copy())
            snr.append(data['snr'].copy())
            flag.append(data['error_flag'].copy())
            e.append(data['e1'] + 1j * data['e2'])
            is_disc.append( (data['bulge_flux'] == 0.) )
            info.append(data['info_flag'].copy())

    rad = numpy.concatenate(rad)
    rgp = numpy.concatenate(rgp)
    snr = numpy.concatenate(snr)
    flag = numpy.concatenate(flag)
    e = numpy.concatenate(e)
    is_disc = numpy.concatenate(is_disc)
    info = numpy.concatenate(info)

    return rad, rgp, snr, flag, e, is_disc, info


def load_ng_data():
    import pyfits
    import glob
    import numpy

    files = sorted(glob.glob('/gpfs01/astro/workarea/esheldon/lensing/great-des/sfit-e02/collated/sfit-e02-*.fits'))

    num = []
    t = []
    flux = []
    snr = []
    tsnr = []
    flag = []
    t_r = []
    snr_r = []
    flag_r = []
    e = []
    e_psf = []
    t_psf = []
    sens = []
    chisq = []
    trcov = []

    for f in files:
        tokens = f.split('-')
        gnum = int(tokens[-1][1:3])
        with pyfits.open(f) as pyf:
            data = pyf[1].data
            print f, len(data)
            fnum = data['fnum']
            xnum = data['number'].copy()
            xnum[:] += fnum * 10000
            xnum[:] += gnum * 10000000
            num.append(xnum)
            t.append(numpy.exp(data['log_T']))
            flux.append(numpy.exp(data['log_flux']))
            snr.append(data['s2n_w'].copy())
            tsnr.append(data['T_s2n'].copy())
            flag.append(data['flags'] | data['psf_flags'])
            t_r.append(numpy.exp(data['log_T_r']))
            snr_r.append(data['s2n_r'].copy())
            flag_r.append(data['flags_r'])
            e1 = data['g'][:,0]
            e2 = data['g'][:,1]
            e.append( e1 + 1j * e2 )
            e_psf.append(data['psf_g'][:,0] + 1j * data['psf_g'][:,1])
            t_psf.append(data['psf_T'].copy())
            sens.append(data['g_sens'][:,0] + 1j * data['g_sens'][:,1])
            chisq.append(data['chi2per'])
            c22 = data['pars_cov'][:,2,2]
            c33 = data['pars_cov'][:,3,3]
            trcov.append(c22 + c33)

    num = numpy.concatenate(num)
    t = numpy.concatenate(t)
    flux = numpy.concatenate(flux)
    snr = numpy.concatenate(snr)
    tsnr = numpy.concatenate(tsnr)
    flag = numpy.concatenate(flag)
    t_r = numpy.concatenate(t_r)
    snr_r = numpy.concatenate(snr_r)
    flag_r = numpy.concatenate(flag_r)
    e = numpy.concatenate(e)
    e_psf = numpy.concatenate(e_psf)
    t_psf = numpy.concatenate(t_psf)
    sens = numpy.concatenate(sens)
    chisq = numpy.concatenate(chisq)
    trcov = numpy.concatenate(trcov)

    return num, t, flux, snr, tsnr, flag, t_r, snr_r, flag_r, e, e_psf, t_psf, sens, chisq, trcov

def load_truth():
    import pyfits
    import numpy
    import os
    import glob

    files = sorted(glob.glob('data/nbc.truth.*.fits'))

    num = []
    rad = []
    flux = []
    e = []
    rawe = []
    g_app = []
    id = []
    srcn = []
    z = []
    use = []

    # Sort first by gnum then fnum
    fnum = [ int(f.split('.')[2]) for f in files ]
    gnum = [ int(f.split('.')[3][1:]) for f in files ]
    flist = sorted( zip(gnum, fnum, files) )

    for gnum, fnum, f in flist:
        #tokens = f.split('.')
        #fnum = int(tokens[2])
        #gnum = int(tokens[3][1:])
        with pyfits.open(f) as pyf:
            data = pyf[1].data
            print f, len(data)
            xnum = data['id'].copy()
            xnum[:] += fnum * 10000
            xnum[:] += gnum * 10000000
            assert all(data['id_shear'] == gnum)
            num.append(xnum)
            rad.append(data['sf_hlr'].copy())
            flux.append(data['flux'].copy())
            orig_e = data['shape_e1'] + 1j * data['shape_e2']
            orig_g = orig_e / (1. + numpy.sqrt(1. - numpy.abs(orig_e)**2))
            rot_g = orig_g * numpy.exp(2j * data['rotation_angle'])
            g = data['g1_true'] + 1j * data['g2_true']
            g_final = (rot_g + g) / (1. + g.conjugate() * rot_g)
            e.append(g_final)
            g_app.append(g)
            rawe.append(numpy.abs(orig_g))
            id.append(data['id_cosmos'].copy())
            srcn.append(data['sf_sersicn'].copy())
            z.append(data['zphot'].copy())
            use.append(data['to_use'] == 1)

    num = numpy.concatenate(num)
    rad = numpy.concatenate(rad)
    flux = numpy.concatenate(flux)
    e = numpy.concatenate(e)
    rawe = numpy.concatenate(rawe)
    g_app = numpy.concatenate(g_app)
    id = numpy.concatenate(id)
    srcn = numpy.concatenate(srcn)
    z = numpy.concatenate(z)
    use = numpy.concatenate(use)

    return num, rad, flux, e, rawe, g_app, id, srcn, z, use


def load_great3(id):
    import pyfits
    with pyfits.open('cosmos_in_great3.fits') as f:
        cosmos_id = f[1].data['cosmos_ident']
        used = f[1].data['used_in_great3']
    with pyfits.open('/astro/u/mjarvis/share/galsim/COSMOS_23.5_training_sample/real_galaxy_catalog_23.5.fits') as f:
        cosmos_id = f[1].data['ident']
        weight = f[1].data['weight']

    return used[id] == 1 , weight[id]


def check_meas(r_true, r_meas, e_true, e_meas, mask, c=None, filename=None):
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    print len(r_true)
    print len(r_true[mask])
    mask2 = mask & (r_meas < 2) & (r_true < 2)
    print len(r_true[mask2])
    if c is not None: 
        print 'c = ',c
        c = c[mask2]
        print 'c => ',c
    if len(r_true[mask2]) <= 100000:
        fig, ax = plt.subplots(1, 2)
        ax[0].scatter(r_meas[mask2], r_true[mask2], c=c, s=0.01, marker='+')
        ax[1].scatter(e_meas[mask2].real, e_true[mask2].real, c=c, s=0.01, marker='+')
        ax[1].scatter(e_meas[mask2].imag, e_true[mask2].imag, c=c, s=0.01, marker='+')
        fig.set_size_inches(8.5,4.0)
        ax[1].set_xlabel(r'Measured e1,e2')
    else:
        fig, ax = plt.subplots(1, 3)
        ax[0].hexbin(r_meas[mask2], r_true[mask2])
        ax[1].hexbin(e_meas[mask2].real, e_true[mask2].real)
        ax[2].hexbin(e_meas[mask2].imag, e_true[mask2].imag)
        fig.set_size_inches(12.5,4.0)
        ax[1].set_xlabel(r'Measured e1')
        ax[2].set_xlabel(r'Measured e2')
        ax[2].set_xlim(-1.0,1.0)
        ax[2].set_ylim(-1.0,1.0)
    ax[0].set_xlabel(r'Measured radius')
    ax[0].set_ylabel(r'"True" value')
    #ax[0].set_ylim(0,2)
    ax[0].set_xlim(0,2)
    ax[1].set_xlim(-1.0,1.0)
    ax[1].set_ylim(-1.0,1.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def check_snr(snr_meas, snr_round, f_meas, r_meas, mask, filename=None):
    import matplotlib.pyplot as plt
    import numpy
    plt.style.use('supermongo')
    if False:
        fig, ax = plt.subplots(4, 4, sharex='col', sharey='row')
        snr_alt = f_meas / r_meas
        mask2 = mask & (snr_meas < 50) & (snr_round < 50) & (r_meas < 2) & (f_meas < 10000)
        ax[3,0].hexbin(snr_meas[mask2], snr_round[mask2])
        ax[2,0].hexbin(snr_meas[mask2], snr_alt[mask2])
        ax[1,0].hexbin(snr_meas[mask2], f_meas[mask2])
        ax[0,0].hexbin(snr_meas[mask2], r_meas[mask2])
        ax[3,1].hexbin(r_meas[mask2], snr_round[mask2])
        ax[2,1].hexbin(r_meas[mask2], snr_alt[mask2])
        ax[1,1].hexbin(r_meas[mask2], f_meas[mask2])
        ax[3,2].hexbin(f_meas[mask2], snr_round[mask2])
        ax[2,2].hexbin(f_meas[mask2], snr_alt[mask2])
        ax[3,3].hexbin(snr_alt[mask2], snr_round[mask2])
        ax[3,0].set_xlabel(r'snr_w')
        ax[3,1].set_xlabel(r'radius')
        ax[3,2].set_xlabel(r'flux')
        ax[3,3].set_xlabel(r'flux / radius')
        ax[3,0].set_ylabel(r'snr_r')
        ax[2,0].set_ylabel(r'flux / radius')
        ax[1,0].set_ylabel(r'flux')
        ax[0,0].set_ylabel(r'radius')
        ax[3,0].set_xlim(0,50)
        ax[3,1].set_xlim(0,2)
        ax[3,2].set_xlim(0,10000)
        ax[3,3].set_xlim(0,10000)
        ax[3,0].set_ylim(0,50)
        ax[2,0].set_ylim(0,10000)
        ax[1,0].set_ylim(0,10000)
        ax[0,0].set_ylim(0,2)
    else:
        fig, ax = plt.subplots(1, 2)
        mask2 = mask & (snr_meas < 50) & (snr_round < 50) 
        ax[0].hexbin(snr_meas[mask2], snr_round[mask2])
        ax[1].hexbin(numpy.log10(snr_meas[mask]), numpy.log10(snr_round[mask]))
        ax[0].set_xlabel(r'snr_w')
        ax[0].set_ylabel(r'snr_r')
        ax[1].set_xlabel(r'log10(snr_w)')
        ax[1].set_ylabel(r'log10(snr_r)')
    fig.set_size_inches(8.5,4.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def compute_empir_cov(e_meas, e_true, g_app, id, mask):
    import numpy 

    id_use = []
    trcov_emp = []
    mean_e = []
    for i in range(numpy.max(id)):
        select = numpy.where((id == i) & mask)
        if numpy.sum(select) < 100:
            continue
        et = e_true[select]
        em = e_meas[select]
        g = g_app[select]
        er = (et-g)/(1.-g.conjugate()*et)
        if numpy.mean(er) == 0.0:
            continue

        id_use.append(i)

        # Undo the applied shear
        em = (em-g)/(1.-g.conjugate()*em)
        # Align shear such that the true shear is a pure positive e1 value.
        em *= er.conjugate()

        c11 = numpy.var(em.real)
        c22 = numpy.var(em.imag)
        print i, numpy.mean(em), c11, c22

        trcov_emp.append(c11+c22)
        mean_e.append(numpy.abs(numpy.mean(em.real)))

    trcov_emp = numpy.array(trcov_emp)
    mean_e = numpy.array(mean_e)
    id_use = numpy.array(id_use)

    return trcov_emp, mean_e, id_use

def compute_trcov_meas(id_use, trcov, e, var_r, snr_meas, snr_round, id, mask):
    import numpy
    ncov = 8
    trcov_meas = [ [] for i in range(ncov) ]
    for i in id_use:
        print i
        select = numpy.where((id == i) & mask)

        factor = numpy.sqrt(1.-numpy.abs(e[select])**2)
        trcov_meas[0].append(numpy.mean(trcov[select]))
        trcov_meas[1].append(numpy.mean(trcov[select]/factor))
        trcov_meas[2].append(numpy.mean(trcov[select]*factor))
        trcov_meas[3].append(numpy.mean(trcov[select]/factor**2))
        trcov_meas[4].append(numpy.mean(trcov[select]/factor**4))
        trcov_meas[5].append(numpy.mean(var_r[select]))
        trcov_meas[6].append(numpy.mean(1./snr_meas[select]**2))
        trcov_meas[7].append(numpy.mean(1./snr_round[select]**2))

    for i in range(ncov):
        trcov_meas[i] = numpy.array(trcov_meas[i])

    return trcov_meas


def check_cov2(trcov_emp, mean_e, id_use, trcov_meas, filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')

    ncov = len(trcov_meas)
    labels = [
        r'$Tr(Cov)$',
        r'$Tr(Cov)/\sqrt{1-|e|^2}$',
        r'$Tr(Cov)*\sqrt{1-|e|^2}$',
        r'$Tr(Cov)/(1-|e|^2)$',
        r'$Tr(Cov)/(1-|e|^2)^2$',
        r'$Var(r)/r^2$',
        r'$(S/N)_w^{-2}$',
        r'$(S/N)_r^{-2}$',
        ]

    fig, axes = plt.subplots(2, ncov)
    axes[0,0].set_ylabel(r'Empirical $Tr(Cov(e1,e2)) \times (S/N)_w^2$')
    axes[1,0].set_ylabel(r'Empirical $log_{10}(Tr(Cov(e1,e2)))$')
    tc = trcov_meas[0]
    snrsq = trcov_meas[6]**-1
    for i in range(ncov):
        if i == 6:
            continue
        print labels[i]
        ax = axes[0,i]
        ax.hist2d(trcov_meas[i]*snrsq, trcov_emp*snrsq, bins=500)
        #ax.hexbin(trcov_meas[i][tc<0.01], trcov_emp[tc<0.01])
        C = numpy.corrcoef(trcov_meas[i]*snrsq, trcov_emp*snrsq)[0,1]
        print 'direct C = ',C
        ax.set_xlabel(labels[i]+r' $\times (S/N)_w^2$')
        ax.set_title('Corr coeff is %.3f'%C)
        ax = axes[1,i]
        ax.hist2d(numpy.log10(trcov_meas[i][tc<1]), numpy.log10(trcov_emp[tc<1]), bins=500)
        #ax.hexbin(trcov_meas[i] , trcov_emp, xscale='log', yscale='log')
        ax.set_xlabel('log10('+labels[i]+')')
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        C = numpy.corrcoef(numpy.log(trcov_meas[i]), numpy.log(trcov_emp))[0,1]
        print 'C for logs = ',C
        ax.set_title('Corr coeff = %.3f'%C)

    #ax.scatter(mean_e, trcov_meas / trcov_emp)
    #ax.set_xlabel(r'$|e|$')
    #ax.set_ylabel(r'TrCov ratio im3shape / Empirical')
    fig.set_size_inches(ncov * 6.0, 12.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)



def check_cov(e_meas, e_true, g_app, trcov, id, mask, filename=None, labels=['Tr(Cov)']):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')

    if not isinstance(trcov,list):
        trcov = [ trcov ]
    ncov = len(trcov)
    if len(labels) != ncov:
        raise ValueError("Wrong number of labels")

    trcov_emp = []
    trcov_meas = [ [] for i in range(ncov) ]
    mean_e = []
    for i in range(numpy.max(id)):
    #for i in range(1000):
        select = (id == i) & mask
        if numpy.sum(select) < 100:
            continue
        et = e_true[select]
        em = e_meas[select]
        g = g_app[select]
        er = (et-g)/(1.-g.conjugate()*et)
        if numpy.mean(er) == 0.0:
            continue

        # Undo the applied shear
        em = (em-g)/(1.-g.conjugate()*em)
        # Align shear such that the true shear is a pure positive e1 value.
        em *= er.conjugate()

        c11 = numpy.var(em.real)
        c22 = numpy.var(em.imag)
        print i, numpy.mean(em), c11, c22

        trcov_emp.append(c11+c22)
        mean_e.append(numpy.abs(numpy.mean(em.real)))
        for i in range(ncov):
            trcov_meas[i].append(numpy.mean(trcov[i][select]))

    trcov_emp = numpy.array(trcov_emp)
    mean_e = numpy.array(mean_e)
    for i in range(ncov):
        trcov_meas[i] = numpy.array(trcov_meas[i])

    fig, axes = plt.subplots(1, ncov, sharex='col', sharey='row')
    axes[0].set_ylabel(r'Empirical TrCov from scatter')
    for i in range(ncov):
        ax = axes[i]
        ax.hexbin(trcov_meas[i] , trcov_emp, xscale='log', yscale='log')
        ax.set_xlabel(labels[i])
        ax.set_xscale('log')
        ax.set_yscale('log')
        C = numpy.corrcoef(numpy.log(trcov_meas[i]), numpy.log(trcov_emp))[0,1]
        ax.set_title('Correlation coeff is %.3f'%C)

    #ax.scatter(mean_e, trcov_meas / trcov_emp)
    #ax.set_xlabel(r'$|e|$')
    #ax.set_ylabel(r'TrCov ratio im3shape / Empirical')
    fig.set_size_inches(ncov * 8.0, 8.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def dist_z(x, y, z, mask, xlabel, ylabel, filename=None, x_corr=None, y_corr=None):
    import matplotlib.pyplot as plt
    import numpy

    #fig, axes = plt.subplots(2, 2, sharex='col', sharey='row')
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)

    zbins = [ 0.3, 0.644, 0.901, 1.3 ]

    for ax, zmin, zmax in ( (axes[0,0], 0.3, 0.644),
                            (axes[0,1], 0.644, 0.901),
                            (axes[1,0], 0.901, 1.3),
                            (axes[1,1], 0.3, 1.3) ):

        mask2 = mask & (z >= zmin) & (z < zmax)
        xx = x[mask2]
        if x_corr is not None:
            xx /= numpy.mean(x_corr[mask2])
        if y is None:
            ax.hist(xx, bins=500, histtype='step', color='black', normed=1,
                    facecolor='red', fill=True)
        else:
            yy = y[mask2]
            if y_corr is not None:
                yy /= numpy.mean(y_corr[mask2])
            ax.hist2d(xx, yy, bins=500)
        ax.set_title('$%f < z < %f$'%(zmin,zmax))
        #ax.plot( (-1,1),(-1,1), color='w')

    axes[1,0].set_xlabel(xlabel)
    axes[1,1].set_xlabel(xlabel)
    if y is None:
        axes[0,0].set_ylabel('relative frequency')
        axes[1,0].set_ylabel('relative frequency')
    else:
        axes[0,0].set_ylabel(ylabel)
        axes[1,0].set_ylabel(ylabel)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def misc_dist(r, rgp, snr_m, snr_r, e, g, mask, filename=None):
    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams.update({'font.size': 6})
    fig, ax = plt.subplots(3, 2)

    zbins = [ 0.3, 0.644, 0.901, 1.3 ]

    ax[0,0].scatter(rgp[mask], snr_r[mask], s=0.1)
    ax[0,0].set_xlabel('rgpp_rp')
    ax[0,0].set_ylabel('snr_r')

    ax[0,1].scatter(rgp[mask], snr_m[mask], s=0.1)
    ax[0,1].set_xlabel('rgpp_rp')
    ax[0,1].set_ylabel('snr_m')

    ax[1,0].scatter(e[mask].real, e[mask].imag, s=0.1)
    ax[1,0].set_xlabel('e1')
    ax[1,0].set_ylabel('e2')

    ax[1,1].scatter(rgp[mask], r[mask], s=0.1)
    ax[1,1].set_xlabel('rgp')
    ax[1,1].set_ylabel('radius')

    ax[2,0].scatter(r[mask], snr_r[mask], s=0.1)
    ax[2,0].set_xlabel('radius')
    ax[2,0].set_ylabel('snr_r')

    ax[2,1].scatter(r[mask], snr_m[mask], s=0.1)
    ax[2,1].set_xlabel('radius')
    ax[2,1].set_ylabel('snr_m')


    fig.suptitle('Galaxies in high z bin with $S/N_r > 15$ but $S/N_w < 15$', fontsize=12)
    fig.set_size_inches(8.0,12.0)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

 
def scatter_r(r_true, r_meas, mask, filename=None):
    import matplotlib.pyplot as plt
    plt.clf()
    plt.scatter(r_meas[mask], r_true[mask])
    plt.ylim(0,3)
    plt.xlim(0,2)
    plt.xlabel('Measured radius')
    plt.ylabel('True radius')
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def scatter_e_r(r, e, mask, c=None, title='size vs shape', filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    fig, ax = plt.subplots(1, 1)
    abse = numpy.abs(e[mask])
    if c is None:
        ax.scatter(abse, r[mask], color='b', s=0.01)
    else:
        ax.scatter(abse, r[mask], c=c[mask], marker='+', s=0.01)
    ax.set_ylabel('measured radius')
    ax.set_xlabel(r'$|e_{true}|$')
    ax.set_ylim(0.,2.)
    ax.set_xlim(0.,1.)
    s = numpy.argsort(abse)
    nn = 300
    nbins = len(s) / nn
    meanr = numpy.array([ numpy.mean(r[mask][s][nn*i:nn*(i+1)]) for i in range(nbins) ])
    meane = numpy.array([ numpy.mean(abse[s][nn*i:nn*(i+1)]) for i in range(nbins) ])
    ax.plot(meane,meanr, color='k')
    #print 'meanr = ',meanr
    #print 'meane = ',meane
    ax.set_title(title)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def scatter_e(e, mask, c=None, title='e1,e2', filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    plt.clf()
    if c is None:
        plt.hexbin(e[mask].real, e[mask].imag)
    else:
        plt.scatter(e[mask].real, e[mask].imag, c=c[mask], s=0.1, marker='+')
    plt.xlabel('e1')
    plt.ylabel('e2')
    if c is not None:
        sortc = numpy.sort(c)
        deciles = numpy.array([ sortc[i*len(sortc)/10] for i in range(10) ] + [ sortc[-1] ])
        print deciles
        index = numpy.digitize(c[mask], deciles)
        meane = numpy.array([ numpy.mean(e[mask][index==i]) for i in range(10) ])
        plt.scatter(meane.real, meane.imag, c=(deciles[:-1] + deciles[1:])/2.)
    plt.title(title)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def hist2d_psf(e_meas, e_psf, rgp, cov, mask, filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    fig, ax = plt.subplots(1,2)

    tmp = (0.121837929+1.45669783*cov[mask]-25.8743626*cov[mask]**2)*e_psf[mask] / (rgp[mask]**2 - 1)
    ax[0].hist2d(e_meas[mask].real,tmp.real,bins=500,range=((-0.75,0.75),(-0.05,0.05)))
    ax[1].hist2d(e_meas[mask].imag,tmp.imag,bins=500,range=((-0.75,0.75),(-0.05,0.05)))
    ax[0].set_xlabel('e1')
    ax[1].set_xlabel('e2')
    ax[0].set_ylim(-0.01,0.01)
    ax[1].set_ylim(-0.01,0.01)
    plt.tight_layout()
    fig.set_size_inches(8.,4.)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def calc_m(e, g, w=None):
    import numpy
    e1 = e.real
    e2 = e.imag
    g1 = g.real
    g2 = g.imag
    n = len(e1)
    if w is None:
        w = numpy.ones(n)
    if (g1.max() - g1.min()) < 0.001:
        # All g1 are equal.  Just do <e>/g-1.
        m1 = numpy.sum(w*e1) / numpy.sum(w*g1) - 1.
        sigm1 = numpy.sqrt(numpy.sum(w**2 * e1**2) / numpy.sum(w*g1)**2)
    else:
        # Fit for e vs g
        # m = (<xy> - <x><y>) / (<x^2> - <x>^2)
        sw = numpy.sqrt(w)
        m1 = (numpy.mean(w*e1*g1) - numpy.mean(sw*g1) * numpy.mean(sw*e1)) / (numpy.mean(w*g1**2)-numpy.mean(sw*g1)**2) - 1.
        sigm1 = numpy.sqrt(numpy.mean(w*e1**2) / ( n * (numpy.mean(w*g1**2) - numpy.mean(sw*g1)**2)))
        #m1 = (numpy.mean(e1*g1) - g1m * numpy.mean(e1)) / (numpy.mean(g1**2)-g1m**2)-1.
        #sigm1 = numpy.std(e1) / numpy.sqrt(n * (numpy.mean(g1**2) - g1m**2))
    if (g2.max() - g2.min()) < 0.001:
        m2 = numpy.sum(w*e2) / numpy.sum(w*g2) - 1.
        sigm2 = numpy.sqrt(numpy.sum(w**2 * e2**2) / numpy.sum(w*g2)**2)
    else:
        sw = numpy.sqrt(w)
        m2 = (numpy.mean(w*e2*g2) - numpy.mean(sw*g2) * numpy.mean(sw*e2)) / (numpy.mean(w*g2**2)-numpy.mean(sw*g2)**2) - 1.
        sigm2 = numpy.sqrt(numpy.mean(w*e2**2) / ( n * (numpy.mean(w*g2**2) - numpy.mean(sw*g2)**2)))

    return m1, m2, sigm1, sigm2

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

def calc_mc(e, g, w=None, mask=None):
    # Fit e1,e2 separately:
    m1, c1, sigm1, sigc1 = linear_fit(e.real-g.real, g.real, w, mask)
    m2, c2, sigm2, sigc2 = linear_fit(e.imag-g.imag, g.imag, w, mask)

    return m1, m2, c1, c2, sigm1, sigm2, sigc1, sigc2

def calc_nbc(e, g, e_psf, rgp, snr, disc, w, mask, e_true=None):
    import numpy
    n = len(e[mask])
    print 'Start calc_nbc. n = ',n
    sw = numpy.sqrt(w)
    # Tp Tg / (Tp^2 + Tg^2) = Tg/Tp / (1 + Tg^2/Tp^2) = (rgp^2-1) / ( (rgp^2-1)^2 + 1)
    Tp_Tg = 1./(rgp**2-1)
    rgpz = Tp_Tg / (Tp_Tg**2 + 1)
    A = numpy.vstack([ 
                       numpy.ones(n),
                       1. / snr[mask]**2,
                       #rgpz[mask] / snr[mask]**2,
                       #rgpz[mask],
                       #disc[mask],
                       #1. / snr[mask]**2 * disc[mask],
                       #rgp[mask]**2/(rgp[mask]**2-1) / snr[mask]**2 * disc[mask],
                       #rgp[mask]**2/(rgp[mask]**2-1) * disc[mask],
                     ])
    print 'A.shape = ',A.shape
    B = numpy.vstack([ 
                       1. / snr[mask]**2,
                       Tp_Tg[mask] / snr[mask]**2,
                       numpy.ones(n),
                       Tp_Tg[mask],
                     ])
    print 'B.shape = ',B.shape
    M = numpy.vstack([ 
                       A * g[mask] * sw[mask], 
                       B * e_psf[mask] * sw[mask],
                     ]).T
    print 'M.shape = ',M.shape
    if e_true is None:
        b = sw[mask] * (e[mask]-g[mask])
    else:
        b = sw[mask] * (e[mask]-e_true[mask])
    print 'b.shape = ',b.shape
    fit = numpy.linalg.lstsq(M,b)[0]
    print 'fit = ',fit
    fit = fit.real
    print 'fit => ',fit
    #m_corr = fit[0] + fit[1]/snr**2 + fit[2]*rgp + fit[3]*rgp/snr**2
    #c_corr = (fit[4] + fit[5]/snr**2 + fit[6]*rgp + fit[7]/snr**2) * e_psf
    m_corr = fit[0] + fit[1]/snr**2 #+ fit[2]*rgpz
    #m_corr += (fit[4] + (fit[5] + fit[6]*rgp**2/(rgp**2-1))/snr**2 + fit[7]*rgp**2/(rgp**2-1)) * disc
    nm = A.shape[0]
    c_corr = (fit[nm+0] + fit[nm+1]*Tp_Tg)/snr**2 * e_psf
    c_corr += ((fit[nm+2] + fit[nm+3]*Tp_Tg)) * e_psf
    return m_corr, c_corr, fit


def nbc(e, e_true, g_app, e_psf, rgp, snr, disc, w, mask, title='NBC solution', filename=None,
        mask_fit=None):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(2, 2, sharey='row', sharex='col')

    rgp_bins = [ 1.2, 1.25, 1.3, 1.4, 1.6, 2.0, 10. ]
    snr_bins = [ 10, 12.5, 15, 17.5, 20, 25, 30, 50, 80, 200, 2000 ]
    nrgp = len(rgp_bins) - 1
    nsnr = len(snr_bins) - 1
    cm = plt.cm.rainbow

    if mask_fit is None:
        mask_fit = mask
    m_corr, c_corr, fit = calc_nbc(e, g_app, e_psf, rgp, snr, disc, w, mask_fit, e_true=e_true)

    rmean = []
    fits = []
    s_fine = numpy.linspace(snr_bins[0], snr_bins[-1], 10000)

    ax = axes[0,0]
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        print 'r min/max = ',rmin,rmax
        # cf. https://geetduggal.wordpress.com/2011/08/22/grabbing-individual-colors-from-color-maps-in-matplotlib/
        color = cm( (i + 0.5) / nrgp )

        s = numpy.zeros(nsnr)
        m = numpy.zeros(nsnr)
        sm = numpy.zeros(nsnr)
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            if numpy.sum(sample) == 0:
                continue
            s[j] = numpy.mean(snr[sample])
            print 'sj = ',s[j]
            e_comp = e[sample] - e_true[sample]
            m1, _, sm1, _ = linear_fit(e_comp.real, g_app[sample].real, w[sample])
            m2, _, sm2, _ = linear_fit(e_comp.imag, g_app[sample].imag, w[sample])

            m[j] = (m1 + m2) / 2.
            sm[j] = (sm1 + sm2) / numpy.sqrt(2.)  # Not really right, but fine for now.
            print 'm,sigm = ',m[j],sm[j]
        print 's = ',s
        print 'm = ',m
        ax.errorbar(s, m, yerr=sm, color=color, fmt='o', label='%.2f < rgpp_rp < %.2f'%(rmin,rmax))

        mask2 = mask & (rgp >= rmin) & (rgp < rmax)
        #meanrgp = numpy.mean(rgp[mask2])
        #meandisc = numpy.mean(disc[mask2])
        #meanrgpsq = numpy.mean(rgp[mask2]**2)
        # Tobs / Tg = (Tobs/Tp) / ((Tobs-Tp)/Tp) = rgp^2 / (rgp^2-1)
        #meanrgpx = numpy.mean(rgp[mask2]**2/(rgp[mask2]**2-1))
        # Tg/Tp = 1 / ((Tobs-Tp)/Tp) = 1 / (rgp^2-1)
        #meanrgpy = numpy.mean(1/(rgp[mask2]**2-1))
        #meanrgpysq = numpy.mean(1/(rgp[mask2]**2-1)**2)
        # Tp Tg / (Tp^2 + Tg^2) = Tp/Tg / (1 + Tp^2/Tg^2) = (rgp^2-1) / ( (rgp^2-1)^2 + 1)
        #meanrgpz = numpy.mean( (rgp[mask2]**2-1) / ( (rgp[mask2]**2-1)**2 + 1) )
        fitline = fit[0] + fit[1]/s_fine**2 #+ fit[2]*meanrgpz
        #fitline += (fit[4] + (fit[5] + fit[6]*meanrgpx)/s_fine**2 + fit[7]*meanrgpx)*meandisc
        print 'fitline = ',fitline
        ax.plot(s_fine, fitline, color=color)
    nm = 2 # How many fit coefficients are in the m part of the fit.
    #ax.set_title(r'fit parameters $\{1, TpTg/(Tp^2+Tg^2)\} \times \{1, \nu^{-2}\}$')
    ax.set_title(r'fit parameters $\{1, \nu^{-2}\}$')

    ax.legend(loc='upper right', fontsize=10)
    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('m')
    ax.set_ylim(-0.3,0.3)
    ax.set_xlim(1.e1, 1.e3)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])

    ax = axes[0,1]

    e_corr = e - c_corr
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        print 'r min/max = ',rmin,rmax
        color = cm( (i + 0.5) / nrgp )

        s = numpy.zeros(nsnr)
        m = numpy.zeros(nsnr)
        sm = numpy.zeros(nsnr)
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            if numpy.sum(sample) == 0:
                continue
            means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
            e_corr = (e[sample] - c_corr[sample]) / means
            s[j] = numpy.mean(snr[sample])
            print 'sj = ',s[j]
            print 'mean c_corr = ',numpy.mean(c_corr[sample])
            print 'mean m_corr = ',means-1
            e_comp = e_corr - e_true[sample]
            m1, _, sm1, _ = linear_fit(e_comp.real, g_app[sample].real, w[sample])
            m2, _, sm2, _ = linear_fit(e_comp.imag, g_app[sample].imag, w[sample])

            m[j] = (m1 + m2) / 2.
            sm[j] = (sm1 + sm2) / numpy.sqrt(2.)  # Not really right, but fine for now.
        ax.errorbar(s, m, yerr=sm, color=color, fmt='o', label='%.2f < rgpp_rp < %.2f'%(rmin,rmax))

    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('m after correction')
    ax.set_xlim(1.e1, 1.e3)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])
    ax.set_title(r'Corrected Shapes')

    ax = axes[1,0]
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        print 'r min/max = ',rmin,rmax
        # cf. https://geetduggal.wordpress.com/2011/08/22/grabbing-individual-colors-from-color-maps-in-matplotlib/
        color = cm( (i + 0.5) / nrgp )

        s = numpy.zeros(nsnr)
        a = numpy.zeros(nsnr)
        sa = numpy.zeros(nsnr)
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            if numpy.sum(sample) == 0:
                continue
            s[j] = numpy.mean(snr[sample])
            print 'sj = ',s[j]
            e_comp = e[sample] - e_true[sample]
            a1, _, sa1, _ = linear_fit(e_comp.real, e_psf[sample].real, w[sample])
            a2, _, sa2, _ = linear_fit(e_comp.imag, e_psf[sample].imag, w[sample])

            a[j] = (a1 + a2) / 2.
            sa[j] = (sa1 + sa2) / numpy.sqrt(2.)  # Not really right, but fine for now.
            print 'a,siga = ',a[j],sa[j]
        print 's = ',s
        print 'alpha = ',a
        ax.errorbar(s, a, yerr=sa, color=color, fmt='o', label='%.2f < rgpp_rp < %.2f'%(rmin,rmax))

        mask2 = mask & (rgp >= rmin) & (rgp < rmax)
        #meanrgp = numpy.mean(rgp[mask2])
        meanrgpx = numpy.mean(1./(rgp[mask2]**2-1))
        #meanrgpsq = numpy.mean(rgp[mask2]**2)
        fitline = (fit[nm+0] + fit[nm+1]*meanrgpx)/s_fine**2
        fitline += (fit[nm+2] + fit[nm+3]*meanrgpx)
        print 'fitline = ',fitline
        ax.plot(s_fine, fitline, color=color)
    ax.set_title(r'fit parameters $\{1, 1/(rgp^2-1)\} \times \{1, \nu^{-2}\}$')

    ax.legend(loc='upper right', fontsize=10)
    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('alpha')
    ax.set_ylim(-0.2,0.5)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])

    ax = axes[1,1]

    e_corr = e - c_corr
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        print 'r min/max = ',rmin,rmax
        color = cm( (i + 0.5) / nrgp )

        s = numpy.zeros(nsnr)
        a = numpy.zeros(nsnr)
        sa = numpy.zeros(nsnr)
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            if numpy.sum(sample) == 0:
                continue
            means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
            e_corr = (e[sample] - c_corr[sample]) / means
            s[j] = numpy.mean(snr[sample])
            print 'sj = ',s[j]
            print 'mean c_corr = ',numpy.mean(c_corr[sample])
            print 'mean m_corr = ',means-1
            e_comp = e_corr - e_true[sample]
            a1, _, sa1, _ = linear_fit(e_comp.real, e_psf[sample].real, w[sample])
            a2, _, sa2, _ = linear_fit(e_comp.imag, e_psf[sample].imag, w[sample])

            a[j] = (a1 + a2) / 2.
            sa[j] = (sa1 + sa2) / numpy.sqrt(2.)  # Not really right, but fine for now.
        ax.errorbar(s, a, yerr=sa, color=color, fmt='o', label='%.2f < rgpp_rp < %.2f'%(rmin,rmax))

    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('alpha after correction')
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])
    ax.set_title(r'Corrected Shapes')

    fig.set_size_inches(12.,12.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

    return m_corr, c_corr

def plot_mean_x_vs_y(ax, x, y, w=None, mask=None, color='red', fmt='o', label=None, line=False, nx=10):
    import numpy

    if w is None:
        w = numpy.ones(len(x))
    if mask is None:
        mask = numpy.ones(len(x),dtype=bool)

    sortx = numpy.sort(x[mask])
    print len(sortx)
    deciles = numpy.array([ sortx[i*len(sortx)/nx] for i in range(nx) ] + [ sortx[-1] ])
    print deciles

    mx = numpy.zeros(nx)
    my = numpy.zeros(nx)
    sy = numpy.zeros(nx)
    use = numpy.ones(nx, dtype=bool)

    for j in range(nx):
        sample = numpy.where(mask & (x >= deciles[j]) & (x < deciles[j+1]))
        if numpy.sum(sample) == 0:
            use[j] = False
            continue
        mx[j] = numpy.sum(w[sample] * x[sample]) / numpy.sum(w[sample])
        my[j] = numpy.sum(w[sample] * y[sample]) / numpy.sum(w[sample])
        sy[j] = numpy.sum(w[sample]**2 * (y[sample]-my[j])**2) / numpy.sum(w[sample])**2

    sy = numpy.sqrt(sy)
    print 'mx = ',mx
    if line:
        ax.plot(mx[use], my[use], color=color, label=label)
    else:
        ax.errorbar(mx[use], my[use], yerr=sy[use], color=color, fmt=fmt, label=label)
    if min(my[use]) < 0 and max(my[use]) > 0:
        ax.plot( [mx[0],mx[-1]], [0.,0.], color='k')


def mean_x_vs_y(x, y, w, mask, xlabel, ylabel, filename, title=None):
    import matplotlib.pyplot as plt
    import matplotlib
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 8})
    fig, ax = plt.subplots(1, 1)

    plot_mean_x_vs_y(ax, x, y, w, mask)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    plt.savefig(filename)


def alpha_vs_x(e, e_true, g_app, e_psf, x, w, mask, label, filename, title=None, m_corr=None, c_corr=None, log=False):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 8})
    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row')

    ax = axes[0]
    nx = 10
    sortx = numpy.sort(x[mask])
    print len(sortx)
    deciles = numpy.array([ sortx[i*len(sortx)/nx] for i in range(nx) ] + [ sortx[-1] ])
    print deciles

    mx = numpy.zeros(nx)
    a1 = numpy.zeros(nx)
    sa1 = numpy.zeros(nx)
    a2 = numpy.zeros(nx)
    sa2 = numpy.zeros(nx)
    a1_t = numpy.zeros(nx)
    sa1_t = numpy.zeros(nx)
    a2_t = numpy.zeros(nx)
    sa2_t = numpy.zeros(nx)
    a1_c = numpy.zeros(nx)
    sa1_c = numpy.zeros(nx)
    a2_c = numpy.zeros(nx)
    sa2_c = numpy.zeros(nx)
    use = numpy.ones(nx, dtype=bool)

    if m_corr is not None or c_corr is not None:
        if m_corr is None:
            m_corr = numpy.zeros(len(e))
        if c_corr is None:
            c_corr = numpy.zeros(len(e))
        corr = True
    else:
        corr = False

    for j in range(nx):
        sample = numpy.where(mask & (x >= deciles[j]) & (x < deciles[j+1]))
        if numpy.sum(sample) == 0:
            use[j] = False
            continue
        mx[j] = numpy.mean(x[sample])
        a1[j], _, sa1[j], _ = linear_fit(e[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2[j], _, sa2[j], _ = linear_fit(e[sample].imag-g_app[sample].imag, e_psf[sample].imag, w[sample])
        a1_t[j], _, sa1_t[j], _ = linear_fit(e_true[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2_t[j], _, sa2_t[j], _ = linear_fit(e_true[sample].imag-g_app[sample].imag, e_psf[sample].imag, w[sample])
        if corr:
            means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
            e_corr = (e[sample] - c_corr[sample]) / means
            a1_c[j], _, sa1_c[j], _ = linear_fit(e_corr.real-g_app[sample].real, e_psf[sample].real, w[sample])
            a2_c[j], _, sa2_c[j], _ = linear_fit(e_corr.imag-g_app[sample].imag, e_psf[sample].imag, w[sample])
            
    print 'mx = ',mx
    print 'alpha1_t = ',a1_t
    print 'alpha2_t = ',a2_t
    if corr:
        print 'alpha1_c = ',a1_c
        print 'alpha2_c = ',a2_c

    if title is None or 'ngmix' not in title:
        ax.errorbar(mx[use], a1[use], yerr=sa1[use], color='red', fmt='o', label=r'$\alpha_1$ measured shapes')
        ax.errorbar(mx[use], a2[use], yerr=sa2[use], color='magenta', fmt='o', label=r'$\alpha_2$ measured shapes')
    ax.errorbar(mx[use], a1_t[use], yerr=sa1_t[use], color='green', fmt='o', label=r'$\alpha_1$ true shapes')
    ax.errorbar(mx[use], a2_t[use], yerr=sa2_t[use], color='LimeGreen', fmt='o', label=r'$\alpha_2$ true shapes')
    if corr:
        ax.errorbar(mx[use], a1_c[use], yerr=sa1_c[use], color='blue', fmt='o', label=r'$\alpha_1$ corrected shapes')
        ax.errorbar(mx[use], a2_c[use], yerr=sa2_c[use], color='cyan', fmt='o', label=r'$\alpha_2$ corrected shapes')

    ax.plot( [mx[0],mx[-1]], [0.,0.], color='k')
    if log:
        ax.set_xscale('log')
    ax.set_ylabel('alpha')
    #ax.set_xlim(sortx[0], sortx[-1])
    if title is not None:
        ax.set_title(title)

    ax = axes[1]
    for j in range(nx):
        if not use[j]:
            continue
        sample = numpy.where(mask & (x >= deciles[j]) & (x < deciles[j+1]))
        mx[j] = numpy.mean(x[sample])
        a1[j], _, sa1[j], _ = linear_fit(e[sample].real-g_app[sample].real, g_app[sample].real, w[sample])
        a2[j], _, sa2[j], _ = linear_fit(e[sample].imag-g_app[sample].imag, g_app[sample].imag, w[sample])
        a1_t[j], _, sa1_t[j], _ = linear_fit(e_true[sample].real-g_app[sample].real, g_app[sample].real, w[sample])
        a2_t[j], _, sa2_t[j], _ = linear_fit(e_true[sample].imag-g_app[sample].imag, g_app[sample].imag, w[sample])
        if corr:
            means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
            e_corr = (e[sample] - c_corr[sample]) / means
            a1_c[j], _, sa1_c[j], _ = linear_fit(e_corr.real-g_app[sample].real, g_app[sample].real, w[sample])
            a2_c[j], _, sa2_c[j], _ = linear_fit(e_corr.imag-g_app[sample].imag, g_app[sample].imag, w[sample])
    print 'mx = ',mx
    print 'm1_t = ',a1_t
    print 'm2_t = ',a2_t
    if corr:
        print 'm1_c = ',a1_c
        print 'm2_c = ',a2_c

    if title is None or 'ngmix' not in title:
        ax.errorbar(mx[use], a1[use], yerr=sa1[use], color='red', fmt='o', label=r'$m_1$ true shapes')
        ax.errorbar(mx[use], a2[use], yerr=sa2[use], color='magenta', fmt='o', label=r'$m_2$ true shapes')
    ax.errorbar(mx[use], a1_t[use], yerr=sa1_t[use], color='green', fmt='o', label=r'$m_1$ true shapes')
    ax.errorbar(mx[use], a2_t[use], yerr=sa2_t[use], color='LimeGreen', fmt='o', label=r'$m_2$ true shapes')
    if corr:
        ax.errorbar(mx[use], a1_c[use], yerr=sa1_c[use], color='blue', fmt='o', label=r'$m_1$ corrected shapes')
        ax.errorbar(mx[use], a2_c[use], yerr=sa2_c[use], color='cyan', fmt='o', label=r'$m_2$ corrected shapes')

    ax.plot( [mx[0],mx[-1]], [0.,0.], color='k')
    ax.set_xlabel(label)
    if log:
        ax.set_xscale('log')
    ax.set_ylabel('m')
    #ax.set_xlim(sortx[0], sortx[-1])

    fig.set_size_inches(5.5,8.0)
    plt.tight_layout()
    plt.savefig(filename)


def alpha_vs_xy(e, e_true, g_app, e_psf, x, y, w, mask, xlabel, ylabel, filename):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 8})
    fig, ax = plt.subplots(2, 3, sharex='col', sharey='row')

    nbin = 30
    sortx = numpy.sort(x[mask])
    print len(sortx)
    xbins = numpy.array([ sortx[i*len(sortx)/nbin] for i in range(nbin) ] + [ sortx[-1] ])
    print xbins

    sorty = numpy.sort(y[mask])
    ybins = numpy.array([ sorty[i*len(sorty)/nbin] for i in range(nbin) ] + [ sorty[-1] ])
    print ybins

    xx,yy = numpy.meshgrid(xbins,ybins)

    alpha_m = numpy.zeros((nbin,nbin))
    alpha_t = numpy.zeros((nbin,nbin))
    m_m = numpy.zeros((nbin,nbin))
    m_t = numpy.zeros((nbin,nbin))
    for i in range(nbin):
      for j in range(nbin):
        sample = numpy.where(mask & (x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1]))
        if numpy.sum(sample) == 0:
            continue
        alpha_m[j,i] = linear_fit(e[sample]-g_app[sample], e_psf[sample], w[sample])[0].real
        alpha_t[j,i] = linear_fit(e_true[sample]-g_app[sample], e_psf[sample], w[sample])[0].real
        m_m[j,i] = linear_fit(e[sample]-g_app[sample], g_app[sample], w[sample])[0].real
        m_t[j,i] = linear_fit(e_true[sample]-g_app[sample], g_app[sample], w[sample])[0].real
        xx[j,i] = numpy.mean(x[sample])
        yy[j,i] = numpy.mean(y[sample])

    amax = numpy.max([alpha_m, alpha_t, alpha_m-alpha_t])
    amin = numpy.min([alpha_m, alpha_t, alpha_m-alpha_t])
    alim = numpy.max([amax, -amin])
    mmax = numpy.max([m_m, m_t, m_m-m_t])
    mmin = numpy.min([m_m, m_t, m_m-m_t])
    mlim = numpy.max([mmax, -mmin])

    ax[0,0].pcolorfast(xx, yy, alpha_m, cmap='RdBu', vmin=-alim, vmax=alim)
    ax[0,1].pcolorfast(xx, yy, alpha_t, cmap='RdBu', vmin=-alim, vmax=alim)
    im = ax[0,2].pcolorfast(xx, yy, alpha_m-alpha_t, cmap='RdBu', vmin=-alim, vmax=alim)
    plt.colorbar(im, ax=ax[0,0])
    plt.colorbar(im, ax=ax[0,1])
    plt.colorbar(im, ax=ax[0,2])

    ax[1,0].pcolorfast(xx, yy, m_m, cmap='RdBu', vmin=-mlim, vmax=mlim)
    ax[1,1].pcolorfast(xx, yy, m_t, cmap='RdBu', vmin=-mlim, vmax=mlim)
    im = ax[1,2].pcolorfast(xx, yy, m_m-m_t, cmap='RdBu', vmin=-mlim, vmax=mlim)
    plt.colorbar(im, ax=ax[1,0])
    plt.colorbar(im, ax=ax[1,1])
    plt.colorbar(im, ax=ax[1,2])

    ax[0,0].set_ylabel(ylabel)
    ax[1,0].set_ylabel(ylabel)
    ax[1,0].set_xlabel(xlabel)
    ax[1,1].set_xlabel(xlabel)
    ax[1,2].set_xlabel(xlabel)

    ax[0,0].set_title('alpha for $e_{meas}$')
    ax[0,1].set_title('alpha for $e_{true}$')
    ax[0,2].set_title('difference')
    ax[1,0].set_title('m for $e_{meas}$')
    ax[1,1].set_title('m for $e_{true}$')
    ax[1,2].set_title('difference')

    fig.set_size_inches(13.0,8.0)
    #plt.tight_layout()
    plt.savefig(filename)


def plot_nbc(ax, rgp, snr, e, g_app, w, mask, rgp_bins, snr_bins, e_true=None):
    import matplotlib.pyplot as plt
    import numpy
    nrgp = len(rgp_bins) - 1
    nsnr = len(snr_bins) - 1
    cm = plt.cm.rainbow

    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        print 'r min/max = ',rmin,rmax
        # cf. https://geetduggal.wordpress.com/2011/08/22/grabbing-individual-colors-from-color-maps-in-matplotlib/
        color = cm( (i + 0.5) / nrgp )
        s = numpy.zeros(nsnr)
        m = numpy.zeros(nsnr)
        sm = numpy.zeros(nsnr)
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            if numpy.sum(sample) == 0:
                continue
            s[j] = numpy.mean(snr[sample])
            #print 'snr = ',s[j]
            #print 'n = ',numpy.sum(sample)
            m1, m2, sigm1, sigm2 = calc_m(e[sample], g_app[sample], w[sample])
            if e_true is not None:
                m1t, m2t, sigm1t, sigm2t = calc_m(e_true[sample], g_app[sample], w[sample])
                m1 -= m1t
                m2 -= m2t

            #print 'm1 = %f +- %f'%(m1,sigm1)
            #print 'm2 = %f +- %f'%(m2,sigm2)
            m[j] = (m1 + m2) / 2.
            sm[j] = (sigm1 + sigm2) / numpy.sqrt(2.)  # Not really right, but fine for now.
        #print 's = ',s
        #print 'm = ',m
        #print 'sm = ',sm
        ax.errorbar(s, m, yerr=sm, color=color, fmt='o', label='%.2f < rgpp_rp < %.2f'%(rmin,rmax))



def old_nbc(e, e_true, g_app, rgp, snr, w, mask, title='nbc', filename=None, mask_fit=None):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(2, 2)

    #rgp_bins = [ 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6 ]
    rgp_bins = [ 1.2, 1.25, 1.3, 1.4, 1.6, 2.0, 3.0 ]
    snr_bins = [ 10, 12.5, 15, 17.5, 20, 25, 30, 50, 80, 200, 2000 ]
    #snr_bins = [ 15, 17.5, 20, 25, 30, 50, 80, 200, 2000 ]
    nrgp = len(rgp_bins) - 1
    nsnr = len(snr_bins) - 1
    cm = plt.cm.rainbow

    rmean = []
    fits = []
    s_fine = numpy.linspace(snr_bins[0], snr_bins[-1], 10000)

    if mask_fit is None:
        mask_fit = mask
    ax = axes[0,0]
    plot_nbc(ax, rgp, snr, e, g_app, w, mask, rgp_bins, snr_bins, e_true=e_true)
    sw = numpy.sqrt(w)
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        color = cm( (i + 0.5) / nrgp )
        mask2 = mask & (rgp >= rmin) & (rgp < rmax)
        if numpy.sum(mask2) == 0:
            continue
        # fit points to m = a/nu^2 + c
        n = numpy.sum(mask2)
        a = (numpy.vstack([ g_app[mask2]/snr[mask2]**2 ]) * sw[mask2]).t
        if e_true is None:
            b = sw[mask2] * (e[mask2]-g_app[mask2])
        else:
            b = sw[mask2] * (e[mask2]-e_true[mask2])
        fit = numpy.linalg.lstsq(a,b)[0]
        print 'fit = ',fit
        a, = fit
        fits.append(a.real)
        rmean.append(numpy.mean(rgp[mask2]))
        ax.plot(s_fine, a/s_fine**2, color=color)

    #ax.set_yscale('log')
    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel('$s/n$')
    ax.set_ylabel('m')
    ax.set_ylim(-0.5,0.5)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])
    ax.set_title(title)
    ax.legend(loc='lower right')

    # look at fit vs rgp:
    ax = axes[0,1]
    ax.scatter(rmean, fits)
    ax.set_xlabel(r'<rgpp_rp>')
    ax.set_ylabel(r'fit coeff of $m$ vs $1/nu^2$')
    ax.set_title(r'dependence of nbc on rgp')

    # this is very linear.  so repeat this as a single fit using the functional form:
    # m = a/nu^2 + b rgp/nu^2
    # e-g = ( a/nu^2 + b rgp/nu^2 + c rgp^2/nu^2 + d ) g
    ax = axes[1,0]
    m_corr, fit = calc_nbc(e, g_app, rgp, snr, w, mask_fit, e_true=e_true)
    #junk = [1./10.**2, 1./10**4, 1./10**3, 1./10**1]
    #print 'junk = ',junk
    #print 'fit.dot(junk) = ',fit.dot(junk)
    #print 'fit.real.dot(junk) = ',fit.real.dot(junk)
    #junk2 = numpy.array([1./s_fine**2, 1./s_fine*4, 1./s_fine*3, 1./s_fine*1])
    #print 'junk2 = ',junk2
    #print 'fit.real.dot(junk2) = ',fit.real.dot(junk)

    plot_nbc(ax, rgp, snr, e, g_app, w, mask, rgp_bins, snr_bins, e_true=e_true)
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        color = cm( (i + 0.5) / nrgp )
        meanrgp = numpy.mean(rgp[mask & (rgp >= rmin) & (rgp < rmax)])
        meanrgpsq = numpy.mean(rgp[mask & (rgp >= rmin) & (rgp < rmax)]**2)
        #fitline = fit[0].real/s_fine**2
        #fitline += fit[1].real/s_fine**4
        #fitline += fit[2].real/s_fine**3
        #fitline += fit[3].real/s_fine**1
        #fitline += fit[1].real
        fitline = (fit[0].real + fit[1].real*meanrgp)/s_fine**2
        #fitline += (fit[2].real + fit[3].real*meanrgp)/s_fine**4
        #fitline += (fit[2].real + fit[3].real*meanrgp)/s_fine**3
        #fitline += (fit[4].real + fit[5].real*meanrgp)/s_fine**1
        #fitline += (fit[6].real + fit[7].real*meanrgp)
        #fitline += fit[2].real + fit[3].real*meanrgp
        #fitline = (fit[0].real + fit[1].real*meanrgp + fit[2].real*meanrgpsq)/s_fine**2
        #fitline += (fit[3].real + fit[4].real*meanrgp + fit[5].real*meanrgpsq)/s_fine**4
        #fitline += (fit[6].real + fit[7].real*meanrgp + fit[8].real*meanrgpsq)/s_fine**3
        #fitline += (fit[9].real + fit[10].real*meanrgp + fit[11].real*meanrgpsq)/s_fine**1
        print 'fitline = ',fitline
        ax.plot(s_fine, fitline, color=color)

    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel('$s/n$')
    ax.set_ylabel('m')
    ax.set_ylim(-0.5,0.5)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])
    ax.set_title(r'fit parameters $\{1,rgp\}\times\{nu^{-0},nu^{-1},nu^{-2},nu^{-3}\}$')
    #ax.set_title(r'fit parameters $\{1,rgp\}\times\{1,nu^{-2}\}$')
    #ax.set_title(r'fit parameters $\{1\}\times\{1,nu^{-2}\}$')

    # finally, repeat with the correction applied directly to the shapes.
    ax = axes[1,1]

    e_corr = e.copy()
    for i in range(nrgp):
        rmin = rgp_bins[i]
        rmax = rgp_bins[i+1]
        for j in range(nsnr):
            smin = snr_bins[j]
            smax = snr_bins[j+1]
            sample = mask & (rgp >= rmin) & (rgp < rmax) & (snr >= smin) & (snr < smax)
            means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
            e_corr[sample] /= means
    plot_nbc(ax, rgp, snr, e_corr, g_app, w, mask, rgp_bins, snr_bins, e_true=e_true)

    ax.plot( [snr_bins[0],snr_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel('$s/n$')
    ax.set_ylabel('m after correction')
    ax.set_ylim(-0.5,0.5)
    ax.set_xscale('log')
    ax.set_xlim(snr_bins[0], snr_bins[-1])
    ax.set_title(r'corrected shapes')

    fig.set_size_inches(12.,12.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

    return m_corr



def calc_nbc_cov(e, g, e_psf, rgp, cov, w, mask):
    import numpy
    n = len(e[mask])
    print 'Start calc_nbc_cov. n = ',n
    sw = numpy.sqrt(w)
    A = numpy.vstack([ 
                       numpy.ones(n),
                       cov[mask],
                       cov[mask]**2,
                     ])
    print 'A.shape = ',A.shape
    M = numpy.vstack([ 
                       A * g[mask] * sw[mask], 
                       A * e_psf[mask] * sw[mask],
                       #A * e_psf[mask] / (rgp[mask]**2-1) * sw[mask],
                       #A * e_psf[mask].conjugate() / (rgp[mask]**2-1) * sw[mask]
                       #A * e_psf[mask].real / (rgp[mask]**2-1) * sw[mask],
                       #A * 1j * e_psf[mask].imag / (rgp[mask]**2-1) * sw[mask]
                     ]).T
    print 'M.shape = ',M.shape
    b = sw[mask] * (e[mask]-g[mask])
    print 'b.shape = ',b.shape
    fit = numpy.linalg.lstsq(M,b)[0]
    print 'fit = ',fit
    fit = fit.real
    print 'fit => ',fit
    m_corr = fit[0] + fit[1]*cov + fit[2]*cov**2
    #c_corr = (fit[3] + fit[4]/(rgp**2-1.)) * e_psf
    c_corr = (fit[3] + fit[4]*cov + fit[5]*cov**2) * e_psf
    #c_corr = (fit[3] + fit[4]*cov + fit[5]*cov**2)/(rgp**2-1.) * e_psf
    #c_corr += (fit[6] + fit[7]*cov + fit[8]*cov**2)/(rgp**2-1.) * e_psf.conjugate()
    #c_corr = (fit[3] + fit[4]*cov + fit[5]*cov**2)/(rgp**2-1.) * (1+0j) * e_psf.real
    #c_corr += (fit[6] + fit[7]*cov + fit[8]*cov**2)/(rgp**2-1.) * 1j * e_psf.imag
    return m_corr, c_corr, fit


def nbc_c(e, e_true, g_app, e_psf, rgp, cov, w, mask, title=None, filename=None, mask_fit=None):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 8})
    #fig, axes = plt.subplots(1, 2, sharex='col', sharey='row')
    fig, axes = plt.subplots(2, 1)

    rgp_bins = numpy.linspace(1.2, 2, num=10)
    nrgp = len(rgp_bins) - 1

    if mask_fit is None:
        mask_fit = mask

    m_corr, c_corr, fit = calc_nbc_cov(e-e_true, g_app, e_psf, rgp, cov, w, mask_fit)
    e_corr = e.copy()
    e_corr -= c_corr

    ax = axes[0]
    r = numpy.zeros(nrgp)
    a1 = numpy.zeros(nrgp)
    sa1 = numpy.zeros(nrgp)
    a2 = numpy.zeros(nrgp)
    sa2 = numpy.zeros(nrgp)
    a1_t = numpy.zeros(nrgp)
    sa1_t = numpy.zeros(nrgp)
    a2_t = numpy.zeros(nrgp)
    sa2_t = numpy.zeros(nrgp)
    a1_c = numpy.zeros(nrgp)
    sa1_c = numpy.zeros(nrgp)
    a2_c = numpy.zeros(nrgp)
    sa2_c = numpy.zeros(nrgp)
    for j in range(nrgp):
        rmin = rgp_bins[j]
        rmax = rgp_bins[j+1]
        sample = mask & (rgp >= rmin) & (rgp < rmax)
        if numpy.sum(sample) == 0:
            continue
        r[j] = numpy.mean(rgp[sample])
        a1[j], _, sa1[j], _ = linear_fit(e[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2[j], _, sa2[j], _ = linear_fit(e[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])
        a1_t[j], _, sa1_t[j], _ = linear_fit(e_true[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2_t[j], _, sa2_t[j], _ = linear_fit(e_true[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])
        means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
        e_corr[sample] /= means
        a1_c[j], _, sa1_c[j], _ = linear_fit(e_corr[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2_c[j], _, sa2_c[j], _ = linear_fit(e_corr[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])
    print 'alpha1 = ',a1_c
    print 'alpha2 = ',a2_c

    ax.errorbar(r, a1, yerr=sa1, color='red', fmt='o', label=r'Uncorrected $\alpha_1$')
    ax.errorbar(r, a2, yerr=sa2, color='magenta', fmt='o', label=r'Uncorrected $\alpha_2$')
    ax.errorbar(r, a1_t, yerr=sa1_t, color='green', fmt='o', label=r'$\alpha_1$ for true shapes')
    ax.errorbar(r, a2_t, yerr=sa2_t, color='LimeGreen', fmt='o', label=r'$\alpha_2$ for true shapes')
    ax.errorbar(r, a1_c, yerr=sa1_c, color='blue', fmt='o', label=r'Corrected $\alpha_1$')
    ax.errorbar(r, a2_c, yerr=sa2_c, color='cyan', fmt='o', label=r'Corrected $\alpha_2$')

    #fitline = fit[0] * c_fine + fit[1] + fit[2]*c_fine**2
    #ax.plot(c_fine, fitline, color='r')
    ax.legend(loc='lower left', fontsize=6)

    ax.plot( [rgp_bins[0],rgp_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel('rgpp_rp')
    ax.set_ylabel(r'$\alpha$')
    ax.set_xlim(rgp_bins[0], rgp_bins[-1])
    if title is None:
        ax.set_title(r'Fit parameters $1$, $(S/N)_r^{-2}$, $(S/N)_r^{-4}$')
    else:
        ax.set_title(title)

    cov_bins = numpy.linspace(0, 0.0045, num=10)
    #cov_bins = numpy.linspace(0, 0.2, num=10)
    ncov = len(cov_bins) - 1

    ax = axes[1]
    c = numpy.zeros(ncov)
    a1 = numpy.zeros(ncov)
    sa1 = numpy.zeros(ncov)
    a2 = numpy.zeros(ncov)
    sa2 = numpy.zeros(ncov)
    a1_c = numpy.zeros(ncov)
    sa1_c = numpy.zeros(ncov)
    a2_c = numpy.zeros(ncov)
    sa2_c = numpy.zeros(ncov)
    a1_t = numpy.zeros(ncov)
    sa1_t = numpy.zeros(ncov)
    a2_t = numpy.zeros(ncov)
    sa2_t = numpy.zeros(ncov)
    for j in range(ncov):
        cmin = cov_bins[j]
        cmax = cov_bins[j+1]
        sample = mask & (cov >= cmin) & (cov < cmax)
        if numpy.sum(sample) == 0:
            continue
        c[j] = numpy.mean(cov[sample])
        a1[j], _, sa1[j], _ = linear_fit(e[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2[j], _, sa2[j], _ = linear_fit(e[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])

        a1_t[j], _, sa1_t[j], _ = linear_fit(e_true[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2_t[j], _, sa2_t[j], _ = linear_fit(e_true[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])

        means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
        e_corr[sample] /= means
        a1_c[j], _, sa1_c[j], _ = linear_fit(e_corr[sample].real-g_app[sample].real, e_psf[sample].real, w[sample])
        a2_c[j], _, sa2_c[j], _ = linear_fit(e_corr[sample].imag-g_app[sample].real, e_psf[sample].imag, w[sample])
    print 'alpha1 = ',a1_c
    print 'alpha2 = ',a2_c

    ax.errorbar(c, a1, yerr=sa1, color='red', fmt='o', label=r'Uncorrected $\alpha_1$')
    ax.errorbar(c, a2, yerr=sa2, color='magenta', fmt='o', label=r'Uncorrected $\alpha_2$')
    ax.errorbar(c, a1_t, yerr=sa1_t, color='green', fmt='o', label=r'$\alpha_1$ for true shapes')
    ax.errorbar(c, a2_t, yerr=sa2_t, color='LimeGreen', fmt='o', label=r'$\alpha_2$ for true shapes')
    ax.errorbar(c, a1_c, yerr=sa1_c, color='b', fmt='o', label=r'Corrected $\alpha_1$')
    ax.errorbar(c, a2_c, yerr=sa2_c, color='cyan', fmt='o', label=r'Corrected $\alpha_2$')

    #fitline = fit[0] * c_fine + fit[1] + fit[2]*c_fine**2
    #ax.plot(c_fine, fitline, color='r')
    ax.legend(loc='upper left', fontsize=6)

    ax.plot( [cov_bins[0],cov_bins[-1]], [0.,0.], color='k')
    #ax.set_xlabel('$C^\prime = Tr(Cov) / \sqrt{1-|e|^2}$')
    ax.set_xlabel(r'$(S/N)_r^{-2}$')
    ax.set_ylabel(r'$\alpha$')
    ax.set_xlim(cov_bins[0], cov_bins[-1])

    fig.set_size_inches(5.5,8.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

    return m_corr, c_corr



def nbc_cov(e, e_true, g_app, e_psf, rgp, cov, w, mask, filename=None, mask_fit=None, title=None):
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy
    plt.style.use('supermongo')
    matplotlib.rcParams.update({'font.size': 12})
    fig, ax = plt.subplots(1, 1)

    cov_bins = numpy.linspace(0, 0.0045, num=10)
    #cov_bins = numpy.linspace(0, 0.2, num=10)
    ncov = len(cov_bins) - 1

    rmean = []
    fits = []
    c_fine = numpy.linspace(cov_bins[0], cov_bins[-1], num=100)

    if mask_fit is None:
        mask_fit = mask
    m_corr, c_corr, fit = calc_nbc_cov(e-e_true, g_app, e_psf, rgp, cov, w, mask_fit)
    e_corr = e.copy()
    e_corr -= c_corr

    c = numpy.zeros(ncov)
    m = numpy.zeros(ncov)
    sm = numpy.zeros(ncov)
    m_c = numpy.zeros(ncov)
    sm_c = numpy.zeros(ncov)
    m_t = numpy.zeros(ncov)
    sm_t = numpy.zeros(ncov)
    for j in range(ncov):
        smin = cov_bins[j]
        smax = cov_bins[j+1]
        sample = mask & (cov >= smin) & (cov < smax)
        if numpy.sum(sample) == 0:
            continue
        c[j] = numpy.mean(cov[sample])
        m1, m2, sigm1, sigm2 = calc_m(e[sample], g_app[sample], w[sample])
        m[j] = (m1 + m2) / 2.
        sm[j] = (sigm1 + sigm2) / numpy.sqrt(2.)  # Not really right, but fine for now.

        m1, m2, sigm1, sigm2 = calc_m(e_true[sample], g_app[sample], w[sample])
        m_t[j] = (m1 + m2) / 2.
        sm_t[j] = (sigm1 + sigm2) / numpy.sqrt(2.)  # Not really right, but fine for now.

        means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
        e_corr[sample] /= means
        m1, m2, sigm1, sigm2 = calc_m(e_corr[sample], g_app[sample], w[sample])
        m_c[j] = (m1 + m2) / 2.
        sm_c[j] = (sigm1 + sigm2) / numpy.sqrt(2.)  # Not really right, but fine for now.

    ax.errorbar(c, m, yerr=sm, color='r', fmt='o', label='Uncorrected')
    ax.errorbar(c, m_c, yerr=sm_c, color='b', fmt='o', label='Corrected')
    ax.errorbar(c, m_t, yerr=sm_t, color='g', fmt='o', label='True shapes')

    fitline = fit[0] + fit[1] * c_fine + fit[2]*c_fine**2
    ax.plot(c_fine, fitline, color='r')
    ax.legend(loc='lower left')

    ax.plot( [cov_bins[0],cov_bins[-1]], [0.,0.], color='k')
    ax.set_xlabel('$(S/N)_r^{-2}$')
    #ax.set_xlabel('$\sqrt{det(Cov(e1,e2))}$')
    ax.set_ylabel('m')
    ax.set_xlim(cov_bins[0], cov_bins[-1])
    if title is None:
        ax.set_title(r'Fit parameters $1$, $(S/N)_r^{-2}$, $(S/N)_r^{-4}$')
    else:
        ax.set_title(title)
    #ax.set_title(r'Fit parameters $1$, $\sqrt{detC}$, $detC$')

    #fig.set_size_inches(12.,6.0)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

    return m_corr, c_corr

def mean_e_vs_r(r, e, g_app, mask, title='e vs radius', filename=None):
    import matplotlib.pyplot as plt
    import numpy
    plt.style.use('supermongo')
    fig, ax = plt.subplots(1, 1)

    index = [ (r[mask] >= 0.0) & (r[mask] < 0.32),
              (r[mask] >= 0.32) & (r[mask] < 0.40),
              (r[mask] >= 0.40) & (r[mask] < 0.50),
              (r[mask] >= 0.50) & (r[mask] < 0.63),
              (r[mask] >= 0.63) & (r[mask] < 11.20) ]
    nbins = len(index)
    m1 = numpy.empty(nbins)
    m2 = numpy.empty(nbins)
    sigm1 = numpy.empty(nbins)
    sigm2 = numpy.empty(nbins)
    for i in range(nbins):
        m1[i], m2[i], sigm1[i], sigm2[i] = calc_m(e[mask][index[i]], g_app[mask][index[i]])
    print 'm1 = ',m1
    print 'm2 = ',m2
    print 'len = ',[ len(e[index[i]]) for i in range(nbins) ]
    print 'sigm1 = ',sigm1
    print 'sigm2 = ',sigm2

    rvals = numpy.array([ numpy.mean(r[mask][index[i]]) for i in range(nbins) ])
    print 'r = ',rvals
    ax.errorbar(rvals, m1, yerr=sigm1, color='blue', fmt='o', label='m1')
    ax.errorbar(rvals+0.01, m2, yerr=sigm2, color='red', fmt='o', label='m2')
    ax.set_xlabel('Radius')
    ax.set_xlim(0.2,0.9)
    ax.set_ylabel('m')
    ax.set_title(title)
    ax.legend(loc='lower right')

    #fig.suptitle(title, fontsize=20)
    #fig.set_size_inches(12.5,4.0)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.9)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def mean_e_vs_z(z, e, m_corr, c_corr, g_app, e_true, e_psf, w, mask, title='e vs z', filename=None):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy
    plt.style.use('supermongo')
    fig, ax = plt.subplots(3, 1, sharex='col', sharey='row')
    matplotlib.rcParams.update({'font.size': 8})

    zbins = [ 0.3, 0.644, 0.901, 1.3 ]
    #zbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3 ]
    nz = len(zbins) - 1
    m1 = numpy.empty(nz)
    m2 = numpy.empty(nz)
    sigm1 = numpy.empty(nz)
    sigm2 = numpy.empty(nz)
    c1 = numpy.empty(nz)
    c2 = numpy.empty(nz)
    sigc1 = numpy.empty(nz)
    sigc2 = numpy.empty(nz)
    zz = numpy.empty(nz)
    nn = numpy.empty(nz)
    alpha1 = numpy.empty(nz)
    alpha2 = numpy.empty(nz)
    sigalpha1 = numpy.empty(nz)
    sigalpha2 = numpy.empty(nz)

    # Plot m vs z
    for i in range(nz):
        mask2 = mask & (z > zbins[i]) & (z < zbins[i+1])
        print 'z = ',zbins[i],zbins[i+1]
        print 'mean m = ',numpy.mean(m_corr[mask2])
        print 'n = ',numpy.sum(mask2)
        mc = calc_mc(e_true, g_app, w, mask2)
        m1[i], m2[i], c1[i], c2[i], sigm1[i], sigm2[i], sigc1[i], sigc2[i] = mc
        zz[i] = numpy.mean(z[mask2])
        nn[i] = numpy.sum(mask2)
        alpha1[i], _, sigalpha1[i], _ = linear_fit(e_true[mask2].real-g_app[mask2].real, e_psf[mask2].real, w[mask2])
        alpha2[i], _, sigalpha2[i], _ = linear_fit(e_true[mask2].imag-g_app[mask2].real, e_psf[mask2].imag, w[mask2])
    print 'zz = ',zz
    print 'nn = ',nn
    print 'm1 = ',m1
    print 'm2 = ',m2
    print 'sigm1 = ',sigm1
    print 'sigm2 = ',sigm2
    ax[0].errorbar(zz, m1, yerr=sigm1, color='green', fmt='o', label='m1, true', lw=0.3)
    ax[0].errorbar(zz+0.01, m2, yerr=sigm2, color='LimeGreen', fmt='o', label='m2, true', lw=0.3)
    ax[1].errorbar(zz, c1, yerr=sigc1, color='green', fmt='o', label='c1, true', lw=0.3)
    ax[1].errorbar(zz+0.01, c2, yerr=sigc2, color='LimeGreen', fmt='o', label='c2, true', lw=0.3)
    ax[2].errorbar(zz, alpha1, yerr=sigalpha1, color='green', fmt='o', label='alpha1, true')
    ax[2].errorbar(zz+0.01, alpha2, yerr=sigalpha2, color='LimeGreen', fmt='o', label='alpha2, true')
    
    for i in range(nz):
        mask2 = mask & (z > zbins[i]) & (z < zbins[i+1])
        mc = calc_mc(e, g_app, w, mask2)
        m1[i], m2[i], c1[i], c2[i], sigm1[i], sigm2[i], sigc1[i], sigc2[i] = mc
        alpha1[i], _, sigalpha1[i], _ = linear_fit(e[mask2].real-g_app[mask2].real, e_psf[mask2].real, w[mask2])
        alpha2[i], _, sigalpha2[i], _ = linear_fit(e[mask2].imag-g_app[mask2].imag, e_psf[mask2].imag, w[mask2])
    print 'm1 = ',m1
    print 'm2 = ',m2
    print 'sigm1 = ',sigm1
    print 'sigm2 = ',sigm2
    if 'ngmix' not in title:
        ax[0].errorbar(zz, m1, yerr=sigm1, color='red', fmt='o', label='m1, raw')
        ax[0].errorbar(zz+0.01, m2, yerr=sigm2, color='magenta', fmt='o', label='m2, raw')
        ax[1].errorbar(zz, c1, yerr=sigc1, color='red', fmt='o', label='c1, raw')
        ax[1].errorbar(zz+0.01, c2, yerr=sigc2, color='magenta', fmt='o', label='c2, raw')
        ax[2].errorbar(zz, alpha1, yerr=sigalpha1, color='red', fmt='o', label='alpha1, raw')
        ax[2].errorbar(zz+0.01, alpha2, yerr=sigalpha2, color='magenta', fmt='o', label='alpha2, raw')

    for i in range(nz):
        mask2 = mask & (z > zbins[i]) & (z < zbins[i+1])
        means = 1. + numpy.sum(w[mask2] * m_corr[mask2]) / numpy.sum(w[mask2])
        ee = (e[mask2] - c_corr[mask2]) / means
        mc = calc_mc(ee, g_app[mask2], w[mask2])
        m1[i], m2[i], c1[i], c2[i], sigm1[i], sigm2[i], sigc1[i], sigc2[i] = mc
        alpha1[i], _, sigalpha1[i], _ = linear_fit(ee.real-g_app[mask2].real, e_psf[mask2].real, w[mask2])
        alpha2[i], _, sigalpha2[i], _ = linear_fit(ee.imag-g_app[mask2].imag, e_psf[mask2].imag, w[mask2])
    print 'm1 = ',m1
    print 'm2 = ',m2
    print 'sigm1 = ',sigm1
    print 'sigm2 = ',sigm2
    ax[0].errorbar(zz, m1, yerr=sigm1, color='blue', fmt='o', label='m1, corr')
    ax[0].errorbar(zz+0.01, m2, yerr=sigm2, color='cyan', fmt='o', label='m2, corr')
    ax[1].errorbar(zz, c1, yerr=sigc1, color='blue', fmt='o', label='c1, corr')
    ax[1].errorbar(zz+0.01, c2, yerr=sigc2, color='cyan', fmt='o', label='c2, corr')
    ax[2].errorbar(zz, alpha1, yerr=sigalpha1, color='blue', fmt='o', label='alpha1, corr')
    ax[2].errorbar(zz+0.01, alpha2, yerr=sigalpha2, color='cyan', fmt='o', label='alpha2, corr')

    ax[0].plot( [zbins[0],zbins[-1]], [0.,0.], color='k')
    ax[1].plot( [zbins[0],zbins[-1]], [0.,0.], color='k')
    ax[2].plot( [zbins[0],zbins[-1]], [0.,0.], color='k')
    ax[2].set_xlabel('zphot')
    ax[1].set_xlim(zbins[0], zbins[-1])
    ax[0].set_ylabel('m')
    ax[1].set_ylabel('c')
    ax[2].set_ylabel('alpha')
    ax[0].set_title(title)
    ax[0].legend(loc='upper left', fontsize=8)

    plt.tight_layout()
    fig.set_size_inches(6.5,12.0)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def mean_e_vs_er(er, e, g_app, mask, title='<e> vs raw |e|', filename=None):
    import matplotlib.pyplot as plt
    import numpy
    plt.style.use('supermongo')
    fig, ax = plt.subplots(1, 1)

    index = [ (er[mask] >= 0.0) & (er[mask] < 0.1),
              (er[mask] >= 0.1) & (er[mask] < 0.2),
              (er[mask] >= 0.2) & (er[mask] < 0.3),
              (er[mask] >= 0.3) & (er[mask] < 0.4),
              (er[mask] >= 0.4) & (er[mask] < 0.5),
              (er[mask] >= 0.5) & (er[mask] < 0.6),
              (er[mask] >= 0.6) & (er[mask] < 1.0) ]
    nbins = len(index)
    m1 = numpy.empty(nbins)
    m2 = numpy.empty(nbins)
    sigm1 = numpy.empty(nbins)
    sigm2 = numpy.empty(nbins)
    for i in range(nbins):
        m1[i], m2[i], sigm1[i], sigm2[i] = calc_m(e[mask][index[i]], g_app[mask][index[i]])
    print 'm1 = ',m1
    print 'm2 = ',m2
    print 'len = ',[ len(e[index[i]]) for i in range(nbins) ]
    print 'sigm1 = ',sigm1
    print 'sigm2 = ',sigm2

    ervals = numpy.array([ numpy.mean(er[mask][index[i]]) for i in range(nbins) ])
    print 'er = ',ervals
    ax.errorbar(ervals, m1, yerr=sigm1, color='blue', fmt='o', label='m1')
    ax.errorbar(ervals+0.01, m2, yerr=sigm2, color='red', fmt='o', label='m2')
    ax.set_xlabel('Raw |e|')
    ax.set_xlim(0.,0.9)
    ax.set_ylabel('m')
    ax.set_title(title)
    ax.legend(loc='lower right')

    #fig.suptitle(title, fontsize=20)
    #fig.set_size_inches(12.5,4.0)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.9)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def e_vs_epsf(e, e_true, g_app, e_psf, w, mask, title='', filename=None, m_corr=None, c_corr=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    fig, axes = plt.subplots(1, 2)

    nx = 10

    me1 = numpy.zeros(nx)
    me2 = numpy.zeros(nx)
    mepsf = numpy.zeros(nx)
    sig1 = numpy.zeros(nx)
    sig2 = numpy.zeros(nx)
    use = numpy.ones(nx, dtype=bool)

    for (epsf, ax, label) in [ (e_psf.real, axes[0], r'PSF $e_1$'), 
                               (e_psf.imag, axes[1], r'PSF $e_2$') ]:
        sorte = numpy.sort(epsf[mask])
        deciles = numpy.array([ sorte[i*len(sorte)/nx] for i in range(nx) ] + [ sorte[-1] ])
        for j in range(nx):
            sample = numpy.where(mask & (epsf >= deciles[j]) & (epsf < deciles[j+1]))
            if numpy.sum(sample) == 0:
                use[j] = False
                continue
            e_corr = e[sample].copy()
            if c_corr is not None:
                e_corr -= c_corr[sample]
            if m_corr is not None:
                means = 1. + numpy.sum(w[sample] * m_corr[sample]) / numpy.sum(w[sample])
                e_corr /= means
            me1[j] = numpy.sum(w[sample] * e_corr.real) / numpy.sum(w[sample])
            me2[j] = numpy.sum(w[sample] * e_corr.imag) / numpy.sum(w[sample])
            mepsf[j] = numpy.sum(w[sample] * epsf[sample]) / numpy.sum(w[sample])
            sig1[j] = numpy.sum(w[sample]**2 * (e_corr.real-me1[j])**2) / numpy.sum(w[sample])**2
            sig2[j] = numpy.sum(w[sample]**2 * (e_corr.imag-me2[j])**2) / numpy.sum(w[sample])**2

        sig1 = numpy.sqrt(sig1)
        sig2 = numpy.sqrt(sig2)

        ax.errorbar(mepsf[use], me1[use], yerr=sig1[use], color='red', fmt='o', label=r'$\langle e_1 \rangle$')
        ax.errorbar(mepsf[use], me2[use], yerr=sig2[use], color='blue', fmt='o', label=r'$\langle e_2 \rangle$')
        ax.plot( [mepsf[0],mepsf[-1]], [0.,0.], color='k')
        ax.plot( [mepsf[0],mepsf[-1]], [0.,0.], color='k')
        ax.set_ylabel(r'$\langle e \rangle$')
        ax.set_xlabel(label)

        alpha, ca, sigalpha, _ = linear_fit(e[mask].real-g_app[mask].real, epsf[mask], w[mask])
        beta, cb, sigbeta, _ = linear_fit(e[mask].imag-g_app[mask].imag, epsf[mask], w[mask])

        ax.plot( [mepsf[0],mepsf[-1]], [alpha*mepsf[0]+ca,alpha*mepsf[-1]+ca], color='red', lw=0.5, ls='-')
        ax.plot( [mepsf[0],mepsf[-1]], [beta*mepsf[0]+cb,beta*mepsf[-1]+cb], color='blue', lw=0.5, ls='-')

        if 'e_1' in label:
            ax.set_title(r'$\alpha = %f \pm %f$, $\beta = %f \pm %f$'%(alpha,sigalpha,beta,sigbeta))
        else:
            ax.set_title(r'$\alpha = %f \pm %f$, $\beta = %f \pm %f$'%(beta,sigbeta,alpha,sigalpha))
        ax.legend(loc='upper left')

    #axes[0].set_title(title, fontsize=20)
    #axes[0].legend(loc='upper left')
    fig.set_size_inches(16.5,8.5)
    plt.tight_layout()

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def niter_vs_r_e(r, e, niter, mask, title='N Levmar Iterations vs r, e', filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    fig, ax = plt.subplots(1, 4, sharex='col', sharey='row')
    abse = numpy.abs(e[mask])
    ax[0].scatter(r[mask], niter[mask], s=0.1)
    ax[1].scatter(e[mask].real, niter[mask], s=0.1, c='b')
    ax[1].scatter(e[mask].imag, niter[mask], s=0.1, c='r')
    ax[2].scatter(abse, niter[mask], s=0.1)
    #ax[0].set_ylim(0, max(niter[mask])*1.05)
    ax[0].set_xlim(0, max(r[mask])*1.05)
    ax[1].set_xlim(-1., 1.)
    e_range = max(abse) - min(abse)
    ax[2].set_xlim(min(abse) - e_range*0.05, max(abse) + e_range*0.05)
    ax[0].set_xlabel('radius')
    ax[1].set_xlabel('$e_1$, $e_2$')
    ax[2].set_xlabel('$|e|$')
    ax[0].set_ylabel('N Levmar Iterations')
    fig.suptitle(title, fontsize=20)
    fig.set_size_inches(12.5,4.0)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def hist_x(x, mask, label, title='', filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    #fig, ax = plt.subplots(1, 1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, p = ax.hist(x[mask], bins=1000)
    ax.set_title(title)
    ax.set_xlabel(label)
    ax.set_ylim(0, numpy.max(n)*1.1)
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


def hist_e(e, niter, mask, title='Distribution of |e|', filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    #fig, ax = plt.subplots(1, 1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, p = ax.hist(numpy.abs(e[mask]), bins=100,
            histtype='step', color='black', normed=1,
            facecolor='cyan', fill=True, label='all')
    print n
    print bins
    if sum(niter[mask] >= 50) > 0:
        n, bins, p = ax.hist(numpy.abs(e[mask & (niter >= 50)]), bins=bins,
                histtype='step', color='black', normed=1, alpha=0.5,
                facecolor='red', fill=True, label=r'niter $>=$ 50')
        ax.legend()
    ax.set_title(title)
    ax.set_xlabel(r'$|e|$')
    ax.set_ylabel('relative frequency')
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def hist_e_z(e, z, mask, filename=None):
    import numpy
    import matplotlib.pyplot as plt
    plt.style.use('supermongo')
    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row')

    bins = 50
    for ax, zmin, zmax in ( (axes[0,0], 0.3, 0.644),
                            (axes[0,1], 0.644, 0.901),
                            (axes[1,0], 0.901, 1.3),
                            (axes[1,1], 0.3, 1.3) ):

        mask2 = mask & (z >= zmin) & (z < zmax)
        n, bins, p = ax.hist(numpy.abs(e[mask2]), bins=bins,
                             histtype='step', color='black', normed=1,
                             facecolor='cyan', fill=True)
        ax.set_title('$%f < z < %f$'%(zmin,zmax))

    axes[1,0].set_xlabel(r'$|e|$')
    axes[1,0].set_ylabel('relative frequency')
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def snr_vs_trcov(snr, trcov, mask, filename=None):
    import matplotlib.pyplot as plt
    import numpy

    fig, ax = plt.subplots(1, 1)

    mask2 = mask & (trcov < 1) & (snr < 1000) & (snr > 3)
    ax.hexbin(snr[mask2], trcov[mask2], xscale='log', yscale='log')
    ax.set_xlabel('$(S/N)_w$')
    ax.set_ylabel('$Tr(Cov(e1,e2))$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def comparison_histograms(e_im, snr_im, rgp_im, r_im, disc_im, e_sv, snr_sv, rgp_sv, r_sv, disc_sv, 
                          filename1, filename2, logsn=True):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy

    fig, axes = plt.subplots(1, 3)

    ax = axes[0]
    e_im = numpy.abs(e_im)
    #e_im = e_im[e_im<0.9]
    e_sv = numpy.abs(e_sv)
    #e_sv = e_sv[e_sv<0.9]
    n, bins, patches = ax.hist(e_im, bins=500, histtype='step', color='red', normed=1, label='GREAT-DES')
    ax.hist(e_sv, bins=bins, histtype='step', color='blue', normed=1, label='SV data')
    ax.set_xlabel(r'$|e|$')
    ax.set_ylabel(r'$N$')
    ax.set_yticks([]) 
    ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8]) 
    ax.set_xlim(0,0.9)
    ax.legend(loc='upper right', fontsize=10)

    ax = axes[1]
    n, bins, patches = ax.hist(rgp_im, bins=500, histtype='step', color='red', normed=1, label='GREAT-DES')
    ax.hist(rgp_sv, bins=bins, histtype='step', color='blue', normed=1, label='SV data')
    ax.set_xlabel(r'$R_{gp}/R_p$')
    ax.set_ylabel(r'$N$')
    ax.set_yticks([]) 
    ax.set_xticks([1.2, 1.4, 1.6, 1.8, 2.0]) 
    ax.set_xlim(1.2,2.)
    #ax.axvline(1.2, color='blue', linestyle='--')
    ax.legend(loc='upper right', fontsize=10)

    ax = axes[2]
    m = snr_im < 100
    n, bins, patches = ax.hist(snr_im[m], bins=500, histtype='step', color='red', normed=1, label='GREAT-DES')
    m = snr_sv < 100
    ax.hist(snr_sv[m], bins=bins, histtype='step', color='blue', normed=1, label='SV data')
    ax.set_xlabel(r'$(S/N)_w$')
    ax.set_ylabel(r'$N$')
    ax.set_yticks([]) 
    if logsn:
        ax.set_xlim(15.,100.)
        ax.set_xscale('log')
    else:
        ax.set_xlim(15.,100.)
    ax.set_xticks([20, 50, 100]) 
    # cf. http://stackoverflow.com/questions/14530113/set-ticks-with-logarithmic-scale
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.axvline(15, color='blue', linestyle='--')
    ax.legend(loc='upper right', fontsize=10)

    #print fig.get_size_inches()  # Default is apparently 8,6
    fig.set_size_inches(12.,4.)
    plt.tight_layout()
    plt.savefig(filename1)
    fig, axes = plt.subplots(1, 3)

    ax = axes[0]
    plot_mean_x_vs_y(ax, snr_im, r_im, color='red', label='GREAT-DES', line=True, nx=50)
    plot_mean_x_vs_y(ax, snr_sv, r_sv, color='blue', label='SV data', line=True, nx=50)
    ax.set_xlabel(r'$(S/N)_w$')
    ax.set_ylabel(r'$\langle R_g \rangle$ (arcsec)')
    ax.set_ylim(0.4,0.9)
    ax.set_yticks([0.4, 0.6, 0.8]) 
    if logsn:
        ax.set_xlim(15.,100.)
        ax.set_xscale('log')
    else:
        ax.set_xlim(15.,100.)
    ax.set_xticks([20, 50, 100]) 
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.axvline(15, color='blue', linestyle='--')
    ax.legend(loc='upper right', fontsize=10)

    ax = axes[1]
    plot_mean_x_vs_y(ax, snr_im, rgp_im, color='red', label='GREAT-DES', line=True, nx=50)
    plot_mean_x_vs_y(ax, snr_sv, rgp_sv, color='blue', label='SV data', line=True, nx=50)
    ax.set_xlabel(r'$(S/N)_w$')
    ax.set_ylabel(r'$\langle R_{gp}/R_p \rangle$')
    ax.set_ylim(1.3,1.6)
    ax.set_yticks([1.3, 1.4, 1.5, 1.6]) 
    if logsn:
        ax.set_xlim(15.,100.)
        ax.set_xscale('log')
    else:
        ax.set_xlim(15.,100.)
    ax.set_xticks([20, 50, 100]) 
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.axvline(15, color='blue', linestyle='--')
    ax.legend(loc='upper right', fontsize=10)

    ax = axes[2]
    plot_mean_x_vs_y(ax, snr_im, 1.-disc_im, color='red', label='GREAT-DES', line=True, nx=50)
    plot_mean_x_vs_y(ax, snr_sv, 1.-disc_sv, color='blue', label='SV data', line=True, nx=50)
    ax.set_xlabel(r'$(S/N)_w$')
    ax.set_ylabel(r'Bulge fraction')
    ax.set_ylim(0.,0.6)
    ax.set_yticks([0., 0.2, 0.4, 0.6]) 
    if logsn:
        ax.set_xlim(15.,100.)
        ax.set_xscale('log')
    else:
        ax.set_xlim(15.,100.)
    ax.set_xticks([20, 50, 100]) 
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.axvline(15, color='blue', linestyle='--')
    ax.legend(loc='upper right', fontsize=10)

    fig.set_size_inches(12.,4.)
    plt.tight_layout()
    plt.savefig(filename2)


def build_m_vs_z_ngmix():
    import numpy
    import tests
    # Pull out the complete script for building the m_vs_z plot for ngmix.

    # Read in the ngmix data
    num_ng, t_ng, flux_ng, snr_ng, tsnr_ng, flag_ng, tr_ng, snrr_ng, flagr_ng, e_ng, epsf_ng, tpsf_ng, sens_ng, chisq_ng, trcov_ng = tests.load_ng_data()

    # Read in the truth data
    num_true, r_true, f_true, e_true, e_raw, g_app, id_cosmos, srcn, z, use = tests.load_truth()

    # Find the index into truth values that refer to the corresponding ngmix values
    # i.e. num_true[ng] == num_ng
    ng = numpy.in1d(num_true, num_ng)
    assert numpy.all(num_true[ng] == num_ng)

    # Select good objects
    good_ng = (flag_ng == 0) & (flagr_ng == 0) & (sens_ng.real > 0) & (sens_ng.imag > 0) & (use[ng])

    # Convert ngmix sensitivity into m,c terminology.
    c0 = numpy.zeros(len(e_ng))
    m_ng = (sens_ng.real + sens_ng.imag)/2. - 1.

    # Calculate the ngmix weights.
    w_ng = 1. / (2 * 0.22**2 + trcov_ng)

    # Define the mask that we apply in the data
    mask_ng = good_ng & (snrr_ng > 15) & (tr_ng/tpsf_ng > 0.15)

    # Plot m vs z.  (Wrongly called e here, but that's what it means.)
    tests.mean_e_vs_z(z[ng], e_ng, m_ng, c0, g_app[ng], e_true[ng], epsf_ng, w_ng, mask_ng, title=r'ngmix $(S/N)_r > 15$, $Tr/Tp > 0.15$, normal $w$', filename='mvsz_ngmix.pdf')

    # Make unweighted version
    w1 = numpy.ones(len(e_ng))
    tests.mean_e_vs_z(z[ng], e_ng, m_ng, c0, g_app[ng], e_true[ng], epsf_ng, w1, mask_ng, title=r'ngmix $(S/N)_r > 15$, $Tr/Tp > 0.15$, $w=1$', filename='mvsz_ngmix_unweighted.pdf')

    # Load im3shape catalog
    num_im, r_im, rgp_im, f_im, snrw_im, flag_im, e_im, snrr_im, disc_im, niter_im, chisq_im, minres_im, maxres_im, info_im, epsf_im, fwhmpsf_im, trcov_im, detcov_im, varr_im, vare_im = tests.load_im_data()
    im = numpy.in1d(num_true, num_im)
    assert numpy.all(num_true[im] == num_im)
    good_im = (flag_im == 0) & (info_im == 0)
    mask_im = good_im & (snrw_im > 15) & (rgp_im > 1.2)

    # Find indices between ng and im catalogs
    # im_ng gives the index into the im catalog to get objects also in ng
    im_ng = numpy.in1d(num_im, num_ng)
    # ng_im gives the index into the ng catalog to get objects also in im
    ng_im = numpy.in1d(num_ng, num_im)
    assert numpy.all(num_im[im_ng] == num_ng[ng_im])

    # mask_c is the combined mask.  Can be used with any of the catalogs by first
    # applying the appropriate mask to get objects that are in both.
    assert len(mask_ng[ng_im]) == len(mask_im[im_ng])
    mask_c = mask_ng[ng_im] & mask_im[im_ng]
    assert len(mask_c) == len(num_true[ng & im])
    assert len(mask_c) == len(num_ng[ng_im])
    assert len(mask_c) == len(num_im[im_ng])

    # m_true for the matched catalog:
    print 'm_true for matched catalog: ',tests.calc_mc(e_true[ng&im][mask_c], g_app[ng&im][mask_c])[:2]
    s_match = 1. + numpy.sum(w_ng[ng_im] * m_ng[ng_im]) / numpy.sum(w_ng[ng_im])
    print 'm_ngmix for matched catalog: ',tests.calc_mc(e_ng[ng_im][mask_c]/s_match, g_app[ng&im][mask_c], w=w_ng[ng_im])[:2]
    s_ave = 1. + numpy.sum(w_ng[mask_ng] * m_ng[mask_ng]) / numpy.sum(w_ng[mask_ng])
    print 'm_ngmix ngmix cuts: ',tests.calc_mc(e_ng[mask_ng]/s_ave, g_app[ng][mask_ng], w_ng[mask_ng])[:2]

    # m vs z for ngmix shapes on the matched catalog
    tests.mean_e_vs_z(z[ng&im], e_ng[ng_im], m_ng[ng_im], c0[ng_im], g_app[ng&im], e_true[ng&im], epsf_ng[ng_im], w_ng[ng_im], mask_c, title=r'ngmix shears on matched selection, normal $w$', filename='mvsz_ngmix_match.pdf')
    tests.mean_e_vs_z(z[ng&im], e_ng[ng_im], m_ng[ng_im], c0[ng_im], g_app[ng&im], e_true[ng&im], epsf_ng[ng_im], w1[ng_im], mask_c, title=r'ngmix shears on matched selection, $w=1$', filename='mvsz_ngmix_match_unweighted.pdf')


def main():
    import numpy
    num_im, r_im, rgp_im, f_im, snrw_im, flag_im, e_im, snrr_im, disc_im, niter_im, chisq_im, minres_im, maxres_im, info_im, epsf_im, fwhmpsf_im, trcov_im, detcov_im, varr_im, vare_im = tests.load_im_data()
    num_ng, t_ng, flux_ng, snr_ng, tsnr_ng, flag_ng, tr_ng, snrr_ng, flagr_ng, e_ng, epsf_ng, tpsf_ng, sens_ng, chisq_ng, trcov_ng = tests.load_ng_data()
    num_true, r_true, f_true, e_true, e_raw, g_app, id_cosmos, srcn, z, use = tests.load_truth()
    great3, w_great3 = tests.load_great3(id_cosmos)
    ng = numpy.in1d(num_true, num_ng)
    im = numpy.in1d(num_true, num_im)

    bade = numpy.abs(e_im) > numpy.max(numpy.abs(e_im)) * 0.999
    good_im = (flag_im == 0) & (info_im == 0) & ~bade
    snr15_im = good_im & (snrr_im > 15)
    tests.check_meas(r_true[im], r_im, e_true[im], e_im, snr15_im, filename='check.pdf')

    good_ng = (flag_ng == 0) & (flagr_ng == 0) & (sens_ng.real > 0) & (sens_ng.imag > 0)
    snr15_ng = good_ng & (snrr_ng > 15)
    tests.check_meas(r_true[ng], t_ng, e_true[ng], e_ng, snr15_ng, filename='check.pdf')

    tests.mean_e_vs_r(r_im, e_im, g_app, snr15_im, title=r'$S/N > 15$, Measured shapes', filename='evsr_snr15_im.pdf')
    tests.mean_e_vs_r(r_im, e_true, g_app, snr15_im, title=r'$S/N > 15$, True shapes', filename='evsr_snr15_true.pdf')

    snr10 = (flag_im == 0) & (snrr_im > 10)
    tests.niter_vs_r_e(r_im, e_im, niter, snr10, r'Measured size, shape, $(S/N)_r > 10$', filename='meas_snrgt10.png')

    info2 = info_im & ~(16 | 256 | 16384) # No rgp or S/N cuts
    mask_im = good_im & (snrw_im > 15) & (rgp_im > 1.2)
    w_im = 1. / ( 2.*0.22**2 + 2./snrr_im**2)

    global_alpha = tests.linear_fit(e_true[mask]-g_app[mask], e_psf[mask])[0].real
    global_m = tests.linear_fit(e_true[mask]-g_app[mask], g_app[mask])[0].real

    m_corr, c_corr = tests.nbc(e_im, e_true[im], g_app[im], epsf_im, rgp_im, snrw_im, disc_im, w_im, mask_im, title='NBC with $(S/N)_w$, with $(S/N)_r$ cut', filename='nbc.pdf')
    tests.mean_e_vs_z(z[im], e_im, m_corr, c_corr, g_app[im], e_true[im], epsf_im, w_im, mask_im, title=r'im3shape $(S/N)_r > 15$, $rgp > 1.2$', filename='evsz.pdf')

    c0 = numpy.zeros(len(e_ng))
    m_ng = (sens_ng.real + sens_ng.imag)/2. - 1.
    w_ng = 1. / (2 * 0.22**2 + trcov_ng)
    mask_ng = good_ng & (snrr_ng > 15) & (tr_ng/tpsf_ng > 0.15)
    tests.mean_e_vs_z(z[ng], e_ng, m_ng, c0, g_app[ng], e_true[ng], epsf_ng, w_ng, mask_ng, title=r'ngmix $(S/N)_r > 15$, $Tr/Tp > 0.15$', filename='evsz.pdf')

    im_ng = numpy.in1d(num_im, num_ng)
    ng_im = numpy.in1d(num_ng, num_im)
    mask_c = mask_ng[ng_im] & mask_im[im_ng]
    # m,c for the matched catalog:
    tests.calc_mc(e_true[ng&im][mask_c], g_app[ng&im][mask_c])
    # As a function of z:
    tests.calc_mc(e_true[ng&im&(z<0.8)][mask_c[z[ng&im]<0.8]], g_app[ng&im&(z<0.8)][mask_c[z[ng&im]<0.8]])

def figure11():
    # Make figure 11
    num_im, r_im, rgp_im, f_im, snrw_im, flag_im, e_im, snrr_im, disc_im, niter_im, chisq_im, minres_im, maxres_im, info_im, epsf_im, fwhmpsf_im, trcov_im, detcov_im, varr_im, vare_im = tests.load_im_data()
    r_sv, rgp_sv, snrw_sv, flag_sv, e_sv, disc_sv, info_sv = tests.load_im_sv()
    f_im = (flag_im == 0) & (info_im == 0)
    f_sv = (flag_sv == 0) & (info_sv == 0)
    f_im = f_im & (snrw_im > 15) & (rgp_im > 1.2)
    f_sv = f_sv & (snrw_sv > 15) & (rgp_sv > 1.2)
    tests.comparison_histograms(e_im[f_im], snrw_im[f_im], rgp_im[f_im], r_im[f_im], disc_im[f_im], e_sv[f_sv], snr_sv[f_sv], rgp_sv[f_sv], r_sv[f_sv], disc_sv[f_sv], 'h1.pdf', 'h2.pdf')



if __name__ == '__main__':
    #main()
    #build_m_vs_z_ngmix()
    figure11()
