#! /usr/bin/env python
from __future__ import print_function
import numpy as np
import os
from read_psf_cats import read_data

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description='Make whisker plots')

    # Drectory arguments
    parser.add_argument('--work', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--file', default='',
                        help='list of exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Options
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--bands', default='grizY', type=str,
                        help='Limit to the given bands')
    parser.add_argument('--use_psfex', default=False, action='store_const', const=True,
                        help='Use PSFEx rather than Piff model')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')

    args = parser.parse_args()
    return args


def get_ccdnums_from_file(file_names):
    return np.array([ 0 ] + [ int(f[-10:-8]) for f in file_names[1:] ])

def get_ngmix_epoch_data(use_gold=True):

    import astropy.io.fits as pyfits
    import glob
    import json

    #band_list = []
    ccd_list = []
    x_list = []
    y_list = []

    e1_list = []
    e2_list = []
    s_list = []
    w_list = []

    if use_gold:
        info_file = '/astro/u/mjarvis/work/des_sv_wl_info.fits'
        with pyfits.open(info_file) as f:
            fcat = f[1].data
        gold = fcat['sva1_gold_flags'] == 0
        spte = fcat['sva1_spte_flags'] == 0
        mag = fcat['sva1_gold_mag_flags'] == 0
        ngmix = fcat['ngmix_flags'] == 0
        #im3shape = fcat['im3shape_flags'] == 0
        all_ng = gold & spte & mag & ngmix
        #all_im = gold & spte & mag & im3shape
        work = '/astro/u/mjarvis/work/ngmix011_collated'
    else:
        work = '/astro/u/mjarvis/work/ngmix011_collated/no_info'

    ngmix_dir = '/astro/u/astrodat/data/DES/wlpipe/ngmix011/collated'
    print('ngmix_dir = ',ngmix_dir)

    file_list = glob.glob(os.path.join(ngmix_dir,'*blind.fits'))

    force = False
    #force = True

    for i,file in enumerate(file_list):

        print('Start work on ',file)
        print(i,'/',len(file_list))
        json_file = os.path.join(work, os.path.basename(file)[:-5] + ".json")
        print('json_file = ',json_file)

        done = False
        if not force and os.path.exists(json_file):
            try:
                with open(json_file,'r') as f:
                    json_data = json.load(f)
                print('loaded from json file')
                print('len(data) = ',len(json_data))
                band, ccdnum, x, y, e1, e2, s, w = json_data
                print('unpacked json data')
                done = True
            except Exception as e:
                print('caught ',e)
                pass

        if not done:

            try:
                print('read ',file)
                with pyfits.open(file) as f:
                    fits = f['model_fits'].data
                    epoch = f['epoch_data'].data
                    meta = f['meta_data'].data
            except:
                continue

            if False:
                try:
                    round_file = os.path.join(ngmix_dir, 'make-round/v6/output',
                                            os.path.basename(file)[:-10] + 'round.fits')
                    print('read ',round_file)
                    with pyfits.open(round_file) as f:
                        round = f[1].data
                except:
                    continue
                assert len(fits) == len(round)

            meds_files = meta['meds_file']
            ccdnums = []
            for m in meds_files:
                with pyfits.open(m) as f:
                    ccdnums.append(get_ccdnums_from_file(f['image_info'].data['image_path']))
            print('ccdnums = ',ccdnums)

            if use_gold:
                use1 = np.in1d(fits['id'], fcat['coadd_objects_id'][all_ng])
            else:
                use1 = ( (fits['flags'] == 0) & (fits['exp_flags'] == 0) )
                #use1 = ( (fits['flags'] == 0) & (round['round_flags'] == 0) &
                         #(fits['exp_flags'] == 0) & (round['exp_round_flags'] == 0) &
                         #(0.4 < fits['exp_arate']) & (fits['exp_arate'] < 0.6) &
                         #(round['exp_s2n_r'] > 15) & (round['exp_T_r']/fits['psfrec_T'] > 0.15) &
                         #(fits['exp_g_sens'][:,0] > 0.0) & (fits['exp_g_sens'][:,1] > 0.0)
                       #)

            print('nobject = ',len(fits))
            print('nuse = ',sum(use1))

            e1 = fits['exp_g'][:,0]
            e1 = -e1 # Basically, that's the effect of the gross properties of the WCS.
                     # cf. http://ast.noao.edu/sites/default/files/NOAO_DHB_v2.2.pdf, p. 4-4.
            e2 = fits['exp_g'][:,1]
            s = fits['exp_T']
            print('e1 = ',e1)

            # w = (2 SN^2 + Tr(cov))^-1
            w = 2. * 0.22**2 + fits['exp_g_cov'][:,0,0] + fits['exp_g_cov'][:,1,1]
            w = 1/w
            print('w = ',w)

            band = epoch['band']
            print('band = ',band)
            band_num = epoch['band_num']
            print('band_num = ',band_num)
            x = epoch['orig_col']
            y = epoch['orig_row']
            print('x,y = ',x,y)
            index = epoch['number']-1
            print('index = ',index)
            ccdnum = np.empty(len(epoch))
            for b in range(len(ccdnums)):
                mask = band_num == b
                ccdnum[mask] = ccdnums[b][epoch['file_id'][mask]]
            print('ccdnum = ',ccdnum)
            use2 = use1[index] & (ccdnum > 0)
            print('use2 = ',use2)

            band = band[use2]
            ccdnum = ccdnum[use2]
            x = x[use2]
            y = y[use2]
            e1 = e1[index[use2]]
            e2 = e2[index[use2]]
            s = s[index[use2]]
            w = w[index[use2]]

            # Write this to a json file for quicker load next time.
            with open(json_file,'w') as f:
                json.dump( (band.tolist(), ccdnum.tolist(), x.tolist(), y.tolist(),
                            e1.tolist(), e2.tolist(), s.tolist(), w.tolist()), f)

        #band_list.append(np.array(band)
        ccd_list.append(np.array(ccdnum,dtype=np.int16))
        x_list.append(np.array(x,dtype=np.float32))
        y_list.append(np.array(y,dtype=np.float32))

        e1_list.append(np.array(e1,dtype=np.float32))
        e2_list.append(np.array(e2,dtype=np.float32))
        s_list.append(np.array(s,dtype=np.float32))
        w_list.append(np.array(w,dtype=np.float32))


    print('\nFinished processing all exposures')

    #band = np.concatenate(band_list)
    ccd = np.concatenate(ccd_list)
    x = np.concatenate(x_list)
    y = np.concatenate(y_list)

    e1 = np.concatenate(e1_list)
    e2 = np.concatenate(e2_list)
    s = np.concatenate(s_list)
    w = np.concatenate(w_list)

    print('Done concatenating.')

    return (ccd, x, y, e1, e2, s, w)


def get_im3shape_epoch_data(use_gold=True):

    import astropy.io.fits as pyfits
    import glob
    import json

    ccd_list = []
    x_list = []
    y_list = []

    e1_list = []
    e2_list = []
    s_list = []
    w_list = []

    if use_gold:
        info_file = '/astro/u/mjarvis/work/des_sv_wl_info.fits'
        with pyfits.open(info_file) as f:
            fcat = f[1].data
        gold = fcat['sva1_gold_flags'] == 0
        spte = fcat['sva1_spte_flags'] == 0
        mag = fcat['sva1_gold_mag_flags'] == 0
        #ngmix = fcat['ngmix_flags'] == 0
        im3shape = fcat['im3shape_flags'] == 0
        print('gold: ',np.sum(gold))
        print('spte: ',np.sum(spte))
        print('mag: ',np.sum(mag))
        print('im3shape: ',np.sum(im3shape))

        #all_ng = gold & spte & mag & ngmix
        all_im = gold & spte & mag & im3shape
        print('all_im: ',np.sum(all_im))
        work = '/astro/u/mjarvis/work/im3shape_v9_epoch'
    else:
        work = '/astro/u/mjarvis/work/im3shape_v9_epoch/no_info'

    im3shape_dir = '/astro/u/astrodat/data/DES/wlpipe/im3shape_v9/full_cats/'
    #epoch_dir = '/astro/u/astrodat/data/DES/wlpipe/im3shape_v9/epoch_cats/'
    epoch_dir = '/astro/u/astrodat/data/DES/wlpipe/im3shape_v8/epoch_cats/'
    print('im3shape_dir = ',im3shape_dir)

    file_list = glob.glob(os.path.join(im3shape_dir,'*.fits'))

    force = False
    #force = True

    for i,file in enumerate(file_list):

        print('Start work on ',file)
        print(i,'/',len(file_list))
        json_file = os.path.join(work, os.path.basename(file)[:-5] + ".json")
        print('json_file = ',json_file)

        epoch_file = os.path.join(epoch_dir, os.path.basename(file))
        print('epoch_file = ',epoch_file)

        done = False
        if not force and os.path.exists(json_file):
            try:
                with open(json_file,'r') as f:
                    json_data = json.load(f)
                print('loaded from json file')
                print('len(data) = ',len(json_data))
                ccdnum, x, y, e1, e2, s, w = json_data
                print('unpacked json data')
                done = True
            except Exception as e:
                print('caught ',e)
                pass

        if not done:

            try:
                with pyfits.open(file) as f:
                    data = f[1].data
                with pyfits.open(epoch_file) as f:
                    epoch = f[1].data
            except:
                continue

            print('data ids = ',data['coadd_objects_id'])
            print('epoch ids = ',epoch['coadd_objects_id'])
            print('num in epoch = ',np.sum(np.in1d(data['coadd_objects_id'], epoch['coadd_objects_id'])))
            print('len(data) = ',len(data))
            print('len(epoch) = ',len(epoch))

            # This is only needed because we are using the v8 epoch catalog with v9 data.
            use_epoch = np.in1d(epoch['coadd_objects_id'], data['coadd_objects_id'])
            epoch = epoch[use_epoch]
            use_data = np.in1d(data['coadd_objects_id'], epoch['coadd_objects_id'])
            data = data[use_data]
            print('after matching data to epoch:')
            print('len(data) => ',len(data))
            print('len(epoch) => ',len(epoch))

            if use_gold:
                print('fcat ids = ',fcat['coadd_objects_id'])
                print('num in fcat = ',np.sum(np.in1d(data['coadd_objects_id'], fcat['coadd_objects_id'][all_im])))
                use1 = np.in1d(data['coadd_objects_id'], fcat['coadd_objects_id'][all_im])
            else:
                use1 = ( (data['error_flag'] == 0) &
                         (data['info_flag'] == 0) &
                         (data['snr'] > 15) &
                         (data['mean_rgpp_rp'] > 1.2)
                       )

            print('nuse = ',np.sum(use1))

            e1 = data['e1']
            e1 = -e1 # WCS
            e2 = data['e2']
            s = data['radius']
            #w = data['w']
            w = np.ones(len(e1))  # No weights in v9.

            x = epoch['orig_col']  # These are currently v8 positions.  v9 was missing these cols.
            y = epoch['orig_row']
            print('x,y = ',x,y)
            ccdnum = epoch['ccd']
            print('ccdnum = ',ccdnum)
            index = np.searchsorted(data['coadd_objects_id'], epoch['coadd_objects_id'])
            print('index = ',index)
            print('len index = ',len(index))
            use2 = use1[index] & (ccdnum > 0)
            print('use2 = ',use2)
            print('len use2 = ',len(use2))
            print('len epoch[use2] = ',len(epoch[use2]))
            # use1 is a mask for data
            # use2 is a mask for epoch
            # data[index] gives the right information for corresponding rows in epoch

            ccdnum = ccdnum[use2]
            x = x[use2]
            y = y[use2]
            e1 = e1[index[use2]]
            e2 = e2[index[use2]]
            s = s[index[use2]]
            w = w[index[use2]]

            # Write this to a json file for quicker load next time.
            with open(json_file,'w') as f:
                json.dump( (ccdnum.tolist(), x.tolist(), y.tolist(),
                            e1.tolist(), e2.tolist(), s.tolist(), w.tolist()), f)

        ccd_list.append(np.array(ccdnum,dtype=np.int16))
        x_list.append(np.array(x,dtype=np.float32))
        y_list.append(np.array(y,dtype=np.float32))

        e1_list.append(np.array(e1,dtype=np.float32))
        e2_list.append(np.array(e2,dtype=np.float32))
        s_list.append(np.array(s,dtype=np.float32))
        w_list.append(np.array(w,dtype=np.float32))


    print('\nFinished processing all exposures')

    ccd = np.concatenate(ccd_list)
    x = np.concatenate(x_list)
    y = np.concatenate(y_list)

    e1 = np.concatenate(e1_list)
    e2 = np.concatenate(e2_list)
    s = np.concatenate(s_list)
    w = np.concatenate(w_list)

    print('Done concatenating.')

    return (ccd, x, y, e1, e2, s, w)



def psf_resid(m, de1, de2, dT, key=None):

    import matplotlib
    matplotlib.use('Agg') # needs to be done before import pyplot
    import matplotlib.pyplot as plt
    plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')

    # Bin by mag
    mag_bins = np.linspace(10,17,71)
    print('mag_bins = ',mag_bins)

    index = np.digitize(m, mag_bins)
    print('len(index) = ',len(index))
    bin_de1 = [de1[index == i].mean() for i in range(1, len(mag_bins))]
    print('bin_de1 = ',bin_de1)
    bin_de2 = [de2[index == i].mean() for i in range(1, len(mag_bins))]
    print('bin_de2 = ',bin_de2)
    bin_dT = [dT[index == i].mean() for i in range(1, len(mag_bins))]
    print('bin_dT = ',bin_dT)
    bin_de1_err = [ np.sqrt(de1[index == i].var() / len(de1[index == i]))
                    for i in range(1, len(mag_bins)) ]
    print('bin_de1_err = ',bin_de1_err)
    bin_de2_err = [ np.sqrt(de2[index == i].var() / len(de2[index == i]))
                    for i in range(1, len(mag_bins)) ]
    print('bin_de2_err = ',bin_de2_err)
    bin_dT_err = [ np.sqrt(dT[index == i].var() / len(dT[index == i]))
                    for i in range(1, len(mag_bins)) ]
    print('bin_dT_err = ',bin_dT_err)

    # Fix up nans
    for i in range(1,len(mag_bins)):
        if i not in index:
            bin_de1[i-1] = 0.
            bin_de2[i-1] = 0.
            bin_dT[i-1] = 0.
            bin_de1_err[i-1] = 0.
            bin_de2_err[i-1] = 0.
            bin_dT_err[i-1] = 0.
    print('fixed nans')

    print('index = ',index)
    print('bin_de1 = ',bin_de1)
    print('bin_de2 = ',bin_de2)
    print('bin_dT = ',bin_dT)

    plt.clf()
    #plt.title('PSF Size residuals')
    plt.xlim(10,16.5)
    if key and 'orig' in key:
        plt.ylim(0.3,0.7)
        plt.ylabel(r'$size_{\rm psf}$ (arcsec)')
    elif key and 'model' in key:
        plt.ylim(0.3,0.7)
        plt.ylabel(r'$size_{\rm model}$ (arcsec)')
    else:
        plt.ylim(-0.002,0.012)
        plt.ylabel(r'$size_{\rm psf} - size_{\rm model}$ (arcsec)')
    plt.plot([10,16.5], [0,0], color='black')
    plt.plot([13,13],[-1,1], color='red')
    plt.fill( [10,10,13,13], [-1,1,1,-1], fill=False, hatch='/', color='red')
    plt.errorbar(mag_bins[:-1], bin_dT, yerr=bin_dT_err, color='blue', fmt='o')
    plt.xlabel('Uncorrected Magnitude')
    #plt.ylabel(r'$size_{\rm psf} - \langle size_{\rm psf} \rangle$ (arcsec)')

    plt.clf()
    #plt.title('PSF Ellipticity residuals')
    plt.xlim(10,16.5)
    plt.ylim(-3.e-4,8.e-4)
    plt.plot([10,16.5], [0,0], color='black')
    plt.plot([13,13],[-1,1], color='red')
    plt.fill( [10,10,13,13], [-1,1,1,-1], fill=False, hatch='/', color='red')
    e1_line = plt.errorbar(mag_bins[:-1], bin_de1, yerr=bin_de1_err, color='red', fmt='o')
    e2_line = plt.errorbar(mag_bins[:-1], bin_de2, yerr=bin_de2_err, color='blue', fmt='o')
    plt.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
    plt.xlabel('Uncorrected Magnitude')
    #plt.ylabel(r'$e_{\rm psf} - \langle e_{\rm psf} \rangle$')
    plt.ylabel(r'$e_{\rm psf} - e_{\rm model}$')
    #plt.tight_layout()
    plt.savefig('fig6.eps')


def bin_by_fov(ccd, x, y, e1, e2, s, w=None, nwhisk=5):
    from toFocal import toFocal

    all_x = np.array([])
    all_y = np.array([])
    all_e1 = np.array([])
    all_e2 = np.array([])
    all_s = np.array([])

    if w is None:
        w = np.ones(len(x))

    x_bins = np.linspace(0,2048,nwhisk+1)
    y_bins = np.linspace(0,4096,2*nwhisk+1)
    print('x_bins = ',x_bins)
    print('y_bins = ',y_bins)

    ccdnums = np.unique(ccd)
    print('ccdnums = ',ccdnums)
    for ccdnum in ccdnums:
        mask = np.where(ccd == ccdnum)[0]
        print('ccdnum = ',ccdnum,', nstar = ',len(mask))
        if mask.sum() < 100: continue
        if ccdnum in [31, 61]: continue

        x_index = np.digitize(x[mask], x_bins)
        y_index = np.digitize(y[mask], y_bins)

        nbins = (len(x_bins)-1) * (len(y_bins)-1)
        bin_e1 = np.empty(nbins)
        bin_e2 = np.empty(nbins)
        bin_s = np.empty(nbins)
        bin_x = np.empty(nbins)
        bin_y = np.empty(nbins)
        bin_nz = np.empty(nbins, dtype=bool)
        sumesq = 0.

        for i in range(1, len(x_bins)):
            for j in range(1, len(y_bins)):
                k = (i-1)*(len(y_bins)-1) + (j-1)
                mask2 = np.where( (x_index == i) & (y_index == j) )[0]
                ww = w[mask][mask2]
                ws = ww.sum()
                bin_e1[k] = (ww * e1[mask][mask2]).sum() / ws
                bin_e2[k] = (ww * e2[mask][mask2]).sum() / ws
                bin_s[k] = (ww * s[mask][mask2]).sum() / ws
                bin_x[k] = (ww * x[mask][mask2]).sum() / ws
                bin_y[k] = (ww * y[mask][mask2]).sum() / ws
                bin_nz[k] = len(mask2) > 0
                print(i,j,k,len(mask2), bin_e1[k], bin_e2[k])
                sumesq += bin_e1[k]**2 + bin_e2[k]**2

        print('rms e = ',np.sqrt(sumesq / nbins))
        print('ccdnum = ',ccdnum)
        print('bin_x,y = ',bin_x,bin_y)
        focal_x, focal_y = toFocal(int(ccdnum), bin_x, bin_y)
        print('x,y = ',focal_x, focal_y)

        print('num with count > 0 = ',bin_nz.sum())
        all_x = np.append(all_x, focal_x[bin_nz])
        all_y = np.append(all_y, focal_y[bin_nz])
        all_e1 = np.append(all_e1, bin_e1[bin_nz])
        all_e2 = np.append(all_e2, bin_e2[bin_nz])
        all_s = np.append(all_s, bin_s[bin_nz])

    return all_x, all_y, all_e1, all_e2, all_s


def make_psf_whiskers(x, y, e1, e2, T, de1, de2, dT):

    import matplotlib
    matplotlib.use('Agg') # needs to be done before import pyplot
    import matplotlib.pyplot as plt

    # I don't know why it isn't finding this correctly.  It's in the right directory according to
    # http://matplotlib.org/users/style_sheets.html
    # But using the explicit path works.
    plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
    fig, ax = plt.subplots(1, 2, sharey=True, subplot_kw={'aspect' : 'equal'})
    print('fig = ',fig)
    print('ax = ',ax)

    theta = np.arctan2(e2,e1)/2.
    r = np.sqrt(e1**2 + e2**2)
    u = r*np.cos(theta)
    v = r*np.sin(theta)
    ax[0].set_xlim(-250,250)
    ax[0].set_ylim(-250,250)
    qv = ax[0].quiver(x,y,u,v, pivot='middle', scale_units='xy',
                        headwidth=0., headlength=0., headaxislength=0.,
                        width=0.001, scale=2.5e-3, color='blue')
    print('qv = ',qv)
    ref_scale = 0.03
    ref_label = 'e = ' + str(ref_scale)
    ax[0].quiverkey(qv, 0.10, 0.10, ref_scale, ref_label,
                      coordinates='axes', color='darkred', labelcolor='darkred',
                      labelpos='S')# fontproperties={'size':'small'})
    ax[0].axis('off')
    print('Done ax[0]')

    theta = np.arctan2(de2,de1)/2.
    r = np.sqrt(de1**2 + de2**2)
    u = r*np.cos(theta)
    v = r*np.sin(theta)
    ax[1].set_xlim(-250,250)
    ax[1].set_ylim(-250,250)
    qv = ax[1].quiver(x,y,u,v, pivot='middle', scale_units='xy',
                        headwidth=0., headlength=0., headaxislength=0.,
                        width=0.001, scale=4.0e-4, color='blue')
    print('qv = ',qv)
    ref_scale = 0.03
    ref_label = 'de = ' + str(ref_scale)
    ax[1].quiverkey(qv, 0.90, 0.10, ref_scale, ref_label,
                      coordinates='axes', color='darkred', labelcolor='darkred',
                      labelpos='S')# fontproperties={'size':'small'})
    ax[1].axis('off')
    print('Done ax[1]')

    fig.set_size_inches(7.5,4.0)
    fig.tight_layout()
    plt.savefig('both_psf_whiskers.eps')

def make_whiskers(x, y, e1, e2, s, filename, scale=1, auto_size=False, title=None, ref=0.01,
                  ref_name='$e$', alt_ref=None):

    import matplotlib
    matplotlib.use('Agg') # needs to be done before import pyplot
    import matplotlib.pyplot as plt

    plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    print('Start make_whiskers')

    print('x = ',x)
    print('y = ',y)
    print('e1 = ',e1)
    print('e2 = ',e2)
    theta = np.arctan2(e2,e1)/2.
    r = np.sqrt(e1**2 + e2**2)
    print('theta = ',theta)
    print('r = ',r)
    u = r*np.cos(theta)
    v = r*np.sin(theta)
    if auto_size:
        ax.set_aspect('equal')
    else:
        ax.set_xlim(-1.1,1.1)
        ax.set_ylim(-1.1,1.1)
        #ax.invert_xaxis()
        #ax.invert_yaxis()
    print('u = ',u)
    print('v = ',v)
    # Units are currently mm.  Switch to degrees
    # 15e-3 mm/pixel / 0.26 arcsec/pixel * 3600 arcsec/degree
    pixel_scale = 15e-3 / 0.26 * 3600.
    # Note: bigger scale means smaller whiskers
    qv = ax.quiver(x/pixel_scale, y/pixel_scale, u, v, pivot='middle', scale_units='xy',
                   headwidth=0., headlength=0., headaxislength=0.,
                   width=0.001, scale=1.0e-3*scale*pixel_scale, color='blue')
    print('qv = ',qv)
    ref_label = (ref_name + '$ = %s$'%ref).replace('$$','')
    ax.quiverkey(qv, 0.15, 0.10, ref, ref_label,
                 coordinates='axes', color='darkred', labelcolor='darkred',
                 labelpos='S', fontproperties={'size':16})

    if alt_ref is not None:
        ref_label = (ref_name + '$ = %s$'%alt_ref).replace('$$','')
        ax.quiverkey(qv, 0.85, 0.10, alt_ref, ref_label,
                     coordinates='axes', color='darkred', labelcolor='darkred',
                     labelpos='S', fontproperties={'size':16})

    ax.axis('off')
    print('Done ax')

    fig.text(0.20, 0.85, title, fontsize=16)
    #ax.set_xlabel('Y Position (degrees)')
    #ax.set_ylabel('X Position (degrees)')

    #fig.set_size_inches(7.5,4.0)
    #fig.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)

    np.savetxt(os.path.splitext(filename)[0] + '.dat',
                  np.array(zip(x, y, u, v, e1, e2, s)), fmt='%r',
                  header='x  y  u (=e cos(theta/2))  v (=e sin(theta/2))  e1  e2  size')

def psf_whiskers(ccd, x, y, e1, e2, T, de1, de2, dT):
    psf_binned_data = bin_by_fov(ccd, x, y, e1, e2, T, nwhisk=4)
    make_whiskers(*psf_binned_data, filename='psf_whiskers.eps', scale=3, title='PSF',
                  ref=0.01, alt_ref=0.03)
    resid_binned_data = bin_by_fov(ccd, x, y, de1, de2, dT, nwhisk=4)
    make_whiskers(*resid_binned_data, filename='resid_whiskers.eps', scale=0.3, title='PSF residual',
                  ref=0.01, alt_ref=0.03, ref_name=r'$\delta e$')
    resid_binned_data = bin_by_fov(ccd, x, y, de1, de2, dT, nwhisk=4)
    make_whiskers(*resid_binned_data, filename='sm_resid_whiskers.eps', scale=3, title='PSF residual',
                  ref=0.01, alt_ref=0.03, ref_name=r'$\delta e$')
    #make_psf_whiskers(x,y,e1,e2,T,de1,de2,dT)

def gal_whiskers(ccd, x, y, e1, e2, s, w, filename, scale=1, title=None):
    binned_data = bin_by_fov(ccd, x, y, e1, e2, s, w=w, nwhisk=4*scale)
    make_whiskers(*binned_data, filename=filename, scale=0.5*scale, title=title)

def wmean(e, w, mask):
    import np
    return np.sum(e[mask] * w[mask]) / np.sum(w[mask])

def evscol(ccd, x, y, e1, e2, s, w, filename, title=None):
    import matplotlib.pyplot as plt
    import np

    # Bin data by column
    ncols_per_bin = 8
    skip = 20
    bins = np.linspace(skip, 2048-skip, (2048-2*skip)/ncols_per_bin + 1)
    print('bins = ',bins)

    # Find mean e1, e2 in each bin
    denom = np.histogram(x, bins, weights=w)[0]
    #index = np.digitize(x, bins)
    e1_bins = np.histogram(x, bins, weights=e1*w)[0] / denom
    #e1_bins = [np.mean(e1[index == i]) for i in range(1, len(bins))]
    print('e1_bins = ',e1_bins)
    e2_bins = np.histogram(x, bins, weights=e2*w)[0] / denom
    #e2_bins = [np.mean(e2[index == i]) for i in range(1, len(bins))]
    print('e2_bins = ',e2_bins)
    x_bins = np.histogram(x, bins, weights=x*w)[0] / denom
    #x_bins = [np.mean(x[index == i]) for i in range(1, len(bins))]
    print('x_bins = ',x_bins)

    # Make 3 axes, where the first one spans the whole row.
    # cf. http://matplotlib.org/users/gridspec.html
    plt.style.use('supermongo')
    fig = plt.figure()
    ax1 = plt.subplot2grid( (2,2), (0,0), colspan=2 )
    ax2 = plt.subplot2grid( (2,2), (1,0) )
    ax3 = plt.subplot2grid( (2,2), (1,1) )

    # Draw a line through the mean e1, e2 to help guide the eye
    sw = np.sum(w)
    mean_e1 = np.sum(e1*w)/sw
    mean_e2 = np.sum(e2*w)/sw
    ax1.plot( [0,2048], [mean_e1, mean_e1], color='red', ls='--', lw=0.5)
    ax2.plot( [0,2048], [mean_e1, mean_e1], color='red', ls='--', lw=0.5)
    ax3.plot( [0,2048], [mean_e1, mean_e1], color='red', ls='--', lw=0.5)
    ax1.plot( [0,2048], [mean_e2, mean_e2], color='blue', ls='--', lw=0.5)
    ax2.plot( [0,2048], [mean_e2, mean_e2], color='blue', ls='--', lw=0.5)
    ax3.plot( [0,2048], [mean_e2, mean_e2], color='blue', ls='--', lw=0.5)

    # Plot data on first three axes:
    ax1.scatter(x_bins, e1_bins, marker='o', color='red', s=1.5, label=r'$\langle e_1 \rangle$')
    ax1.scatter(x_bins, e2_bins, marker='o', color='blue', s=1.5, label=r'$\langle e_2 \rangle$')

    ax2.scatter(x_bins, e1_bins, marker='o', color='red', s=3.0)
    ax2.scatter(x_bins, e2_bins, marker='o', color='blue', s=3.0)

    ax3.scatter(x_bins, e1_bins, marker='o', color='red', s=3.0)
    ax3.scatter(x_bins, e2_bins, marker='o', color='blue', s=3.0)

    # Draw the legend only once
    ax1.legend(loc=(0.05, 0.04), fontsize=10, frameon=True)#, fancybox=True)

    # Limit the axis ranges
    ax1.set_xlim(0, 2048)    # Show the whole range, but no points in first or last 20.
    ax2.set_xlim(0, 75)
    ax3.set_xlim(1973,2048)

    ymin = -0.004
    ymax = 0.003
    ax1.set_ylim(ymin, ymax)
    ax2.set_ylim(ymin, ymax)
    ax3.set_ylim(ymin, ymax)

    # Make ax2, ax3 look like a single plot with a broken x-axis
    ax2.spines['right'].set_visible(False) # Hide the right spine
    ax2.yaxis.tick_left()                  # Only put ticks on the left.

    ax3.spines['left'].set_visible(False)  # Hide the left spine
    ax3.yaxis.tick_right()                 # Only put ticks on the right.
    ax3.yaxis.set_ticklabels([])           # Don't label the y ticks.

    # It wants to label ever 0.001, but I think it looks nicer every 0.002.
    ax1.yaxis.set_ticks( np.arange(-0.004, 0.004, 0.002) )
    ax2.yaxis.set_ticks( np.arange(-0.004, 0.004, 0.002) )
    ax3.yaxis.set_ticks( np.arange(-0.004, 0.004, 0.002) )

    ax2.xaxis.set_ticks( np.arange(0, 80, 20) )
    ax3.xaxis.set_ticks( np.arange(1980, 2060, 20) )

    # Make little diagonal cuts to make the discontinuity clearer:
    # cf. http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis
    d = .03 # how big to make the diagonal lines in axes coordinates
    m = 0.6  # slope of the cut lines
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    ax2.plot((1-m*d,1+m*d),(-d,+d), **kwargs) # top-left diagonal
    ax2.plot((1-m*d,1+m*d),(1-d,1+d), **kwargs) # bottom-left diagonal

    kwargs.update(transform=ax3.transAxes) # switch to the bottom axes
    ax3.plot((-m*d,m*d),(-d,+d), **kwargs) # top-right diagonal
    ax3.plot((-m*d,m*d),(1-d,1+d), **kwargs) # bottom-right diagonal

    # Squeeze them a bit closer together
    plt.subplots_adjust(wspace=0.06)

    # Put the y label on ax1, ax2
    ax1.set_ylabel(r'$\langle e \rangle$')
    ax2.set_ylabel(r'$\langle e \rangle$')
    # Make a little more room for the label.
    plt.subplots_adjust(left=0.16)

    # For the x axis, we want it to be centered.  Easiest to just do this by hand.
    fig.text(0.5, 0.04, 'X Position on Chip', ha='center', va='center')
    # Make room.
    plt.subplots_adjust(bottom=0.12)

    # Shade in a grey region where the chips are masked
    ax1.fill_between([0,15],[ymin,ymin], [ymax,ymax], color='LightGrey')
    ax1.fill_between([2033,2048],[ymin,ymin], [ymax,ymax], color='LightGrey')
    ax2.fill_between([0,15],[ymin,ymin], [ymax,ymax], color='LightGrey')
    ax3.fill_between([2033,2048],[ymin,ymin], [ymax,ymax], color='LightGrey')

    if title is not None:
        fig.text(0.77, 0.57, title, fontsize=16)

    fig.tight_layout()
    plt.savefig(filename)


def main():

    args = parse_args()
    print('args = ',args)

    # Make the work directory if it does not exist yet.
    work = os.path.expanduser(args.work)
    print('work dir = ',work)
    try:
        if not os.path.isdir(work):
            os.makedirs(work)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    if args.use_psfex:
        prefix='psfex'
    else:
        prefix='piff'


    if True:
        if args.file != '':
            print('Read file ',args.file)
            with open(args.file) as fin:
                exps = [ line.strip() for line in fin if line[0] != '#' ]
            print('File included %d exposures'%len(exps))
        else:
            exps = args.exps
            print('Explicit listing of %d exposures'%len(exps))
        exps = sorted(exps)

        keys = ['ra', 'dec', 'x', 'y', 'mag', 'obs_e1', 'obs_e2', 'obs_T',
                prefix+'_e1', prefix+'_e2', prefix+'_T']
        data, bands, tilings = read_data(exps, work, keys, limit_bands=args.bands, prefix=prefix,
                                         use_reserved=args.use_reserved, frac=args.frac)
        e1 = data['obs_e1']
        e2 = data['obs_e2']
        T = data['obs_T']
        p_e1 = data[prefix+'_e1']
        p_e2 = data[prefix+'_e2']
        p_T = data[prefix+'_T']

        de1 = e1-p_e1
        de2 = e2-p_e2
        dT = (T-p_T)/T

        #psf_resid(data['mag'], de1, de2, dT)

        psf_whiskers(data['ccd'], data['x'], data['y'], e1, e2, T, de1, de2, dT)

    if False:
        ngmix_data = get_ngmix_epoch_data()
        gal_whiskers(*ngmix_data, filename='ngmix_whiskers.eps', title='ngmix')

        ngmix_data2 = get_ngmix_epoch_data(use_gold=False)
        evscol(*ngmix_data2, filename='ngmix_evscol.eps', title='ngmix')

        ccd31 = ccd==31
        ccd, x, y, e1, e2, s, w = ngmix_data
        gal_whiskers(ccd[ccd31], x[ccd31], y[ccd31], e1[ccd31], e2[ccd31], s[ccd31], w[ccd31],
                     filename='ngmix_ccd31.eps', scale=5, auto_size=True, title='ccd31')

    if False:
        im3shape_data = get_im3shape_epoch_data()
        gal_whiskers(*im3shape_data, filename='im3shape_whiskers.eps', title='im3shape')

        im3shape_data2 = get_im3shape_epoch_data(use_gold=False)
        evscol(*im3shape_data2, filename='im3shape_evscol.eps', title='im3shape')

if __name__ == "__main__":
    main()
