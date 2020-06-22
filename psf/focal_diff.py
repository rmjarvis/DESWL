#! /usr/bin/env python
from __future__ import print_function
import numpy as np
import os
from read_psf_cats import read_data
from toFocal import toFocal, toFocalArcmin

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


def bin_by_fov(ccd, x, y, e1, e2, s, w=None, nwhisk=5):

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
    make_whiskers(*psf_binned_data, filename='psf_whiskers.pdf', scale=3, title='PSF',
                  ref=0.01, alt_ref=0.03)
    resid_binned_data = bin_by_fov(ccd, x, y, de1, de2, dT, nwhisk=4)
    make_whiskers(*resid_binned_data, filename='resid_whiskers.pdf', scale=0.3, title='PSF residual',
                  ref=0.01, alt_ref=0.03, ref_name=r'$\delta e$')
    resid_binned_data = bin_by_fov(ccd, x, y, de1, de2, dT, nwhisk=4)
    make_whiskers(*resid_binned_data, filename='sm_resid_whiskers.pdf', scale=3, title='PSF residual',
                  ref=0.01, alt_ref=0.03, ref_name=r'$\delta e$')

def psf_hex(ccd, x, y, e1, e2, T, de1, de2, dT):
    import matplotlib
    matplotlib.use('Agg') # needs to be done before import pyplot
    import matplotlib.pyplot as plt

    plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    u = np.zeros_like(x)
    v = np.zeros_like(y)
    for ccdnum in np.unique(ccd):
        uu, vv = toFocalArcmin(int(ccdnum), x[ccd==ccdnum], y[ccd==ccdnum])
        u[ccd==ccdnum] = uu
        v[ccd==ccdnum] = vv

    ngrid = 300
    vmax = 0.01
    fig, axs = plt.subplots(ncols=3, sharey=True, figsize=(15, 4))
    fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
    ax = axs[0]
    hb = ax.hexbin(u, v, de1, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\delta e_1$')
    ax.set_aspect('equal')
    ax.set_xlim(-75, 75)
    ax.set_ylim(-75, 75)
    #cb = fig.colorbar(hb, ax=ax)
    #cb.set_label('counts')

    ax = axs[1]
    hb = ax.hexbin(u, v, de2, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\delta e_2$')
    ax.set_aspect('equal')
    ax.set_xlim(-75, 75)
    ax.set_ylim(-75, 75)
    #cb = fig.colorbar(hb, ax=ax)

    ax = axs[2]
    hb = ax.hexbin(u, v, dT, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\delta T/T$')
    ax.set_aspect('equal')
    ax.set_xlim(-75, 75)
    ax.set_ylim(-75, 75)
    cb = fig.colorbar(hb, ax=ax)

    filename = 'psf_resid_hex.pdf'
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)


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

    psf_whiskers(data['ccd'], data['x'], data['y'], e1, e2, T, de1, de2, dT)
    psf_hex(data['ccd'], data['x'], data['y'], e1, e2, T, de1, de2, dT)

if __name__ == "__main__":
    main()
