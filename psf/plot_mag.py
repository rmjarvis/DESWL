#! /usr/bin/env python

from __future__ import print_function
import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os
import glob
import galsim
import json
import numpy as np
from read_psf_cats import read_data, band_combinations

plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')

RESERVED = 64

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--file', default='',
                        help='list of exposures (in lieu of exps)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Options
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--use_psfex', default=False, action='store_const', const=True,
                        help='Use PSFEx rather than Piff model')
    parser.add_argument('--bands', default='grizY', type=str,
                        help='Limit to the given bands')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')

    args = parser.parse_args()
    return args


def break_axes(ax, ax2):
    """Do what is needed to axes to make them have a broken x axis and share the y axis
    """
    # cf. http://matplotlib.org/examples/pylab_examples/broken_axis.html
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')

    d = 0.015
    kwargs = dict(transform=ax.transAxes, color='black', clip_on=False)
    ax.plot((1-d,1+d),(-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d,+d),(-d,+d), **kwargs)
    ax2.plot((-d,+d),(1-d,1+d), **kwargs)


def add_to_list(band, vlist, value):
    vlist['griz'].append(value)
    if band not in vlist.keys():
        vlist[band] = []
    if 'g' not in band:
        vlist['riz'].append(value)
    if 'g' not in band and 'z' not in band:
        vlist['ri'].append(value)
    vlist[band].append(value)


def bin_by_mag(m, dT, de1, de2, min_mused, bands):

    import numpy as np
    import matplotlib.pyplot as plt

    min_mag = 13.5
    max_mag = 21
    #min_mused = 15

    # Bin by mag
    mag_bins = np.linspace(min_mag,max_mag,71)
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

    fig, axes = plt.subplots(2,1, sharex=True)

    ax = axes[0]
    ax.set_ylim(-0.002,0.012)
    ax.plot([min_mag,max_mag], [0,0], color='black')
    ax.plot([min_mused,min_mused],[-1,1], color='Grey')
    ax.fill( [min_mag,min_mag,min_mused,min_mused], [-1,1,1,-1], fill=False, hatch='/', color='Grey')
    t_line = ax.errorbar(mag_bins[:-1], bin_dT, yerr=bin_dT_err, color='green', fmt='o')
    ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r'$T_{\rm PSF} - T_{\rm model} ({\rm arcsec}^2)$')

    ax = axes[1]
    ax.set_ylim(-3.e-4,6.5e-4)
    ax.plot([min_mag,max_mag], [0,0], color='black')
    ax.plot([min_mused,min_mused],[-1,1], color='Grey')
    ax.fill( [min_mag,min_mag,min_mused,min_mused], [-1,1,1,-1], fill=False, hatch='/', color='Grey')
    e1_line = ax.errorbar(mag_bins[:-1], bin_de1, yerr=bin_de1_err, color='red', fmt='o')
    e2_line = ax.errorbar(mag_bins[:-1], bin_de2, yerr=bin_de2_err, color='blue', fmt='o')
    ax.legend([e1_line, e2_line], [r'$\delta e_1$', r'$\delta e_2$'])
    ax.set_ylabel(r'$e_{\rm PSF} - e_{\rm model}$')

    ax.set_xlim(min_mag,max_mag)
    ax.set_xlabel('Magnitude')

    fig.set_size_inches(7.0,10.0)
    plt.tight_layout()
    plt.savefig('dpsf_mag_' + bands + '.pdf')

    if True:
        cols = np.array((mag_bins[:-1],
                         bin_dT, bin_dT_err,
                         bin_de1, bin_de1_err,
                         bin_de2, bin_de2_err))
        outfile = 'dpsf_mag_' + bands + '.dat'
        np.savetxt(outfile, cols, fmt='%.6e',
                   header='mag_bins  '+
                          'dT  sig_dT '+
                          'de1  sig_de1 '+
                          'de2  sig_de2 ')
        print('wrote',outfile)


def bin_by_chip_pos(x, dT, de1, de2, bands, xy):

    import numpy as np
    import matplotlib.pyplot as plt

    # Bin by chip x or y position
    if xy == 'x':
        xmax = 2048
    else:
        xmax = 4096

    x_bins = np.linspace(0,xmax,129)
    print('x_bins = ',x_bins)
    index = np.digitize(x, x_bins)
    print('len(index) = ',len(index))
    bin_de1 = [de1[index == i].mean() for i in range(1, len(x_bins))]
    print('bin_de1 = ',bin_de1)
    bin_de2 = [de2[index == i].mean() for i in range(1, len(x_bins))]
    print('bin_de2 = ',bin_de2)
    bin_dT = [dT[index == i].mean() for i in range(1, len(x_bins))]
    print('bin_dT = ',bin_dT)
    bin_de1_err = [ np.sqrt(de1[index == i].var() / len(de1[index == i]))
                    for i in range(1, len(x_bins)) ]
    print('bin_de1_err = ',bin_de1_err)
    bin_de2_err = [ np.sqrt(de2[index == i].var() / len(de2[index == i]))
                    for i in range(1, len(x_bins)) ]
    print('bin_de2_err = ',bin_de2_err)
    bin_dT_err = [ np.sqrt(dT[index == i].var() / len(dT[index == i]))
                    for i in range(1, len(x_bins)) ]
    print('bin_dT_err = ',bin_dT_err)

    # Fix up nans
    for i in range(1,len(x_bins)):
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
    plt.title('PSF Size residuals')
    plt.xlim(0,xmax)
    plt.ylim(0,0.0025)
    plt.plot([0,xmax], [0,0], color='black')
    plt.errorbar(x_bins[2:-3], bin_dT[2:-2], yerr=bin_dT_err[2:-2], color='blue', fmt='o')
    plt.xlabel('Chip '+xy+' position')
    plt.ylabel('$size_{psf} - size_{model}$')
    plt.savefig('dsize_'+xy+'_' + bands + '.png')
    plt.savefig('dsize_'+xy+'_' + bands + '.pdf')
    plt.savefig('dsize_'+xy+'_' + bands + '.eps')

    plt.clf()
    plt.title('PSF Ellipticity residuals')
    plt.xlim(0,xmax)
    plt.ylim(-0.002,0.002)
    plt.plot([0,xmax], [0,0], color='black')
    e1_line = plt.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
    e2_line = plt.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
    plt.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
    plt.xlabel('Chip '+xy+' position')
    plt.ylabel('$e_{psf} - e_{model}$')
    plt.savefig('de_'+xy+'_' + bands + '.png')
    plt.savefig('de_'+xy+'_' + bands + '.pdf')
    plt.savefig('de_'+xy+'_' + bands + '.eps')

    # Make broken x-axis.
    plt.clf()
    f, (ax,ax2) = plt.subplots(1,2,sharey=True)
    ax.set_xlim(0,200)
    ax.set_ylim(0,0.0025)
    ax.plot([0,xmax], [0,0], color='black')
    ax.errorbar(x_bins[2:-3], bin_dT[2:-2], yerr=bin_dT_err[2:-2], color='blue', fmt='o')

    ax2.set_xlim(xmax-200,xmax)
    ax2.set_ylim(0,0.0025)
    ax2.plot([0,xmax], [0,0], color='black')
    ax2.errorbar(x_bins[2:-3], bin_dT[2:-2], yerr=bin_dT_err[2:-2], color='blue', fmt='o')

    break_axes(ax,ax2)

    ax.set_title('PSF Size residuals', x=1.08, y=1.03)
    ax.set_xlabel('Chip '+xy+' position')
    ax.xaxis.set_label_coords(1.08,-0.05)
    ax.set_ylabel('$size_{psf} - size_{model}$')

    plt.savefig('dsize_'+xy+'2_' + bands + '.png')
    plt.savefig('dsize_'+xy+'2_' + bands + '.pdf')
    plt.savefig('dsize_'+xy+'2_' + bands + '.eps')

    plt.clf()
    f, (ax,ax2) = plt.subplots(1,2,sharey=True)
    ax.set_xlim(0,200)
    ax.set_ylim(-0.002,0.002)
    ax.plot([0,xmax], [0,0], color='black')
    ax.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
    ax.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')

    ax2.set_xlim(xmax-200,xmax)
    ax2.set_ylim(-0.002,0.002)
    ax2.plot([0,xmax], [0,0], color='black')
    e1_line = ax2.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
    e2_line = ax2.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
    ax2.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])

    break_axes(ax,ax2)

    ax.set_title('PSF Ellipticity residuals', x=1.08, y=1.03)
    ax.set_xlabel('Chip '+xy+' position')
    ax.xaxis.set_label_coords(1.08,-0.05)
    ax.set_ylabel('$e_{psf} - e_{model}$')

    plt.savefig('de_'+xy+'2_' + bands + '.png')
    plt.savefig('de_'+xy+'2_' + bands + '.pdf')
    plt.savefig('de_'+xy+'2_' + bands + '.eps')


def bin_by_fov(ccd, x, y, dT, de1, de2, bands):

    import numpy as np
    import matplotlib.pyplot as plt
    from toFocal import toFocal

    all_x = np.array([])
    all_y = np.array([])
    all_e1 = np.array([])
    all_e2 = np.array([])
    all_T = np.array([])

    nwhisk = 5
    x_bins = np.linspace(0,2048,nwhisk+1)
    y_bins = np.linspace(0,4096,2*nwhisk+1)

    ccdnums = np.unique(ccd)
    for ccdnum in ccdnums:
        mask = np.where(ccd == ccdnum)[0]
        print('ccdnum = ',ccdnum,', nstar = ',mask.sum())
        if mask.sum() < 100: continue

        x_index = np.digitize(x[mask], x_bins)
        y_index = np.digitize(y[mask], y_bins)

        bin_de1 = np.array([ de1[mask][(x_index == i) & (y_index == j)].mean()
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_de1 = ',bin_de1)
        bin_de2 = np.array([ de2[mask][(x_index == i) & (y_index == j)].mean()
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_de2 = ',bin_de2)
        bin_dT = np.array([ dT[mask][(x_index == i) & (y_index == j)].mean()
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_dT = ',bin_dT)
        bin_x = np.array([ x[mask][(x_index == i) & (y_index == j)].mean()
                  for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_x = ',bin_x)
        bin_y = np.array([ y[mask][(x_index == i) & (y_index == j)].mean()
                  for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_y = ',bin_y)
        bin_count = np.array([ ((x_index == i) & (y_index == j)).sum()
                      for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print('bin_count = ',bin_count)

        focal_x, focal_y = toFocal(ccdnum, bin_x, bin_y)

        mask2 = np.where(bin_count > 0)[0]
        print('num with count > 0 = ',mask2.sum())
        all_x = np.append(all_x, focal_x[mask2])
        all_y = np.append(all_y, focal_y[mask2])
        all_e1 = np.append(all_e1, bin_de1[mask2])
        all_e2 = np.append(all_e2, bin_de2[mask2])
        all_T = np.append(all_T, bin_dT[mask2])


    plt.clf()
    #plt.title('PSF Ellipticity residuals in DES focal plane')
    theta = np.arctan2(all_e2,all_e1)/2.
    r = np.sqrt(all_e1**2 + all_e2**2)
    u = r*np.cos(theta)
    v = r*np.sin(theta)
    plt.xlim(-250,250)
    plt.ylim(-250,250)
    print('all_x = ',all_x)
    print('len(all_x) = ',len(all_x))
    qv = plt.quiver(all_x,all_y,u,v, pivot='middle', scale_units='xy',
                    headwidth=0., headlength=0., headaxislength=0.,
                    width=0.001, scale=1.e-3, color='blue')
    ref_scale = 0.01
    if 'raw' in bands:
        ref_label = 'e = ' + str(ref_scale)
    else:
        ref_label = 'de = ' + str(ref_scale)
    plt.quiverkey(qv, 0.10, 0.08, ref_scale, ref_label,
                  coordinates='axes', color='darkred', labelcolor='darkred',
                  labelpos='E', fontproperties={'size':'x-small'})
    plt.axis('off')
    plt.savefig('de_fov_' + bands + '.png')
    plt.savefig('de_fov_' + bands + '.pdf')
    plt.savefig('de_fov_' + bands + '.eps')

def make_hist(dT, T, de1, de2, bands):

    import numpy as np
    import matplotlib.pyplot as plt

    nbins = 1000
    range = (-0.1,0.1)

    fig, ax = plt.subplots(1,4, sharey=True)

    ax[0].hist(dT/T, bins=nbins, histtype='step', range=range, fill=True)
    ax[0].set_xlabel('dT/T')
    ax[1].hist(de1, bins=nbins, histtype='step', range=range, fill=True)
    ax[1].set_xlabel('de1')
    ax[2].hist(de2, bins=nbins, histtype='step', range=range, fill=True)
    ax[2].set_xlabel('de2')

    fig.set_size_inches(15.0,7.0)
    plt.tight_layout()
    plt.savefig('dpsf_hist_' + bands + '.pdf')

def main():
    matplotlib.use('Agg') # needs to be done before import pyplot

    args = parse_args()

    work = os.path.expanduser(args.work)
    print('work dir = ',work)

    if args.file != '':
        print('Read file ',args.file)
        with open(args.file) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
    else:
        exps = args.exps
    exps = sorted(exps)

    if args.use_psfex:
        prefix = 'psfex'
    else:
        prefix = 'piff'

    keys = ['ra', 'dec', 'x', 'y', 'mag', 'obs_e1', 'obs_e2', 'obs_T', 'obs_flag',
            prefix+'_e1', prefix+'_e2', prefix+'_T', prefix+'_flag']

    data, bands, tilings = read_data(exps, work, keys, limit_bands=args.bands, prefix=prefix,
                                     use_reserved=args.use_reserved, frac=args.frac)

    use_bands = band_combinations(args.bands)
    for bands in use_bands:
        this_data = data[np.in1d(data['band'], list(bands))]
        if len(this_data) == 0:
            print('No files with bands ',bands)
            continue

        used = this_data[prefix+'_flag'] & ~RESERVED == 0

        #airmass = this_data['airmass']
        #sky = this_data['sky']
        #sigsky = this_data['sigsky']
        #fwhm = this_data['fwhm']

        #med_airmass = np.median(airmass)
        #med_sky = np.median(sky)
        #med_sigsky = np.median(sigsky)
        #med_fwhm = np.median(fwhm)
        #print('airmass: ',min(airmass),med_airmass,max(airmass))
        #print('sky: ',min(sky),med_sky,max(sky))
        #print('sigsky: ',min(sigsky),med_sigsky,max(sigsky))
        #print('fwhm: ',min(fwhm),med_fwhm,max(fwhm))

        ra = this_data['ra']
        dec = this_data['dec']
        x = this_data['x']
        y = this_data['y']
        m = this_data['mag']
        print('full mag range = ',np.min(m),np.max(m))
        print('used mag range = ',np.min(m[used]),np.max(m[used]))

        e1 = this_data['obs_e1']
        print('mean e1 = ',np.mean(e1))
        e2 = this_data['obs_e2']
        print('mean e2 = ',np.mean(e2))
        T = this_data['obs_T']
        print('mean s = ',np.mean(T))
        pe1 = this_data[prefix+'_e1']
        print('mean pe1 = ',np.mean(pe1))
        pe2 = this_data[prefix+'_e2']
        print('mean pe2 = ',np.mean(pe2))
        pT = this_data[prefix+'_T']
        print('mean pT = ',np.mean(pT))

        print('min mag = ',np.min(m))
        print('max mag = ',np.max(m))
        print('mean T (used) = ',np.mean(T[used]))
        print('mean e1 (used) = ',np.mean(e1[used]))
        print('mean e2 (used) = ',np.mean(e2[used]))
        print('mean pT (used) = ',np.mean(pT[used]))
        print('mean pe1 (used) = ',np.mean(pe1[used]))
        print('mean pe2 (used) = ',np.mean(pe2[used]))
        de1 = e1 - pe1
        de2 = e2 - pe2
        dT = T - pT
        print('mean dT (used) = ',np.mean(dT[used]))
        print('mean de1 (used) = ',np.mean(de1[used]))
        print('mean de2 (used) = ',np.mean(de2[used]))

        if args.use_psfex:
            min_mused = 0
        else:
            min_mused = np.min(m[used])
        print('min_mused = ',min_mused)
        bin_by_mag(m, dT, de1, de2, min_mused, bands)
        bin_by_mag(m[used], dT[used], de1[used], de2[used], min_mused, bands+'_used')

        make_hist(dT, T, de1, de2, bands)
        make_hist(dT[used], T[used], de1[used], de2[used], bands+'_used')

        #bin_by_mag(m, T, e1, e2, min_mused, bands+'_orig')
        #bin_by_mag(m, pT, pe1, pe2, min_mused, bands+'_model')

        #mask = sky<med_sky
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_lowsky')
        #mask = sky>med_sky
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_highsky')

        #mask = sigsky<med_sigsky
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_lowsigsky')
        #mask = sigsky>med_sigsky
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_highsigsky')

        #mask = mask & (airmass<med_airmass)
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_lowairmass')
        #mask = mask & (airmass>med_airmass)
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_highairmass')

        #mask = mask & (fwhm<med_fwhm)
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_lowfwhm')
        #mask = mask & (fwhm>med_fwhm)
        #bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, bands+'_highfwhm')

        #bin_by_mag(m, T-np.mean(T), e1-np.mean(e1), e2-np.mean(e2), min_mused, bands)

        #bin_by_chip_pos(x[used], dT[used], de1[used], de2[used], bands, 'x')
        #bin_by_chip_pos(y[used], dT[used], de1[used], de2[used], bands, 'y')

        #bin_by_fov(ccd[used], x[used], y[used], dT[used], de1[used], de2[used], bands)
        #bin_by_fov(ccd[used], x[used], y[used], T[used], e1[used], e2[used], bands + 'raw')


if __name__ == "__main__":
    main()
