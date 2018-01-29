#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')

    # Options
    parser.add_argument('--use_psfex', default=False, action='store_const', const=True,
                        help='Use PSFEx rather than Piff model')
    parser.add_argument('--use_reserved', default=False, action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')


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

def get_data(exps, work,
             mask_list, used_list, ccd_list, 
             airmass_list, sky_list, sigsky_list, fwhm_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, T_list, pe1_list, pe2_list, pT_list,
             tag, prefix='piff', reserved=False):

    import fitsio
    import numpy as np
    import os

    BAD_CCDS = [2, 31, 61]

    for exp in exps:

        print('Start work on exp = ',exp)
        expnum = int(exp)
        print('expnum = ',expnum)
        expinfo = fitsio.read(os.path.join(work, exp, 'exp_info_%d.fits'%expnum))

        if expnum not in expinfo['expnum']:
            print('expnum is not in expinfo!')
            print('expinfo[expnum] = ',expinfo['expnum'])
            #raise RuntimeError("Could not find information about this expnum")
            print('Skipping ',run,exp)
            continue

        i = np.nonzero(expinfo['expnum'] == expnum)[0][0]
        #print('i = ',i)
        band = expinfo['band'][i]
        print('band = ',band)
        
        for k in range(len(expinfo)):
            ccdnum = expinfo[k]['ccdnum']
            if expinfo[k]['flag'] != 0:
                print('Skipping ccd %d because it is blacklisted: '%ccdnum, expinfo[k]['flag'])
                continue
            if ccdnum in BAD_CCDS:
                print('Skipping ccd %d because it is BAD'%ccdnum)
                continue

            cat_file = os.path.join(work, exp, "psf_cat_%d_%d.fits"%(expnum,ccdnum))
            #print('cat_file = ',cat_file)
            try:
                data = fitsio.read(cat_file)
                flag = data[prefix+'_flag']
            except (OSError, IOError):
                print('Unable to open cat_file %s.  Skipping this file.'%cat_file)
                continue

            #print('max flag = ',max(flag))
            #print('min flag = ',min(flag))
            if len(flag) == 0:
                print('No data')
                continue
            if min(flag) < 0:
                print('!!!Flag is negative!!!')
            if reserved:
                mask = (flag == 64) | (flag == 65)
            else:
                mask = (flag == 0) | (flag == 1)

            if mask.sum() == 0:
                print('All objects in this exposure are flagged.')
                print('Probably due to astrometry flags. Skip this exposure.')
                continue

            used = (flag == 0)

            #print('nobject = ',sum(mask))
            #print('nused = ',sum(used))
            #print('full mag range = ',min(data['mag']),max(data['mag']))
            #print('masked mag range = ',min(data['mag'][mask]),max(data['mag'][mask]))
            #print('used mag range = ',min(data['mag'][used]),max(data['mag'][used]))
            #print('alt used mag range = ',min(data[used]['mag']),max(data[used]['mag']))

            T = data['obs_T']
            dT = (data[prefix + '_T'] - data['obs_T'])
            de1 = (data[prefix + '_e1'] - data['obs_e1'])
            de2 = (data[prefix + '_e2'] - data['obs_e2'])
            print(expnum, ccdnum, len(dT), band)
            print('dT = ',np.mean(dT[used]),np.std(dT[used]))
            print('de1 = ',np.mean(de1[used]),np.std(de1[used]))
            print('de2 = ',np.mean(de2[used]),np.std(de2[used]))
            if np.std(dT[used]/T[used]) > 0.03:
                print('high std(dT)')
                #continue
            if np.std(de1[used]) > 0.02:
                print('high std(de1)')
                #continue
            if np.std(de2[used]) > 0.02:
                print('high std(de2)')
                #continue
            if abs(np.mean(dT[used]/T[used])) > 0.03:
                print('high mean(dT)')
                #continue
            if abs(np.mean(de1[used])) > 0.02:
                print('high mean(de1)')
                #continue
            if abs(np.mean(de2[used])) > 0.02:
                print('high mean(de2)')
                #continue

            good = (abs(dT/data['obs_T']) < 0.1) & (abs(de1) < 0.1) & (abs(de2) < 0.1)

            add_to_list(band, mask_list, mask & good)
            add_to_list(band, used_list, used & good)
            add_to_list(band, ccd_list, [ccdnum]*len(used))

            add_to_list(band, ra_list, data['ra'])
            add_to_list(band, dec_list, data['dec'])
            add_to_list(band, x_list, data['x'])
            add_to_list(band, y_list, data['y'])
            add_to_list(band, m_list, data['mag'])

            add_to_list(band, e1_list, data['obs_e1'])
            add_to_list(band, e2_list, data['obs_e2'])
            add_to_list(band, T_list, data['obs_T'])

            add_to_list(band, pe1_list, data[prefix + '_e1'])
            add_to_list(band, pe2_list, data[prefix + '_e2'])
            add_to_list(band, pT_list, data[prefix + '_T'])

            n = len(mask)
            #add_to_list(band, airmass_list, [expinfo['airmass'][k]] * n)
            #add_to_list(band, sky_list, [expinfo['sky'][k]] * n)
            #add_to_list(band, sigsky_list, [expinfo['sigsky'][k]] * n)
            #add_to_list(band, fwhm_list, [expinfo['fwhm'][k]] * n)

    print('\nFinished processing all exposures')


def bin_by_mag(m, dT, de1, de2, min_mused, key):

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
    plt.savefig('dpsf_mag_' + key + '.pdf')

    if True:
        cols = numpy.array((mag_bins[:-1],
                            bin_dT, bin_dT_err,
                            bin_de1, bin_de1_err,
                            bin_de2, bin_de2_err))
        outfile = 'dpsf_mag_' + key + '.dat'
        numpy.savetxt(outfile, cols, fmt='%.6e',
                      header='mag_bins  '+
                             'dT  sig_dT '+
                             'de1  sig_de1 '+
                             'de2  sig_de2 ')
        print('wrote',outfile)
 

def bin_by_chip_pos(x, dT, de1, de2, key, xy):

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
    plt.savefig('dsize_'+xy+'_' + key + '.png')
    plt.savefig('dsize_'+xy+'_' + key + '.pdf')
    plt.savefig('dsize_'+xy+'_' + key + '.eps')

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
    plt.savefig('de_'+xy+'_' + key + '.png')
    plt.savefig('de_'+xy+'_' + key + '.pdf')
    plt.savefig('de_'+xy+'_' + key + '.eps')

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

    plt.savefig('dsize_'+xy+'2_' + key + '.png')
    plt.savefig('dsize_'+xy+'2_' + key + '.pdf')
    plt.savefig('dsize_'+xy+'2_' + key + '.eps')

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

    plt.savefig('de_'+xy+'2_' + key + '.png')
    plt.savefig('de_'+xy+'2_' + key + '.pdf')
    plt.savefig('de_'+xy+'2_' + key + '.eps')


def bin_by_fov(ccd, x, y, dT, de1, de2, key):

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
    if 'raw' in key:
        ref_label = 'e = ' + str(ref_scale)
    else:
        ref_label = 'de = ' + str(ref_scale)
    plt.quiverkey(qv, 0.10, 0.08, ref_scale, ref_label,
                  coordinates='axes', color='darkred', labelcolor='darkred',
                  labelpos='E', fontproperties={'size':'x-small'})
    plt.axis('off')
    plt.savefig('de_fov_' + key + '.png')
    plt.savefig('de_fov_' + key + '.pdf')
    plt.savefig('de_fov_' + key + '.eps')

def make_hist(dT, T, de1, de2, key):

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
    plt.savefig('dpsf_hist_' + key + '.pdf')

def main():
    import os
    import glob
    import galsim
    import json
    import numpy as np
    import matplotlib
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

    mask_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    used_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    ccd_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    airmass_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    sky_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    sigsky_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    fwhm_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    ra_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    dec_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    x_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    y_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    m_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    e1_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    e2_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    T_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    pe1_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    pe2_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    pT_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    if args.use_psfex:
        prefix = 'psfex'
    else:
        prefix = 'piff'

    get_data(exps, work,
             mask_list, used_list, ccd_list,
             airmass_list, sky_list, sigsky_list, fwhm_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, T_list, pe1_list, pe2_list, pT_list, args.tag,
             prefix=prefix, reserved=args.use_reserved)

    #for key in ra_list.keys():
    for key in ['r', 'i', 'z', 'riz']:
        if key not in mask_list: continue
        mask = np.concatenate(mask_list[key])
        used = np.concatenate(used_list[key])
        ccd = np.concatenate(ccd_list[key])

        #airmass = np.concatenate(airmass_list[key])
        #sky = np.concatenate(sky_list[key])
        #sigsky = np.concatenate(sigsky_list[key])
        #fwhm = np.concatenate(fwhm_list[key])

        #med_airmass = np.median(airmass)
        #med_sky = np.median(sky)
        #med_sigsky = np.median(sigsky)
        #med_fwhm = np.median(fwhm)
        #print('airmass: ',min(airmass),med_airmass,max(airmass))
        #print('sky: ',min(sky),med_sky,max(sky))
        #print('sigsky: ',min(sigsky),med_sigsky,max(sigsky))
        #print('fwhm: ',min(fwhm),med_fwhm,max(fwhm))

        ra = np.concatenate(ra_list[key])
        dec = np.concatenate(dec_list[key])
        x = np.concatenate(x_list[key])
        y = np.concatenate(y_list[key])
        m = np.concatenate(m_list[key])
        print('full mag range = ',np.min(m),np.max(m))
        print('masked mag range = ',np.min(m[mask]),np.max(m[mask]))
        print('used mag range = ',np.min(m[used]),np.max(m[used]))
        print('len(m) = ',len(m))
        print('len(mask) = ',len(mask))
        print('len(used) = ',len(used))
        print('sum(mask) = ',np.sum(mask))
        print('sum(used) = ',np.sum(used))

        e1 = np.concatenate(e1_list[key])
        print('mean e1 = ',np.mean(e1[mask]))
        e2 = np.concatenate(e2_list[key])
        print('mean e2 = ',np.mean(e2[mask]))
        T = np.concatenate(T_list[key])
        print('mean s = ',np.mean(T[mask]))
        pe1 = np.concatenate(pe1_list[key])
        print('mean pe1 = ',np.mean(pe1[mask]))
        pe2 = np.concatenate(pe2_list[key])
        print('mean pe2 = ',np.mean(pe2[mask]))
        pT = np.concatenate(pT_list[key])
        print('mean pT = ',np.mean(pT[mask]))

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
        bin_by_mag(m[mask], dT[mask], de1[mask], de2[mask], min_mused, key)
        bin_by_mag(m[used], dT[used], de1[used], de2[used], min_mused, key+'_used')

        make_hist(dT[mask], T[mask], de1[mask], de2[mask], key)
        make_hist(dT[used], T[used], de1[used], de2[used], key+'_used')

        #bin_by_mag(m[mask], T[mask], e1[mask], e2[mask], min_mused, key+'_orig')
        #bin_by_mag(m[mask], pT[mask], pe1[mask], pe2[mask], min_mused, key+'_model')

        #mask2 = mask & (sky<med_sky)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowsky')
        #mask2 = mask & (sky>med_sky)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_highsky')

        #mask2 = mask & (sigsky<med_sigsky)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowsigsky')
        #mask2 = mask & (sigsky>med_sigsky)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_highsigsky')

        #mask2 = mask & (airmass<med_airmass)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowairmass')
        #mask2 = mask & (airmass>med_airmass)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_highairmass')

        #mask2 = mask & (fwhm<med_fwhm)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowfwhm')
        #mask2 = mask & (fwhm>med_fwhm)
        #bin_by_mag(m[mask2], dT[mask2], de1[mask2], de2[mask2], min_mused, key+'_highfwhm')

        #bin_by_mag(m[mask], T[mask] - np.mean(T[mask]),
        #           e1[mask] - np.mean(e1[mask]), e2[mask] - np.mean(e2[mask]), min_mused, key)

        #bin_by_chip_pos(x[used], dT[used], de1[used], de2[used], key, 'x')
        #bin_by_chip_pos(y[used], dT[used], de1[used], de2[used], key, 'y')

        #bin_by_fov(ccd[used], x[used], y[used], dT[used], de1[used], de2[used], key)
        #bin_by_fov(ccd[used], x[used], y[used], T[used], e1[used], e2[used], key + 'raw')


if __name__ == "__main__":
    main()
