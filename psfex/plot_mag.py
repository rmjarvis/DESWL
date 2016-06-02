#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

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
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits.fz',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')

    # Options
    parser.add_argument('--single_ccd', default=False, action='store_const', const=True,
                        help='Only do 1 ccd per exposure (used for debugging)')

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

def add_to_list(filter, vlist, value):
    vlist['griz'].append(value)
    if filter not in vlist.keys():
        vlist[filter] = []
    if 'g' not in filter:
        vlist['riz'].append(value)
    if 'g' not in filter and 'z' not in filter:
        vlist['ri'].append(value)
    vlist[filter].append(value)

def get_data(runs, exps, work,
             mask_list, used_list, ccd_list, 
             airmass_list, sky_list, sigsky_list, fwhm_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, s_list, pe1_list, pe2_list, ps_list,
             tag, psfex = 'psfex'):

    import astropy.io.fits as pyfits
    import numpy
    import os

    expinfo_file = '/astro/u/mjarvis/work/exposure_info_' + tag + '.fits'
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    cat_dir = os.path.join(work,'psf_cats')

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        #print 'expnum = ',expnum

        if expnum not in expinfo['expnum']:
            print 'expnum is not in expinfo!'
            print 'expinfo[expnum] = ',expinfo['expnum']
            #raise RuntimeError("Could not find information about this expnum")
            print 'Skipping ',run,exp
            continue
        k = numpy.nonzero(expinfo['expnum'] == expnum)[0][0]
        #print 'k = ',k
        filter = expinfo['filter'][k]
        print 'filter = ',filter

        cat_file = os.path.join(cat_dir, exp + "_psf.fits")
        try:
            with pyfits.open(cat_file) as pyf:
                data = pyf[1].data.copy()
        except:
            print 'Could not open cat_file %s.  Skipping this one.'%cat_file
            continue
            

        #print 'max flag = ',max(data['flag'])
        #print 'min flag = ',min(data['flag'])
        if len(data['flag']) == 0:
            print 'No data'
            continue
        if min(data['flag']) < 0:
            print '!!!Flag is negative!!!'
        mask = (data['flag'] == 0) | (data['flag'] == 1)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        used = (data['flag'] == 0)
        bright = mask & (data['mag'] < 10.2)

        #print 'nobject = ',sum(mask)
        #print 'nused = ',sum(used)
        #print 'full mag range = ',min(data['mag']),max(data['mag'])
        #print 'masked mag range = ',min(data['mag'][mask]),max(data['mag'][mask])
        #print 'used mag range = ',min(data['mag'][used]),max(data['mag'][used])
        #print 'alt used mag range = ',min(data[used]['mag']),max(data[used]['mag'])

        add_to_list(filter, mask_list, mask)
        add_to_list(filter, used_list, used)
        add_to_list(filter, ccd_list, data['ccdnum'])

        add_to_list(filter, ra_list, data['ra'])
        add_to_list(filter, dec_list, data['dec'])
        add_to_list(filter, x_list, data['x'])
        add_to_list(filter, y_list, data['y'])
        add_to_list(filter, m_list, data['mag'])

        add_to_list(filter, e1_list, data['e1'])
        add_to_list(filter, e2_list, data['e2'])
        add_to_list(filter, s_list, data['size'])

        add_to_list(filter, pe1_list, data[psfex + '_e1'])
        add_to_list(filter, pe2_list, data[psfex + '_e2'])
        add_to_list(filter, ps_list, data[psfex + '_size'])

        n = len(mask)
        add_to_list(filter, airmass_list, [expinfo['airmass'][k]] * n)
        add_to_list(filter, sky_list, [expinfo['sky'][k]] * n)
        add_to_list(filter, sigsky_list, [expinfo['sigsky'][k]] * n)
        add_to_list(filter, fwhm_list, [expinfo['fwhm'][k]] * n)

    print '\nFinished processing all exposures'


def bin_by_mag(m, dt, de1, de2, min_mused, key):

    import numpy
    import matplotlib.pyplot as plt

    # Adjust for wrong zero point in final cut images.
    # From Eli:
    #   The typical r-band zeropoint for photometric data, top-of-the-atmosphere is 30.3. 
    #   The raw zeropoint in the finalcut catalogs is 25.0. 
    #   So just add 5.3 and you're in the ballpark.
    m = m + 5.3
    min_mused += 5.3

    # Bin by mag
    mag_bins = numpy.linspace(15.3,22,71)
    print 'mag_bins = ',mag_bins

    index = numpy.digitize(m, mag_bins)
    print 'len(index) = ',len(index)
    bin_de1 = [de1[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_de1 = ',bin_de1
    bin_de2 = [de2[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_de2 = ',bin_de2
    bin_dt = [dt[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_dt = ',bin_dt
    bin_de1_err = [ numpy.sqrt(de1[index == i].var() / len(de1[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_de1_err = ',bin_de1_err
    bin_de2_err = [ numpy.sqrt(de2[index == i].var() / len(de2[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_de2_err = ',bin_de2_err
    bin_dt_err = [ numpy.sqrt(dt[index == i].var() / len(dt[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_dt_err = ',bin_dt_err

    # Fix up nans
    for i in range(1,len(mag_bins)):
        if i not in index:
            bin_de1[i-1] = 0.
            bin_de2[i-1] = 0.
            bin_dt[i-1] = 0.
            bin_de1_err[i-1] = 0.
            bin_de2_err[i-1] = 0.
            bin_dt_err[i-1] = 0.
    print 'fixed nans'

    print 'index = ',index
    print 'bin_de1 = ',bin_de1
    print 'bin_de2 = ',bin_de2
    print 'bin_dt = ',bin_dt

    fig, axes = plt.subplots(2,1, sharex=True)

    ax = axes[0]
    ax.set_ylim(-0.002,0.012)
    ax.plot([15.3,21.5], [0,0], color='black')
    ax.plot([min_mused,min_mused],[-1,1], color='Grey')
    ax.fill( [15.3,15.3,min_mused,min_mused], [-1,1,1,-1], fill=False, hatch='/', color='Grey')
    t_line = ax.errorbar(mag_bins[:-1], bin_dt, yerr=bin_dt_err, color='green', fmt='o')
    ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r'$T_{\rm PSF} - T_{\rm model} ({\rm arcsec}^2)$')

    ax = axes[1]
    ax.set_ylim(-3.e-4,6.5e-4)
    ax.plot([15.3,21.5], [0,0], color='black')
    ax.plot([min_mused,min_mused],[-1,1], color='Grey')
    ax.fill( [15.3,15.3,min_mused,min_mused], [-1,1,1,-1], fill=False, hatch='/', color='Grey')
    e1_line = ax.errorbar(mag_bins[:-1], bin_de1, yerr=bin_de1_err, color='red', fmt='o')
    e2_line = ax.errorbar(mag_bins[:-1], bin_de2, yerr=bin_de2_err, color='blue', fmt='o')
    ax.legend([e1_line, e2_line], [r'$\delta e_1$', r'$\delta e_2$'])
    ax.set_ylabel(r'$e_{\rm PSF} - e_{\rm model}$')

    ax.set_xlim(15.3,21.5)
    ax.set_xlabel('Magnitude')

    fig.set_size_inches(7.0,10.0)
    plt.tight_layout()
    plt.savefig('dpsf_mag_' + key + '.pdf')


def bin_by_chip_pos(x, ds, de1, de2, key, xy):

    import numpy
    import matplotlib.pyplot as plt

    # Bin by chip x or y position
    if xy == 'x':
        xmax = 2048
    else:
        xmax = 4096

    x_bins = numpy.linspace(0,xmax,129)
    print 'x_bins = ',x_bins
    index = numpy.digitize(x, x_bins)
    print 'len(index) = ',len(index)
    bin_de1 = [de1[index == i].mean() for i in range(1, len(x_bins))]
    print 'bin_de1 = ',bin_de1
    bin_de2 = [de2[index == i].mean() for i in range(1, len(x_bins))]
    print 'bin_de2 = ',bin_de2
    bin_ds = [ds[index == i].mean() for i in range(1, len(x_bins))]
    print 'bin_ds = ',bin_ds
    bin_de1_err = [ numpy.sqrt(de1[index == i].var() / len(de1[index == i])) 
                    for i in range(1, len(x_bins)) ]
    print 'bin_de1_err = ',bin_de1_err
    bin_de2_err = [ numpy.sqrt(de2[index == i].var() / len(de2[index == i])) 
                    for i in range(1, len(x_bins)) ]
    print 'bin_de2_err = ',bin_de2_err
    bin_ds_err = [ numpy.sqrt(ds[index == i].var() / len(ds[index == i])) 
                    for i in range(1, len(x_bins)) ]
    print 'bin_ds_err = ',bin_ds_err

    # Fix up nans
    for i in range(1,len(x_bins)):
        if i not in index:
            bin_de1[i-1] = 0.
            bin_de2[i-1] = 0.
            bin_ds[i-1] = 0.
            bin_de1_err[i-1] = 0.
            bin_de2_err[i-1] = 0.
            bin_ds_err[i-1] = 0.
    print 'fixed nans'

    print 'index = ',index
    print 'bin_de1 = ',bin_de1
    print 'bin_de2 = ',bin_de2
    print 'bin_ds = ',bin_ds

    plt.clf()
    plt.title('PSF Size residuals')
    plt.xlim(0,xmax)
    plt.ylim(0,0.0025)
    plt.plot([0,xmax], [0,0], color='black')
    plt.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')
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
    ax.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

    ax2.set_xlim(xmax-200,xmax)
    ax2.set_ylim(0,0.0025)
    ax2.plot([0,xmax], [0,0], color='black')
    ax2.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

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


def bin_by_fov(ccd, x, y, ds, de1, de2, key):

    import numpy
    import matplotlib.pyplot as plt
    from toFocal import toFocal

    all_x = numpy.array([])
    all_y = numpy.array([])
    all_e1 = numpy.array([])
    all_e2 = numpy.array([])
    all_s = numpy.array([])

    nwhisk = 5
    x_bins = numpy.linspace(0,2048,nwhisk+1)
    y_bins = numpy.linspace(0,4096,2*nwhisk+1)

    ccdnums = numpy.unique(ccd)
    for ccdnum in ccdnums:
        mask = numpy.where(ccd == ccdnum)[0]
        print 'ccdnum = ',ccdnum,', nstar = ',mask.sum()
        if mask.sum() < 100: continue

        x_index = numpy.digitize(x[mask], x_bins)
        y_index = numpy.digitize(y[mask], y_bins)

        bin_de1 = numpy.array([ de1[mask][(x_index == i) & (y_index == j)].mean() 
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_de1 = ',bin_de1
        bin_de2 = numpy.array([ de2[mask][(x_index == i) & (y_index == j)].mean() 
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_de2 = ',bin_de2
        bin_ds = numpy.array([ ds[mask][(x_index == i) & (y_index == j)].mean() 
                    for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_ds = ',bin_ds
        bin_x = numpy.array([ x[mask][(x_index == i) & (y_index == j)].mean() 
                  for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_x = ',bin_x
        bin_y = numpy.array([ y[mask][(x_index == i) & (y_index == j)].mean() 
                  for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_y = ',bin_y
        bin_count = numpy.array([ ((x_index == i) & (y_index == j)).sum()
                      for i in range(1, len(x_bins)) for j in range(1, len(y_bins)) ])
        print 'bin_count = ',bin_count

        focal_x, focal_y = toFocal(ccdnum, bin_x, bin_y)

        mask2 = numpy.where(bin_count > 0)[0]
        print 'num with count > 0 = ',mask2.sum()
        all_x = numpy.append(all_x, focal_x[mask2])
        all_y = numpy.append(all_y, focal_y[mask2])
        all_e1 = numpy.append(all_e1, bin_de1[mask2])
        all_e2 = numpy.append(all_e2, bin_de2[mask2])
        all_s = numpy.append(all_s, bin_ds[mask2])


    plt.clf()
    #plt.title('PSF Ellipticity residuals in DES focal plane')
    theta = numpy.arctan2(all_e2,all_e1)/2.
    r = numpy.sqrt(all_e1**2 + all_e2**2)
    u = r*numpy.cos(theta)
    v = r*numpy.sin(theta)
    plt.xlim(-250,250)
    plt.ylim(-250,250)
    print 'all_x = ',all_x
    print 'len(all_x) = ',len(all_x)
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

 

def main():
    import os
    import glob
    import galsim
    import json
    import numpy
    import matplotlib
    matplotlib.use('Agg') # needs to be done before import pyplot

    args = parse_args()

    work = os.path.expanduser(args.work)
    print 'work dir = ',work

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

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
    s_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    pe1_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    pe2_list = { 'griz' : [], 'riz' : [], 'ri' : [] }
    ps_list = { 'griz' : [], 'riz' : [], 'ri' : [] }

    get_data(runs, exps, work,
             mask_list, used_list, ccd_list,
             airmass_list, sky_list, sigsky_list, fwhm_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, s_list, pe1_list, pe2_list, ps_list, args.tag, psfex='psfex')
             #e1_list, e2_list, s_list, pe1_list, pe2_list, ps_list, psfex='erin')

    #for key in ra_list.keys():
    for key in ['r', 'ri', 'riz']:
        mask = numpy.concatenate(mask_list[key])
        used = numpy.concatenate(used_list[key])
        ccd = numpy.concatenate(ccd_list[key])

        airmass = numpy.concatenate(airmass_list[key])
        sky = numpy.concatenate(sky_list[key])
        sigsky = numpy.concatenate(sigsky_list[key])
        fwhm = numpy.concatenate(fwhm_list[key])

        med_airmass = numpy.median(airmass)
        med_sky = numpy.median(sky)
        med_sigsky = numpy.median(sigsky)
        med_fwhm = numpy.median(fwhm)
        print 'airmass: ',min(airmass),med_airmass,max(airmass)
        print 'sky: ',min(sky),med_sky,max(sky)
        print 'sigsky: ',min(sigsky),med_sigsky,max(sigsky)
        print 'fwhm: ',min(fwhm),med_fwhm,max(fwhm)

        ra = numpy.concatenate(ra_list[key])
        dec = numpy.concatenate(dec_list[key])
        x = numpy.concatenate(x_list[key])
        y = numpy.concatenate(y_list[key])
        m = numpy.concatenate(m_list[key])
        print 'full mag range = ',numpy.min(m),numpy.max(m)
        print 'masked mag range = ',numpy.min(m[mask]),numpy.max(m[mask])
        print 'used mag range = ',numpy.min(m[used]),numpy.max(m[used])
        print 'len(m) = ',len(m)
        print 'len(mask) = ',len(mask)
        print 'len(used) = ',len(used)
        print 'sum(mask) = ',numpy.sum(mask)
        print 'sum(used) = ',numpy.sum(used)

        e1 = numpy.concatenate(e1_list[key])
        print 'mean e1 = ',numpy.mean(e1[mask])
        e2 = numpy.concatenate(e2_list[key])
        print 'mean e2 = ',numpy.mean(e2[mask])
        s = numpy.concatenate(s_list[key])
        print 'mean s = ',numpy.mean(s[mask])
        pe1 = numpy.concatenate(pe1_list[key])
        print 'mean pe1 = ',numpy.mean(pe1[mask])
        pe2 = numpy.concatenate(pe2_list[key])
        print 'mean pe2 = ',numpy.mean(pe2[mask])
        ps = numpy.concatenate(ps_list[key])
        print 'mean ps = ',numpy.mean(ps[mask])

        print 'min mag = ',numpy.min(m)
        print 'max mag = ',numpy.max(m)
        print 'mean s (used) = ',numpy.mean(s[used])
        print 'mean e1 (used) = ',numpy.mean(e1[used])
        print 'mean e2 (used) = ',numpy.mean(e2[used])
        print 'mean ps (used) = ',numpy.mean(ps[used])
        print 'mean pe1 (used) = ',numpy.mean(pe1[used])
        print 'mean pe2 (used) = ',numpy.mean(pe2[used])
        de1 = e1 - pe1
        de2 = e2 - pe2
        ds = s - ps
        dt = s**2 - ps**2
        print 'mean ds (used) = ',numpy.mean(ds[used])
        print 'mean de1 (used) = ',numpy.mean(de1[used])
        print 'mean de2 (used) = ',numpy.mean(de2[used])

        min_mused = numpy.min(m[used])
        bin_by_mag(m[mask], dt[mask], de1[mask], de2[mask], min_mused, key)
        bin_by_mag(m[used], dt[used], de1[used], de2[used], min_mused, key+'_used')
        #bin_by_mag(m[mask], s[mask], e1[mask], e2[mask], min_mused, key+'_orig')
        #bin_by_mag(m[mask], ps[mask], pe1[mask], pe2[mask], min_mused, key+'_model')

        #mask2 = mask & (sky<med_sky)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowsky')
        #mask2 = mask & (sky>med_sky)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_highsky')

        #mask2 = mask & (sigsky<med_sigsky)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowsigsky')
        #mask2 = mask & (sigsky>med_sigsky)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_highsigsky')

        #mask2 = mask & (airmass<med_airmass)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowairmass')
        #mask2 = mask & (airmass>med_airmass)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_highairmass')

        #mask2 = mask & (fwhm<med_fwhm)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_lowfwhm')
        #mask2 = mask & (fwhm>med_fwhm)
        #bin_by_mag(m[mask2], ds[mask2], de1[mask2], de2[mask2], min_mused, key+'_highfwhm')

        #bin_by_mag(m[mask], s[mask] - numpy.mean(s[mask]),
        #           e1[mask] - numpy.mean(e1[mask]), e2[mask] - numpy.mean(e2[mask]), min_mused, key)

        #bin_by_chip_pos(x[used], ds[used], de1[used], de2[used], key, 'x')
        #bin_by_chip_pos(y[used], ds[used], de1[used], de2[used], key, 'y')

        #bin_by_fov(ccd[used], x[used], y[used], ds[used], de1[used], de2[used], key)
        #bin_by_fov(ccd[used], x[used], y[used], s[used], e1[used], e2[used], key + 'raw')


if __name__ == "__main__":
    main()
