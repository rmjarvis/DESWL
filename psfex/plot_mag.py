#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits*',
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

def get_data(runs, exps, work,
             mask_list, used_list, ccd_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, s_list, pe1_list, pe2_list, ps_list):

    import astropy.io.fits as pyfits
    import numpy
    import os

    expinfo_file = 'exposure_info.fits'
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    for run,exp in zip(runs,exps):

        print 'Start work on run, exp = ',run,exp
        expnum = int(exp[6:])
        print 'expnum = ',expnum

        if expnum not in expinfo['expnum']:
            print 'expnum is not in expinfo!'
            print 'expinfo[expnum] = ',expinfo['expnum']
            raise RuntimeError("Could not find information about this expnum")
        k = numpy.nonzero(expinfo['expnum'] == expnum)[0][0]
        print 'k = ',k
        filter = expinfo['filter'][k]
        print 'filter[k] = ',filter

        exp_dir = os.path.join(work,exp)
        print 'exp_dir = ',exp_dir

        cat_file = os.path.join(exp_dir, exp + "_psf.fits")
        with pyfits.open(cat_file) as pyf:
            data = pyf[1].data.copy()

        mask = (data['flag'] <= 1)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        used = (data['flag'] == 0)

        mask_list['griz'].append(mask)
        used_list['griz'].append(used)
        ccd_list['griz'].append(data['ccdnum'])

        ra_list['griz'].append(data['ra'])
        dec_list['griz'].append(data['dec'])
        x_list['griz'].append(data['x'])
        y_list['griz'].append(data['y'])
        m_list['griz'].append(data['mag'])

        e1_list['griz'].append(data['e1'])
        e2_list['griz'].append(data['e2'])
        s_list['griz'].append(data['size'])

        pe1_list['griz'].append(data['psfex_e1'])
        pe2_list['griz'].append(data['psfex_e2'])
        ps_list['griz'].append(data['psfex_size'])

        if filter not in ra_list.keys():
            mask_list[filter] = []
            used_list[filter] = []
            ccd_list[filter] = []
            ra_list[filter] = []
            dec_list[filter] = []
            x_list[filter] = []
            y_list[filter] = []
            m_list[filter] = []
            e1_list[filter] = []
            e2_list[filter] = []
            s_list[filter] = []
            pe1_list[filter] = []
            pe2_list[filter] = []
            ps_list[filter] = []

        if 'g' not in filter:
            mask_list['riz'].append(mask)
            used_list['riz'].append(used)
            ccd_list['riz'].append(data['ccdnum'])

            ra_list['riz'].append(data['ra'])
            dec_list['riz'].append(data['dec'])
            x_list['riz'].append(data['x'])
            y_list['riz'].append(data['y'])
            m_list['riz'].append(data['mag'])

            e1_list['riz'].append(data['e1'])
            e2_list['riz'].append(data['e2'])
            s_list['riz'].append(data['size'])

            pe1_list['riz'].append(data['psfex_e1'])
            pe2_list['riz'].append(data['psfex_e2'])
            ps_list['riz'].append(data['psfex_size'])

        mask_list[filter].append(mask)
        used_list[filter].append(used)
        ccd_list[filter].append(data['ccdnum'])

        ra_list[filter].append(data['ra'])
        dec_list[filter].append(data['dec'])
        x_list[filter].append(data['x'])
        y_list[filter].append(data['y'])
        m_list[filter].append(data['mag'])

        e1_list[filter].append(data['e1'])
        e2_list[filter].append(data['e2'])
        s_list[filter].append(data['size'])

        pe1_list[filter].append(data['psfex_e1'])
        pe2_list[filter].append(data['psfex_e2'])
        ps_list[filter].append(data['psfex_size'])

    print '\nFinished processing all exposures'


def bin_by_mag(m, ds, de1, de2, key):

    import numpy
    import matplotlib.pyplot as plt

    # Bin by mag
    mag_bins = numpy.linspace(10,17,71)
    print 'mag_bins = ',mag_bins

    index = numpy.digitize(m, mag_bins)
    print 'len(index) = ',len(index)
    bin_de1 = [de1[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_de1 = ',bin_de1
    bin_de2 = [de2[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_de2 = ',bin_de2
    bin_ds = [ds[index == i].mean() for i in range(1, len(mag_bins))]
    print 'bin_ds = ',bin_ds
    bin_de1_err = [ numpy.sqrt(de1[index == i].var() / len(de1[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_de1_err = ',bin_de1_err
    bin_de2_err = [ numpy.sqrt(de2[index == i].var() / len(de2[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_de2_err = ',bin_de2_err
    bin_ds_err = [ numpy.sqrt(ds[index == i].var() / len(ds[index == i])) 
                    for i in range(1, len(mag_bins)) ]
    print 'bin_ds_err = ',bin_ds_err

    # Fix up nans
    for i in range(1,len(mag_bins)):
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
    plt.xlim(10,16.5)
    plt.ylim(-0.002,0.012)
    plt.plot([10,16.5], [0,0], color='black')
    plt.plot([13,13],[-1,1], color='red')
    plt.fill( [10,10,13,13], [-1,1,1,-1], fill=False, hatch='/', color='red')
    plt.errorbar(mag_bins[:-1], bin_ds, yerr=bin_ds_err, color='blue', fmt='o')
    plt.xlabel('Uncorrected Magnitude')
    plt.ylabel('$size_{psf} - size_{model}$')
    plt.savefig('dsize_mag_' + key + '.png')
    plt.savefig('dsize_mag_' + key + '.pdf')

    plt.clf()
    plt.title('PSF Ellipticity residuals')
    plt.xlim(10,16.5)
    plt.ylim(-3.e-4,8.e-4)
    plt.plot([10,16.5], [0,0], color='black')
    plt.plot([13,13],[-1,1], color='red')
    plt.fill( [10,10,13,13], [-1,1,1,-1], fill=False, hatch='/', color='red')
    e1_line = plt.errorbar(mag_bins[:-1], bin_de1, yerr=bin_de1_err, color='blue', fmt='o')
    e2_line = plt.errorbar(mag_bins[:-1], bin_de2, yerr=bin_de2_err, color='green', fmt='o')
    plt.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
    plt.xlabel('Uncorrected Magnitude')
    plt.ylabel('$e_{psf} - e_{model}$')
    plt.savefig('de_mag_' + key + '.png')
    plt.savefig('de_mag_' + key + '.pdf')


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


def bin_by_fov(ccd, x, y, ds, de1, de2, key):

    import numpy
    import matplotlib.pyplot as plt
    from toFocal import toFocal

    all_x = numpy.array([])
    all_y = numpy.array([])
    all_e1 = numpy.array([])
    all_e2 = numpy.array([])
    all_s = numpy.array([])

    x_bins = numpy.linspace(0,2048,9)
    y_bins = numpy.linspace(0,4096,17)

    ccdnums = numpy.unique(ccd)
    for ccdnum in ccdnums:
        mask = (ccd == ccdnum)
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

        mask2 = (bin_count > 0)
        print 'num with count > 0 = ',mask2.sum()
        all_x = numpy.append(all_x, focal_x[mask2])
        all_y = numpy.append(all_y, focal_y[mask2])
        all_e1 = numpy.append(all_e1, bin_de1[mask2])
        all_e2 = numpy.append(all_e2, bin_de2[mask2])
        all_s = numpy.append(all_s, bin_ds[mask2])


    plt.clf()
    plt.title('PSF Ellipticity residuals in DES focal plane')
    theta = numpy.arctan2(all_e2,all_e1)/2.
    r = numpy.sqrt(all_e1**2 + all_e2**2)
    u = r*numpy.cos(theta)
    v = r*numpy.sin(theta)
    plt.xlim(-250,250)
    plt.ylim(-250,250)
    qv = plt.quiver(all_x,all_y,u,v, pivot='middle', scale_units='xy',
                    headwidth=0., headlength=0., headaxislength=0.,
                    width=0.001, scale=6.e-4, color='blue')
    plt.quiverkey(qv, 0.10, 0.08, 0.01, "de = " + str(0.01),
                  coordinates='axes', color='darkred', labelcolor='darkred',
                  labelpos='E', fontproperties={'size':'x-small'})
    plt.axis('off')
    plt.savefig('de_fov_' + key + '.png')
    plt.savefig('de_fov_' + key + '.pdf')

 

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

    mask_list = { 'griz' : [], 'riz' : [] }
    used_list = { 'griz' : [], 'riz' : [] }
    ccd_list = { 'griz' : [], 'riz' : [] }

    ra_list = { 'griz' : [], 'riz' : [] }
    dec_list = { 'griz' : [], 'riz' : [] }
    x_list = { 'griz' : [], 'riz' : [] }
    y_list = { 'griz' : [], 'riz' : [] }
    m_list = { 'griz' : [], 'riz' : [] }

    e1_list = { 'griz' : [], 'riz' : [] }
    e2_list = { 'griz' : [], 'riz' : [] }
    s_list = { 'griz' : [], 'riz' : [] }

    pe1_list = { 'griz' : [], 'riz' : [] }
    pe2_list = { 'griz' : [], 'riz' : [] }
    ps_list = { 'griz' : [], 'riz' : [] }

    get_data(runs, exps, work,
             mask_list, used_list, ccd_list,
             ra_list, dec_list, x_list, y_list, m_list,
             e1_list, e2_list, s_list, pe1_list, pe2_list, ps_list)

    for key in ra_list.keys():
    #for key in ['riz']:
        mask = numpy.concatenate(mask_list[key])
        used = numpy.concatenate(used_list[key])
        ccd = numpy.concatenate(ccd_list[key])

        ra = numpy.concatenate(ra_list[key])
        dec = numpy.concatenate(dec_list[key])
        x = numpy.concatenate(x_list[key])
        y = numpy.concatenate(y_list[key])
        m = numpy.concatenate(m_list[key])

        e1 = numpy.concatenate(e1_list[key])
        e2 = numpy.concatenate(e2_list[key])
        s = numpy.concatenate(s_list[key])
        pe1 = numpy.concatenate(pe1_list[key])
        pe2 = numpy.concatenate(pe2_list[key])
        ps = numpy.concatenate(ps_list[key])

        print 'min mag = ',numpy.min(m)
        print 'max mag = ',numpy.max(m)
        print 'mean s = ',numpy.mean(s[mask])
        print 'mean e1 = ',numpy.mean(e1[mask])
        print 'mean e2 = ',numpy.mean(e2[mask])
        print 'mean s (used) = ',numpy.mean(s[used])
        print 'mean e1 (used) = ',numpy.mean(e1[used])
        print 'mean e2 (used) = ',numpy.mean(e2[used])
        de1 = e1 - pe1
        de2 = e2 - pe2
        ds = s - ps

        bin_by_mag(m[mask], ds[mask], de1[mask], de2[mask], key)

        bin_by_chip_pos(x[used], ds[used], de1[used], de2[used], key, 'x')
        bin_by_chip_pos(y[used], ds[used], de1[used], de2[used], key, 'y')

        bin_by_fov(ccd[used], x[used], y[used], ds[used], de1[used], de2[used], key)


if __name__ == "__main__":
    main()
