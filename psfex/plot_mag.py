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


def main():
    import os
    import glob
    import galsim
    import json
    import numpy
    import astropy.io.fits as pyfits
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

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

    expinfo_file = 'exposure_info.fits'
    with pyfits.open(expinfo_file) as pyf:
        expinfo = pyf[1].data

    mask_list = { 'griz' : [], 'riz' : [] }
    used_list = { 'griz' : [], 'riz' : [] }

    ra_list = { 'griz' : [], 'riz' : [] }
    dec_list = { 'griz' : [], 'riz' : [] }
    m_list = { 'griz' : [], 'riz' : [] }
    x_list = { 'griz' : [], 'riz' : [] }
    y_list = { 'griz' : [], 'riz' : [] }

    e1_list = { 'griz' : [], 'riz' : [] }
    e2_list = { 'griz' : [], 'riz' : [] }
    s_list = { 'griz' : [], 'riz' : [] }

    pe1_list = { 'griz' : [], 'riz' : [] }
    pe2_list = { 'griz' : [], 'riz' : [] }
    ps_list = { 'griz' : [], 'riz' : [] }

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

        ccdnums = numpy.unique(data['ccdnum'])
        #print 'ccdnums = ',ccdnums

        stats = []

        mask = (data['flag'] <= 1)
        if mask.sum() == 0:
            print 'All objects in this exposure are flagged.'
            print 'Probably due to astrometry flags. Skip this exposure.'
            continue

        used = (data['flag'] == 0)

        mask_list['griz'].append(mask)
        used_list['griz'].append(used)

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

    #for key in ra_list.keys():
    for key in ['riz']:
        mask = numpy.concatenate(mask_list[key])
        used = numpy.concatenate(used_list[key])

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

        # Bin by mag
        mag_bins = numpy.linspace(10,17,71)
        print 'mag_bins = ',mag_bins

        index = numpy.digitize(m[mask], mag_bins)
        print 'len(index) = ',len(index)
        bin_de1 = [de1[mask][index == i].mean() for i in range(1, len(mag_bins))]
        print 'bin_de1 = ',bin_de1
        bin_de2 = [de2[mask][index == i].mean() for i in range(1, len(mag_bins))]
        print 'bin_de2 = ',bin_de2
        bin_ds = [ds[mask][index == i].mean() for i in range(1, len(mag_bins))]
        print 'bin_ds = ',bin_ds
        bin_de1_err = [ numpy.sqrt(de1[mask][index == i].var() / len(de1[mask][index == i])) 
                        for i in range(1, len(mag_bins)) ]
        print 'bin_de1_err = ',bin_de1_err
        bin_de2_err = [ numpy.sqrt(de2[mask][index == i].var() / len(de2[mask][index == i])) 
                        for i in range(1, len(mag_bins)) ]
        print 'bin_de2_err = ',bin_de2_err
        bin_ds_err = [ numpy.sqrt(ds[mask][index == i].var() / len(ds[mask][index == i])) 
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

        # Bin by chip x position
        x_bins = numpy.linspace(0,2048,129)
        print 'x_bins = ',x_bins
        index = numpy.digitize(x[used], x_bins)
        print 'len(index) = ',len(index)
        bin_de1 = [de1[used][index == i].mean() for i in range(1, len(x_bins))]
        print 'bin_de1 = ',bin_de1
        bin_de2 = [de2[used][index == i].mean() for i in range(1, len(x_bins))]
        print 'bin_de2 = ',bin_de2
        bin_ds = [ds[used][index == i].mean() for i in range(1, len(x_bins))]
        print 'bin_ds = ',bin_ds
        bin_de1_err = [ numpy.sqrt(de1[used][index == i].var() / len(de1[used][index == i])) 
                        for i in range(1, len(x_bins)) ]
        print 'bin_de1_err = ',bin_de1_err
        bin_de2_err = [ numpy.sqrt(de2[used][index == i].var() / len(de2[used][index == i])) 
                        for i in range(1, len(x_bins)) ]
        print 'bin_de2_err = ',bin_de2_err
        bin_ds_err = [ numpy.sqrt(ds[used][index == i].var() / len(ds[used][index == i])) 
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
        plt.xlim(0,2048)
        plt.ylim(0,0.0025)
        plt.plot([0,2048], [0,0], color='black')
        plt.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')
        plt.xlabel('Chip x position')
        plt.ylabel('$size_{psf} - size_{model}$')
        plt.savefig('dsize_x_' + key + '.png')
        plt.savefig('dsize_x_' + key + '.pdf')

        plt.clf()
        plt.title('PSF Ellipticity residuals')
        plt.xlim(0,2048)
        plt.ylim(-0.002,0.002)
        plt.plot([0,2048], [0,0], color='black')
        e1_line = plt.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        e2_line = plt.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
        plt.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
        plt.xlabel('Chip x position')
        plt.ylabel('$e_{psf} - e_{model}$')
        plt.savefig('de_x_' + key + '.png')
        plt.savefig('de_x_' + key + '.pdf')

        # Make broken x-axis. 
        plt.clf()
        f, (ax,ax2) = plt.subplots(1,2,sharey=True)
        ax.set_xlim(0,200)
        ax.set_ylim(0,0.0025)
        ax.plot([0,2048], [0,0], color='black')
        ax.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

        ax2.set_xlim(2048-200,2048)
        ax2.set_ylim(0,0.0025)
        ax2.plot([0,2048], [0,0], color='black')
        ax2.errorbar(x_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

        break_axes(ax,ax2)

        ax.set_title('PSF Size residuals', x=1.08, y=1.03)
        ax.set_xlabel('Chip x position')
        ax.xaxis.set_label_coords(1.08,-0.05)
        ax.set_ylabel('$size_{psf} - size_{model}$')

        plt.savefig('dsize_x2_' + key + '.png')
        plt.savefig('dsize_x2_' + key + '.pdf')

        plt.clf()
        f, (ax,ax2) = plt.subplots(1,2,sharey=True)
        ax.set_xlim(0,200)
        ax.set_ylim(-0.002,0.002)
        ax.plot([0,2048], [0,0], color='black')
        ax.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        ax.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')

        ax2.set_xlim(2048-200,2048)
        ax2.set_ylim(-0.002,0.002)
        ax2.plot([0,2048], [0,0], color='black')
        e1_line = ax2.errorbar(x_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        e2_line = ax2.errorbar(x_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
        ax2.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])

        break_axes(ax,ax2)

        ax.set_title('PSF Ellipticity residuals', x=1.08, y=1.03)
        ax.set_xlabel('Chip x position')
        ax.xaxis.set_label_coords(1.08,-0.05)
        ax.set_ylabel('$e_{psf} - e_{model}$')

        plt.savefig('de_x2_' + key + '.png')
        plt.savefig('de_x2_' + key + '.pdf')


        # Bin by chip y position
        y_bins = numpy.linspace(0,4096,257)
        print 'y_bins = ',y_bins
        index = numpy.digitize(y[used], y_bins)
        print 'len(index) = ',len(index)
        bin_de1 = [de1[used][index == i].mean() for i in range(1, len(y_bins))]
        print 'bin_de1 = ',bin_de1
        bin_de2 = [de2[used][index == i].mean() for i in range(1, len(y_bins))]
        print 'bin_de2 = ',bin_de2
        bin_ds = [ds[used][index == i].mean() for i in range(1, len(y_bins))]
        print 'bin_ds = ',bin_ds
        bin_de1_err = [ numpy.sqrt(de1[used][index == i].var() / len(de1[used][index == i])) 
                        for i in range(1, len(y_bins)) ]
        print 'bin_de1_err = ',bin_de1_err
        bin_de2_err = [ numpy.sqrt(de2[used][index == i].var() / len(de2[used][index == i])) 
                        for i in range(1, len(y_bins)) ]
        print 'bin_de2_err = ',bin_de2_err
        bin_ds_err = [ numpy.sqrt(ds[used][index == i].var() / len(ds[used][index == i])) 
                       for i in range(1, len(y_bins)) ]
        print 'bin_ds_err = ',bin_ds_err

        # Fix up nans
        for i in range(1,len(y_bins)):
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
        plt.xlim(0,4096)
        plt.ylim(0,0.002)
        plt.plot([0,4096], [0,0], color='black')
        plt.errorbar(y_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')
        plt.ylabel('Chip y position')
        plt.ylabel('$size_{psf} - size_{model}$')
        plt.savefig('dsize_y_' + key + '.png')
        plt.savefig('dsize_y_' + key + '.pdf')

        plt.clf()
        plt.title('PSF Ellipticity residuals')
        plt.xlim(0,4096)
        plt.ylim(-0.002,0.002)
        plt.plot([0,4096], [0,0], color='black')
        e1_line = plt.errorbar(y_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        e2_line = plt.errorbar(y_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
        plt.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
        plt.ylabel('Chip y position')
        plt.ylabel('$e_{psf} - e_{model}$')
        plt.savefig('de_y_' + key + '.png')
        plt.savefig('de_y_' + key + '.pdf')

        # Make broken x-axis. 
        plt.clf()
        f, (ax,ax2) = plt.subplots(1,2,sharey=True)
        ax.set_xlim(0,200)
        ax.set_ylim(0,0.0025)
        ax.plot([0,4096], [0,0], color='black')
        ax.errorbar(y_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

        ax2.set_xlim(4096-200,4096)
        ax2.set_ylim(0,0.0025)
        ax2.plot([0,4096], [0,0], color='black')
        ax2.errorbar(y_bins[2:-3], bin_ds[2:-2], yerr=bin_ds_err[2:-2], color='blue', fmt='o')

        break_axes(ax,ax2)

        ax.set_title('PSF Size residuals', x=1.08, y=1.03)
        ax.set_xlabel('Chip y position')
        ax.xaxis.set_label_coords(1.08,-0.05)
        ax.set_ylabel('$size_{psf} - size_{model}$')

        plt.savefig('dsize_y2_' + key + '.png')
        plt.savefig('dsize_y2_' + key + '.pdf')

        plt.clf()
        f, (ax,ax2) = plt.subplots(1,2,sharey=True)
        ax.set_xlim(0,200)
        ax.set_ylim(-0.002,0.002)
        ax.plot([0,4096], [0,0], color='black')
        ax.errorbar(y_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        ax.errorbar(y_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')

        ax2.set_xlim(4096-200,4096)
        ax2.set_ylim(-0.002,0.002)
        ax2.plot([0,4096], [0,0], color='black')
        e1_line = ax2.errorbar(y_bins[2:-3], bin_de1[2:-2], yerr=bin_de1_err[2:-2], color='blue', fmt='o')
        e2_line = ax2.errorbar(y_bins[2:-3], bin_de2[2:-2], yerr=bin_de2_err[2:-2], color='green', fmt='o')
        ax2.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])

        break_axes(ax,ax2)

        ax.set_title('PSF Ellipticity residuals', x=1.08, y=1.03)
        ax.set_xlabel('Chip y position')
        ax.xaxis.set_label_coords(1.08,-0.05)
        ax.set_ylabel('$e_{psf} - e_{model}$')

        plt.savefig('de_y2_' + key + '.png')
        plt.savefig('de_y2_' + key + '.pdf')


if __name__ == "__main__":
    main()
