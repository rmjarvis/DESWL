#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import os
import numpy as np
import fitsio
from toFocal import toFocal

def read_data(exps, work, keys, limit_bands=None, prefix='piff', use_reserved=False):

    RESERVED = 64
    BAD_CCDS = [2, 31, 61]
    MAX_TILING = 10

    all_keys = keys

    if 'x' in keys:
        all_keys += ['fov_x', 'fov_y']

    all_keys += ['exp', 'ccd', 'band', 'tiling']

    all_data = { key : [] for key in all_keys }

    bands = set()   # This is the set of all bands being used
    tilings = set()   # This is the set of all tilings being used

    n_reject_mean_dt = 0
    n_reject_mean_de1 = 0
    n_reject_mean_de2 = 0
    n_reject_mean_e1 = 0
    n_reject_mean_e2 = 0
    n_reject_std_dt = 0
    n_reject_std_de1 = 0
    n_reject_std_de2 = 0
    n_reject_rho2 = 0

    nrows = 0

    for exp in exps:

        expnum = int(exp)
        try:
            expinfo = fitsio.read(os.path.join(work, exp, 'exp_info_%d.fits'%expnum))
        except Exception as e:
            #print('Caught: ',e)
            #print('Skip this exposure')
            continue

        if expnum not in expinfo['expnum']:
            print('expnum is not in expinfo!')
            print('expinfo[expnum] = ',expinfo['expnum'])
            print('Could not find information about this expnum.  Skipping ',run,exp)
            continue
        i = np.nonzero(expinfo['expnum'] == expnum)[0][0]
        #print('i = ',i)
        band = expinfo['band'][i]
        if (limit_bands is not None) and (band not in limit_bands):
            #print('Not doing band = %s.'%band)
            continue

        print('Start work on exp = ',exp)
        print('band = ',band)

        if 'tiling' in expinfo:
            tiling = int(expinfo['tiling'][i])
            if tiling == 0:
                # This shouldn't happen, but it did for a few exposures.  Just skip them, since this
                # might indicate some kind of problem.
                print('tiling == 0.  Skip this exposure.')
                continue
            if tiling > MAX_TILING
                print('tiling is > %d.  Skip this exposure.'%MAX_TILING)
                continue
            print('tiling = ',tiling)
        else:
            tiling = 0

        for k in range(len(expinfo)):
            ccdnum = expinfo[k]['ccdnum']
            if expinfo[k]['flag'] != 0:
                #print('Skipping ccd %d because it is blacklisted: '%ccdnum, expinfo[k]['flag'])
                continue
            if ccdnum in BAD_CCDS:
                #print('Skipping ccd %d because it is BAD'%ccdnum)
                continue

            cat_file = os.path.join(work, exp, "psf_cat_%d_%d.fits"%(expnum,ccdnum))
            #print('cat_file = ',cat_file)
            try:
                data = fitsio.read(cat_file)
                flag = data[prefix+'_flag']
            except (OSError, IOError):
                #print('Unable to open cat_file %s.  Skipping this file.'%cat_file)
                continue

            ntot = len(data)
            nused = np.sum((flag & 1) != 0)
            nreserved = np.sum((flag & RESERVED) != 0)
            ngood = np.sum(flag == 0)
            #print('nused = ',nused)
            #print('nreserved = ',nreserved)
            #print('ngood = ',ngood)

            if use_reserved:
                mask = flag & ~(RESERVED+1) == RESERVED
            else:
                mask = flag & ~1 == 0
            used = flag == 0
            #print('mask = ',mask)

            T = data['obs_T']
            e1 = data['obs_e1']
            e2 = data['obs_e2']
            dT = data['obs_T'] - data[prefix + '_T']
            de1 = data['obs_e1'] - data[prefix + '_e1']
            de2 = data['obs_e2'] - data[prefix + '_e2']
            print(expnum, ccdnum, len(dT), band)
            #print('T = ',np.mean(T[used]),np.std(T[used]))
            #print('e1 = ',np.mean(e1[used]),np.std(e1[used]))
            #print('e2 = ',np.mean(e2[used]),np.std(e2[used]))
            #print('dT/T = ',np.mean(dT[used]/T[used]),np.std(dT[used]/T[used]))
            #print('de1 = ',np.mean(de1[used]),np.std(de1[used]))
            #print('de2 = ',np.mean(de2[used]),np.std(de2[used]))
            rho2 = (e1 - 1j*e2) * (de1 + 1j*de2)
            #print('mean rho2 = ',np.mean(rho2))
            if abs(np.mean(dT[used]/T[used])) > 0.01:
                print('mean dT/T = %f on ccd %d.'%(np.mean(dT[used]/T[used]),ccdnum))
                n_reject_mean_dt += 1
                #continue
            if abs(np.mean(de1[used])) > 0.01:
                print('mean de1 = %f on ccd %d.'%(np.mean(de1[used]),ccdnum))
                n_reject_mean_de1 += 1
                #continue
            if abs(np.mean(de2[used])) > 0.01:
                print('mean de2 = %f on ccd %d.'%(np.mean(de2[used]),ccdnum))
                n_reject_mean_de2 += 1
                #continue
            if abs(np.std(dT[used]/T[used])) > 0.1:
                print('std dT/T = %f on ccd %d.'%(np.std(dT[used]/T[used]),ccdnum))
                n_reject_std_dt += 1
                #continue
            if abs(np.std(de1[used])) > 0.1:
                print('std de1 = %f on ccd %d.'%(np.std(de1[used]),ccdnum))
                n_reject_std_de1 += 1
                #continue
            if abs(np.std(de2[used])) > 0.1:
                print('std de2 = %f on ccd %d.'%(np.std(de2[used]),ccdnum))
                n_reject_std_de2 += 1
                #continue
            if abs(np.mean(rho2)) > 5.e-4:
                print('mean rho2 = %s on ccd %d.'%(np.mean(rho2),ccdnum))
                n_reject_rho2 += 1
                #continue
            if abs(np.mean(e1[used])) > 0.03:
                print('mean e1 = %f on ccd %d.'%(np.mean(e1[used]),ccdnum))
                n_reject_mean_e1 += 1
                #continue
            if abs(np.mean(e2[used])) > 0.03:
                print('mean e2 = %f on ccd %d.'%(np.mean(e2[used]),ccdnum))
                n_reject_mean_e2 += 1
                #continue

            # Filter out egregiously bad values.  Just in case.
            good = (abs(dT/T) < 0.1) & (abs(de1) < 0.1) & (abs(de2) < 0.1)
            mask = mask & good

            ngood = np.sum(mask)
            #print('ngood = ',ngood,'/',len(data))
            assert ngood == len(data[mask])
            if ngood == 0:
                print('All objects in ccd %d are flagged.'%ccdnum)
                print('Probably due to astrometry flags. Skip this exposure.')
                continue

            # Start with just the input keys, which should be columns in data.
            for key in keys:
                all_data[key].append(data[key][mask])

            # Now add the extra ones we added to all_keys
            if 'x' in keys:
                # Convert to focal position.
                x,y = toFocal(ccdnum, data['x'][mask], data['y'][mask])
                # This comes back in units of mm.  Convert to arcsec.
                # 1 pixel = 15e-3 mm = 0.263 arcsec
                x *= 0.263/15e-3
                y *= 0.263/15e-3
                all_data['fov_x'].append(x)
                all_data['fov_y'].append(y)

            all_data['exp'].append([expnum] * ngood)
            all_data['ccd'].append([ccdnum] * ngood)
            all_data['band'].append([band] * ngood)
            all_data['tiling'].append([tiling] * ngood)
            bands.add(band)
            tilings.add(tiling)
            nrows += ngood

    print('\nFinished processing %d exposures'%len(exp))
    print('bands = ',bands)
    print('tilings = ',tilings)
    print('total good stars = ',nrows)

    print('Potential rejections: (not enabled):')
    print('n_reject_mean_dt = ',n_reject_mean_dt)
    print('n_reject_mean_de1 = ',n_reject_mean_de1)
    print('n_reject_mean_de2 = ',n_reject_mean_de2)
    print('n_reject_mean_e1 = ',n_reject_mean_de1)
    print('n_reject_mean_e2 = ',n_reject_mean_de2)
    print('n_reject_std_dt = ',n_reject_std_dt)
    print('n_reject_std_de1 = ',n_reject_std_de1)
    print('n_reject_std_de2 = ',n_reject_std_de2)
    print('n_reject_rho2 = ',n_reject_rho2)

    # Turn the data into a recarray
    # Pick appropriate formats for each kind of data
    formats = []
    for key in all_keys:
        if key == 'ccd' or key == 'tiling':
            formats.append('i2')
        elif key == 'exp' or 'flag' in key:
            formats.append('i4')
        elif key == 'band':
            formats.append('a1')
        else:
            formats.append('f8')
    data = np.recarray(shape=(nrows,), formats=formats, names=all_keys)
    #print('data.dtype = ',data.dtype)
    for key in all_keys:
        data[key] = np.concatenate(all_data[key])
    print('made recarray')

    return data, bands, tilings
