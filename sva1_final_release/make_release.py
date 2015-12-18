#!/usr/bin/env python
"""
Based on Matt's make_flatcats.py, but changed the names of some of the column
names and added/subtracted a few columns.
"""

import os
import sys
import numpy as np
import fitsio
import numpy_util
import unblind

def make_info(config):
    """Make the WL Info file

    I define the flags a bit differently here than Matt did.  He only wrote out objects
    that had sva1_gold_flags <= 3, so he didn't have all the possible flag values for sva1_flag
    that we claimed in the paper.

    Also, rather than keep sva1_gold_mag_flags as a separate column with values of 1 or 0,
    I add that on as another bit onto sva1_gold_flags.  sva1_flag = 4096 is what Matt had
    declared sva1_gold_mag_flags = 1.
    """

    # Read in the gold_flags file.  This has most of what we need for the info file.
    flags = fitsio.read(config['gold_flags'])

    # Apparently, we don't have any objects outside the spte region.  So ignore this one
    # from here on out.
    assert(np.all(flags['sva1_spte_flags']==0))

    # Make sva1_gold_mag_flags the 4096 bit in sva1_flags
    sva1_flag = flags['sva1_gold_flags']
    sva1_flag[flags['sva1_gold_mag_flags'] == 1] |= 4096

    # We add in columns to tell when im3shape and ngmix have valid measurements.
    # So read in those files.
    im3shape_flag = fitsio.read(config['im3shape_flags'])['im3shape_flags']
    ngmix_flag = fitsio.read(config['ngmix_flags'])['ngmix_flags']

    # These flag values are either 0 or 4.  Switch to 0 or 1.
    im3shape_flag = im3shape_flag > 0
    ngmix_flag = ngmix_flag > 0

    # Setup the dtype:
    dt = np.dtype([
        ('coadd_objects_id', 'i8'),
        ('ra', 'f8'),
        ('dec', 'f8'),
        ('mag_auto_g', 'f4'),
        ('mag_auto_r', 'f4'),
        ('mag_auto_i', 'f4'),
        ('mag_auto_z', 'f4'),
        ('photoz_bin', 'i4'),
        ('mean_photoz', 'f4'),
        ('sva1_flag', 'i4'),
        ('im3shape_flag', 'i4'),
        ('ngmix_flag', 'i4') ])

    # Copy in the columns from the flags file
    data = np.zeros(len(flags), dtype=dt)
    for col in ['coadd_objects_id', 'ra', 'dec', 
                'mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z']:
        data[col] = flags[col]

    # Copy in the other flag columns
    data['sva1_flag'] = sva1_flag
    data['im3shape_flag'] = im3shape_flag
    data['ngmix_flag'] = ngmix_flag

    # Read the photoz binning file:
    pzbin = np.load(config['photoz_binning_ngmix'])
    pz_data = np.empty(len(pzbin),
            dtype=[('coadd_objects_id','i8'),('photoz_bin','i4'),('mean_photoz','f8')])    
    pz_data['coadd_objects_id'] = pzbin[:,0]
    pz_data['photoz_bin'] = pzbin[:,1]
    pz_data['mean_photoz'] = pzbin[:,2]

    # Make sure all pz ids are in the main file:
    assert np.all(np.in1d(pz_data['coadd_objects_id'], data['coadd_objects_id']))
    # Also, both id lists should already be sorted.
    assert np.all(pz_data['coadd_objects_id'] == sorted(pz_data['coadd_objects_id']))
    assert np.all(data['coadd_objects_id'] == sorted(data['coadd_objects_id']))

    # Set reasonable defaults for photoz_bin and mean_photoz
    data['photoz_bin'] = -99
    data['mean_photoz'] = -99.0

    pz_index = np.in1d(data['coadd_objects_id'], pz_data['coadd_objects_id'])
    data['photoz_bin'][pz_index] = pz_data['photoz_bin']
    data['mean_photoz'][pz_index] = pz_data['mean_photoz']
    
    # Before writing out, make all the names upper case
    data.dtype.names = tuple([name.upper() for name in data.dtype.names])

    fitsio.write(config['release_info'],data,clobber=True)
    return data
    

def make_ngmix(config, info):
    """Make the ngmix catalog file
    """

    # Setup the dtype:
    dt = np.dtype([
        ('coadd_objects_id', 'i8'),
        ('e_1', 'f8'),
        ('e_2', 'f8'),
        ('sens_avg', 'f8'),
        ('w', 'f8'),
        ('e_cov_1_1', 'f8'),
        ('e_cov_1_2', 'f8'),
        ('e_cov_2_1', 'f8'),
        ('e_cov_2_2', 'f8'),
        ('error_flag', 'i4'),
        ('snr_w', 'f4'),
        ('snr_r', 'f4'),
        ('flux_i', 'f4'),
        ('mag_i', 'f4'),
        ('t', 'f4'),
        ('t_err', 'f4'),
        ('t_r', 'f4'),
        ('snr_t', 'f4'),
        ('snr_t_r', 'f4'),
        ('log10_sb_i', 'f4'),
        ('mean_psf_e1', 'f8'),
        ('mean_psf_e2', 'f8'),
        ('mean_psf_t', 'f4'),
        ('sens_1', 'f4'),
        ('sens_2', 'f4'),
        ('arate', 'f4'),
        ('stamp_size', 'i4'),
        ('mask_frac', 'f4') ])

    # Read the source catalogs
    ngcat = fitsio.read(config['ngmix'])

    # ngmix catalogs are not sorted, so need to use match function
    # Also, ngmix has some object that are not in info catalog.
    ng, dng = numpy_util.match(ngcat['coadd_objects_id'], info['COADD_OBJECTS_ID'])

    # Set reasonable defaults for float fields:
    data = np.zeros(len(info), dtype=dt)
    for col, t in dt.descr:
        if 'f' in t:
            data[col] = -9999

    # The coadd_objects_id column can match for all rows
    data['coadd_objects_id'] = info['COADD_OBJECTS_ID']

    # Check that the overlap with ngcat is correct
    assert np.all(data['coadd_objects_id'][dng] == ngcat['coadd_objects_id'][ng])

    # Copy in the columns from the source catalog that keep the same name:
    for col in ['mask_frac']:
        data[col][dng] = ngcat[col][ng]

    # Some columns just need to drop the 'exp_' prefix:
    for col in ['e_1', 'e_2', 'flux_i', 'mag_i', 'arate', 't', 't_err',
                'e_cov_1_1', 'e_cov_1_2', 'e_cov_2_2']:
        data[col][dng] = ngcat['exp_' + col][ng]

    # Some need to be renamed:
    data['snr_w'][dng] = ngcat['exp_s2n_w'][ng]
    data['snr_t'][dng] = ngcat['exp_t_s2n'][ng]
    data['sens_1'][dng] = ngcat['exp_e_sens_1'][ng]
    data['sens_2'][dng] = ngcat['exp_e_sens_2'][ng]
    data['stamp_size'][dng] = ngcat['box_size'][ng]

    # Combine the flags we have so far:
    print 'flags range from %d to %d'%(
        np.min(ngcat['flags'][ngcat['flags']>0]),
        np.max(ngcat['flags']))
    print 'exp_flags range from %d to %d'%(
        np.min(ngcat['exp_flags'][ngcat['exp_flags']>0]),
        np.max(ngcat['exp_flags']))
    data['error_flag'] = 2**30
    data['error_flag'][dng] = ngcat['flags'][ng]
    data['error_flag'][dng] |= ngcat['exp_flags'][ng]

    # Calculate mean sensitivity
    data['sens_avg'][dng] = (ngcat['exp_e_sens_1'][ng] + ngcat['exp_e_sens_2'][ng]) / 2.

    # Calculate the recommended weight.
    data['w'][dng] = 1.0/(2.0*0.22*0.22 + ngcat['exp_e_cov_1_1'][ng] + ngcat['exp_e_cov_2_2'][ng])

    # Calculate log10(sb)
    data['log10_sb_i'][dng] = np.log10(np.abs(ngcat['exp_flux_i'][ng]/ngcat['exp_t'][ng]))

    # swap e1 signs
    for tag in ['e_1', 'e_cov_1_2', 'e_cov_2_1']:
        data[tag][dng] *= -1.0
    data['e_cov_2_1'][dng] = data['e_cov_1_2'][dng]

    # unblind
    for tag in ['e_1','e_2']:
        data[tag][dng] /= unblind.get_factor()

    # Bring in the round columns from the s2n catalog
    scat = fitsio.read(config['ngmix_s2n'])
    sc, dsc = numpy_util.match(scat['id'], info['COADD_OBJECTS_ID'])
    data['snr_r'][dsc] = scat['exp_s2n_r'][sc]
    data['t_r'][dsc] = scat['exp_T_r'][sc]
    data['snr_t_r'][dsc] = scat['exp_T_s2n_r'][sc]

    # Check that round_flags are consisten with current error_flag
    print 'round_flags range from %d to %d'%(
        np.min(scat['round_flags'][scat['round_flags']>0]),
        np.max(scat['round_flags']))
    print 'exp_round_flags range from %d to %d'%(
        np.min(scat['exp_round_flags'][scat['exp_round_flags']>0]),
        np.max(scat['exp_round_flags'][scat['exp_round_flags']<2**30]))
    assert np.all(data['error_flag'][dsc][ scat['round_flags'][sc] > 0 ] > 0)
    assert np.all(data['error_flag'][dsc][ scat['exp_round_flags'][sc] == 2**30 ] > 0)

    # Combine the round flags into the error_flag column
    data['error_flag'][dsc] |= scat['round_flags'][sc]
    data['error_flag'][dsc] |= scat['exp_round_flags'][sc]

    # Bring in the mean psf information from the psf catalog
    pcat = fitsio.read(config['ngmix_psfs'])
    pc, dpc = numpy_util.match(pcat['id'], info['COADD_OBJECTS_ID'])
    data['mean_psf_e1'][dpc] = pcat['psfrec_e'][pc,0]
    data['mean_psf_e2'][dpc] = pcat['psfrec_e'][pc,1]
    data['mean_psf_t'][dpc] = pcat['psfrec_T'][pc]

    # swap e1 signs
    for tag in ['mean_psf_e1']:
        data[tag][dng] *= -1.0

    # Before writing out, make all the names upper case
    data.dtype.names = tuple([name.upper() for name in data.dtype.names])

    fitsio.write(config['release_ngmix'],data,clobber=True)
 

def make_im3shape(config, info):
    """Make the im3shape catalog file
    """

    # Setup the dtype:
    dt = np.dtype([
        ('coadd_objects_id', 'i8'),
        ('e_1', 'f8'),
        ('e_2', 'f8'),
        ('nbc_m', 'f8'),
        ('nbc_c1', 'f8'),
        ('nbc_c2', 'f8'),
        ('w', 'f8'),
        ('error_flag', 'i4'),
        ('info_flag', 'i4'),
        ('snr_w', 'f4'),
        ('snr_r', 'f4'),
        ('flux_r', 'f4'),
        ('radius', 'f4'),
        ('is_bulge', 'i4'),
        ('mean_rgpp_rp', 'f4'),
        ('mean_psf_e1', 'f8'),
        ('mean_psf_e2', 'f8'),
        ('mean_psf_fwhm', 'f4'),
        ('ra_shift', 'f4'),
        ('dec_shift', 'f4'),
        ('chi2', 'f4'),
        ('likelihood', 'f4'),
        ('stamp_size', 'i4'),
        ('n_exposure', 'i4') ])

    # Read the source catalog
    imcat = fitsio.read(config['im3shape'])

    # im3shape objects are sorted, but has some that are not in info catalog.
    im, dim = numpy_util.match(imcat['coadd_objects_id'], info['COADD_OBJECTS_ID'])
    assert np.all(imcat['coadd_objects_id'] == sorted(imcat['coadd_objects_id']))

    # Set reasonable defaults for float fields:
    data = np.zeros(len(info), dtype=dt)
    for col, t in dt.descr:
        if 'f' in t:
            data[col] = -9999

    # The coadd_objects_id column can match for all rows
    data['coadd_objects_id'] = info['COADD_OBJECTS_ID']

    # Check that the overlap with imcat is correct
    assert np.all(data['coadd_objects_id'][dim] == imcat['coadd_objects_id'][im])

    # Default error_flag is NO_ATTEMPT
    data['error_flag'] = 2**30
    data['info_flag'] = 2**25

    # Copy in the columns from the source catalog that keep the same name:
    for col in ['nbc_m', 'nbc_c1', 'nbc_c2', 'w', 'error_flag', 'info_flag',
                'radius', 'mean_rgpp_rp', 'mean_psf_fwhm', 'likelihood',
                'stamp_size', 'n_exposure']:
        data[col][dim] = imcat[col][im]

    # Some get a new name:
    data['e_1'][dim] = imcat['e1'][im]
    data['e_2'][dim] = imcat['e2'][im]
    data['snr_w'][dim] = imcat['snr'][im]
    data['snr_r'][dim] = imcat['round_snr'][im]
    data['mean_psf_e1'][dim] = imcat['mean_psf_e1_sky'][im]
    data['mean_psf_e2'][dim] = imcat['mean_psf_e2_sky'][im]
    data['ra_shift'][dim] = imcat['ra_as'][im]
    data['dec_shift'][dim] = imcat['dec_as'][im]
    data['chi2'][dim] = imcat['chi2_pixel'][im]

    # Do a calculation to get the flux from separate bulge/disc fluxes.
    data['is_bulge'][dim] = imcat['bulge_flux'][im] > 0.
    # Only one is non-zero:
    assert np.all((imcat['bulge_flux'] > 0) != (imcat['disc_flux'] > 0))
    data['flux_r'][dim] = imcat['mean_flux'][im] * (imcat['bulge_flux'][im]+imcat['disc_flux'][im])

    # clip the weights
    data['w'][dim] = np.clip(data['w'][dim], 0.0, 0.24**(-2.0))

    # unblind
    for tag in ['e_1','e_2']:
        data[tag][dim] /= unblind.get_factor()
    
    # Before writing out, make all the names upper case
    data.dtype.names = tuple([name.upper() for name in data.dtype.names])

    fitsio.write(config['release_im3shape'],data,clobber=True)

def verify_info(config):
    """Verify that the info file we made is consistent with what Matt's script builds.
    """
    print 'Verifying info catalog...'
    rel = fitsio.read(config['release_info'])
    print 'Full catalog has %d rows'%len(rel)
    v18 = fitsio.read(os.path.join('../v18',config['flatcats_info']))
    print 'v18 catalog has %d rows'%len(v18)
    q = np.where(rel['SVA1_FLAG'] <= 3)[0]
    print 'mask has %d rows'%len(q)
    assert len(q) == len(v18)

    for col in ['COADD_OBJECTS_ID', 'RA', 'DEC',
                'MAG_AUTO_G', 'MAG_AUTO_R', 'MAG_AUTO_I', 'MAG_AUTO_Z',
                'PHOTOZ_BIN', 'MEAN_PHOTOZ']:
        print 'Test rel[%r] == v18[%r]'%(col,col.lower())
        assert np.all(rel[col][q] == v18[col.lower()])

    # v18 used 4 for im3shape, ngmix flags rather than 1.
    v18['im3shape_flags'] /= 4
    v18['ngmix_flags'] /= 4

    for col1, col2 in [ ('IM3SHAPE_FLAG', 'im3shape_flags'),
                        ('NGMIX_FLAG', 'ngmix_flags'),
                        ('SVA1_FLAG', 'sva1_gold_flags') ]:
        print 'Test rel[%r] == v18[%r]'%(col1,col2)
        assert np.all(rel[col1][q] == v18[col2])

    print 'info file passed verification tests'
    return q


def verify_ngmix(config, q, ng):
    """Verify that the ngmix file we made is consistent with what Matt's script builds.
    """
    print 'Verifying ngmix catalog...'
    rel = fitsio.read(config['release_ngmix'])
    print 'Full catalog has %d rows'%len(rel)
    v18 = fitsio.read(os.path.join('../v18',config['flatcats_ngmix']))
    print 'v18 catalog has %d rows'%len(v18)
    print 'mask has %d rows'%len(q)
    assert len(q) == len(v18)

    relq = rel[q]
    for col1, col2 in [ ('COADD_OBJECTS_ID', 'coadd_objects_id'),
                        ('E_1', 'exp_e_1'), 
                        ('E_2', 'exp_e_2'),
                        ('SENS_AVG', 'exp_e_sens_avg'),
                        ('W', 'exp_w'),
                        ('E_COV_1_1', 'exp_e_cov_1_1'),
                        ('E_COV_1_2', 'exp_e_cov_1_2'),
                        ('E_COV_2_1', 'exp_e_cov_2_1'),
                        ('E_COV_2_2', 'exp_e_cov_2_2'),
                        ('MEAN_PSF_E1', 'psfrec_e_1'),
                        ('MEAN_PSF_E2', 'psfrec_e_2'),
                        ('STAMP_SIZE', 'box_size') ]:
        print 'Test rel[%r] == v18[%r]'%(col1,col2)
        assert np.all(relq[col1] == v18[col2])

    for col1, col2 in [ ('SNR_W', 'exp_s2n_w'),
                        ('SNR_R', 'exp_s2n_r'),
                        ('FLUX_I', 'exp_flux_i'),
                        ('MAG_I', 'exp_mag_i'),
                        ('T', 'exp_T'),
                        ('T_ERR', 'exp_T_err'),
                        ('T_R', 'exp_T_r'),
                        ('SNR_T', 'exp_T_s2n'),
                        ('SNR_T_R', 'exp_T_s2n_r'),
                        ('LOG10_SB_I', 'exp_log10sb_i'),
                        ('MEAN_PSF_T', 'psfrec_T'),
                        ('SENS_1', 'exp_e_sens_1'),
                        ('SENS_2', 'exp_e_sens_2'),
                        ('ARATE', 'exp_arate'),
                        ('MASK_FRAC', 'mask_frac') ]:
        print 'Test rel[%r] ~= v18[%r]'%(col1,col2)
        assert np.all(np.abs(relq[col1] - v18[col2]) <= 1.e-7 * np.abs(v18[col2]))

    assert np.all(relq['ERROR_FLAG'] == (v18['flags'] | v18['exp_flags']))
    
    print 'ngmix file passed verification tests'


def verify_im3shape(config, q, im):
    """Verify that the im3shape file we made is consistent with what Matt's script builds.
    """
    print 'Verifying im3shape catalog...'
    rel = fitsio.read(config['release_im3shape'])
    print 'Full catalog has %d rows'%len(rel)
    v18 = fitsio.read(os.path.join('../v18',config['flatcats_im3shape']))
    print 'v18 catalog has %d rows'%len(v18)
    print 'mask has %d rows'%len(q)
    assert len(q) == len(v18)

    relq = rel[q][im]
    v18 = v18[im]
    for col1, col2 in [ ('COADD_OBJECTS_ID', 'coadd_objects_id'),
                        ('E_1', 'e1'), 
                        ('E_2', 'e2'),
                        ('NBC_M', 'nbc_m'),
                        ('NBC_C1', 'nbc_c1'),
                        ('NBC_C2', 'nbc_c2'),
                        ('W', 'w'),
                        ('ERROR_FLAG', 'error_flag'),
                        ('INFO_FLAG', 'info_flag'),
                        ('MEAN_PSF_E1', 'mean_psf_e1_sky'),
                        ('MEAN_PSF_E2', 'mean_psf_e2_sky'),
                        ('STAMP_SIZE', 'stamp_size'),
                        ('N_EXPOSURE', 'n_exposure') ]:
        print 'Test rel[%r] == v18[%r]'%(col1,col2)
        assert np.all(relq[col1] == v18[col2])

    for col1, col2 in [ ('SNR_W', 'snr'),
                        ('SNR_R', 'round_snr'),
                        ('RADIUS', 'radius'),
                        ('MEAN_RGPP_RP', 'mean_rgpp_rp'),
                        ('MEAN_PSF_FWHM', 'mean_psf_fwhm'),
                        ('RA_SHIFT', 'ra_as'),
                        ('DEC_SHIFT', 'dec_as'),
                        ('CHI2', 'chi2_pixel'),
                        ('LIKELIHOOD', 'likelihood') ]:
        print 'Test rel[%r] ~= v18[%r]'%(col1,col2)
        assert np.all(np.abs(relq[col1] - v18[col2]) <= 1.e-7 * np.abs(v18[col2]))

    assert np.all(relq['IS_BULGE'] == (v18['bulge_flux'] > 0.))
    assert np.all(relq['IS_BULGE'] == (v18['disc_flux'] == 0.))
    flux1 = relq['FLUX_R']
    flux2 = v18['mean_flux'] * (v18['bulge_flux'] + v18['disc_flux'])
    assert np.all(np.abs(flux1 - flux2) <= 1.e-7 * np.abs(flux2))
    
    print 'im3shape file passed verification tests'

if __name__ == "__main__":
    import yaml
    with open(sys.argv[1],'r') as fp:
        config = yaml.load(fp)        
    info = make_info(config)
    make_ngmix(config, info)
    make_im3shape(config, info)
    q = verify_info(config)
    verify_ngmix(config, q, info['NGMIX_FLAG'][q]==0)
    verify_im3shape(config, q, info['IM3SHAPE_FLAG'][q]==0)

