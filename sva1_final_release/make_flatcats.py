#!/usr/bin/env python
"""
This is code to make the flat catalogs for the final DES SVA1 WL analysis.

People to bug if things are weird
=================================

Matthew Becker 
Eli Rykoff (SVA1 Gold guru)
Michael Troxel (for photo-z binning)
im3shape team (Troxel, Tomek, Joe, Sarah, Niall, ...)
ngmix team (Erin Sheldon)

Input Data
==========
I use the WL flags files along with the ngmix and im3shape outputs.
I also use the photo-z binning from Troxel and the NBC from Tomek.

Usgae and Outputs
=================
The selection in these files is 

sva1_spte_flags == 0 &
sva1_gold_flags <= 3 &
sva1_gold_mag_flags == 0 &
(im3shape_flags == 0 | ngmix_flags == 0)

The three sva1 flags fields get combined into a single flag field, sva1_flags, 
with the first two bits of the sva1_gold_flags, 
followed by a third bit set if any higher bits are nonzero, 
followed by two bits set to 0 or 1 if sva1_spte_flags == 0 and sva1_gold_mag_flags == 0, in that order.

To get the final ngmix or im3shape catalogs, you then make your own selection. The 
ngmix_flags and im3shape_flags have less conservative shear source selections built in. 

Three output files are made which can all be read in parallel

1) *_info.py - catalog with all flags, ra, dec, photo-z bins, etc.
2) *_im3shape.py - catalog with common im3shape parameters
3) *_ngmix.py - catalog with common ngmix011 parameters

"""

import os
import sys
import numpy as np
import fitsio
import numpy_util
import unblind

def main(config):
    ###############################
    # info file

    # do the flags
    flags = fitsio.read(config['gold_flags'])
    imflags = fitsio.read(config['im3shape_flags'])
    ngflags = fitsio.read(config['ngmix_flags'])    
    if False:
        q, = np.where((flags['sva1_gold_flags'] <= 3) &
                      (flags['sva1_spte_flags'] == 0) &
                      (flags['sva1_gold_mag_flags'] == 0) &
                      ((imflags['im3shape_flags'] == 0) | 
                       (ngflags['ngmix_flags'] == 0)))
    else:
        q, = np.where((flags['sva1_gold_flags'] <= 3) &
                      (flags['sva1_spte_flags'] == 0) &
                      (flags['sva1_gold_mag_flags'] == 0))
    flags = flags[q]
    imflags = imflags[q]
    ngflags = ngflags[q]
    
    # print number of sources
    for code,cflags,tag in zip(['ngmix','im3shape'],[ngflags,imflags],['ngmix_flags','im3shape_flags']):
        q, = np.where((flags['sva1_gold_flags'] == 0) & (cflags[tag] == 0))
        print "%s:" % code,len(q)
    
    # read in the photo-z bins - match
    dtot = []
    dz = []
    if 'photoz_binning_im3shape' in config:
        dim = np.load(config['photoz_binning_im3shape']).astype(int)        
        for i in xrange(len(dim)):
            dtot.append(tuple(dim[i]))
        dim = np.load(config['photoz_binning_im3shape']).astype(float)
        for i in xrange(len(dim)):
            dz.append(tuple(dim[i])[-1])        
        del dim
    if 'photoz_binning_ngmix' in config:
        dng = np.load(config['photoz_binning_ngmix']).astype(int)
        for i in xrange(len(dng)):
            dtot.append(tuple(dng[i]))    
        dng = np.load(config['photoz_binning_ngmix']).astype(float)
        for i in xrange(len(dng)):
            dz.append(tuple(dng[i])[-1])        
        del dng
    if 'photoz_binning_extra' in config:
        de = np.load(config['photoz_binning_extra']).astype(int)
        for i in xrange(len(de)):
            dtot.append(tuple(de[i]))
        de = np.load(config['photoz_binning_extra']).astype(float)
        for i in xrange(len(de)):
            dz.append(tuple(de[i])[-1])        
        del de

    if len(dtot) > 0:
        dtot = np.array(dtot,dtype=[('coadd_objects_id','i8'),('photoz_bin','i4'),('mean_photoz','f8')])    
        dz = np.array(dz,dtype='f8')
        dtot['mean_photoz'][:] = dz[:]
        u, indices = np.unique(dtot['coadd_objects_id'], return_index=True)
        dtot = dtot[indices]
        del u,indices

        # make sure all in file
        fids = set(flags['coadd_objects_id'])
        pzids = set(dtot['coadd_objects_id'])
        if not fids <= pzids:
            print "not all gals got photoz binning indexes!"
            if True:
                fids = flags['coadd_objects_id']
                pzids = dtot['coadd_objects_id']    
                finds,pzinds = numpy_util.match(fids,pzids)
                finds = set(finds)
                tinds = set(np.arange(len(flags)))
                minds = tinds - finds
                minds = np.array(list(minds))
                #print set(flags['coadd_objects_id'][minds])
                print '# of missing gals is',len(minds)
                #assert False

        fids = flags['coadd_objects_id']
        pzids = dtot['coadd_objects_id']    
        finds,pzinds = numpy_util.match(fids,pzids)
    
    # make final catalog
    ftags = ['mag_auto_g','mag_auto_r','mag_auto_i','mag_auto_z',
             'ra','dec','mean_photoz']
    itags = ['coadd_objects_id','sva1_flags','im3shape_flags','ngmix_flags',
             'sva1_gold_flags','sva1_spte_flags','sva1_gold_mag_flags',
             'photoz_bin']
    dlist = []
    for ftag in ftags:
        dlist.append((ftag,'f8'))
    for itag in itags:
        dlist.append((itag,'i8'))
    d = np.zeros(len(flags),dtype=dlist)
    d['photoz_bin'][:] = -99
    d['mean_photoz'][:] = -99.0
    for tag in flags.dtype.names:
        d[tag] = flags[tag]
    for tag in imflags.dtype.names:
        d[tag] = imflags[tag]
    for tag in ngflags.dtype.names:
        d[tag] = ngflags[tag]
    if len(dtot) > 0:
        d['photoz_bin'][finds] = dtot['photoz_bin'][pzinds]    
        d['mean_photoz'][finds] = dtot['mean_photoz'][pzinds]    
            
    q, = np.where(flags['sva1_gold_flags'] & 1)
    d['sva1_flags'][q] |= 1
    
    q, = np.where(flags['sva1_gold_flags'] & 2)
    d['sva1_flags'][q] |= 2

    q, = np.where((flags['sva1_gold_flags'] >> 2) != 0)
    d['sva1_flags'][q] |= 4
    
    q, = np.where(flags['sva1_spte_flags'] != 0)
    d['sva1_flags'][q] |= 8
    
    q, = np.where(flags['sva1_gold_mag_flags'] != 0)
    d['sva1_flags'][q] |= 16
    
    fitsio.write(config['flatcats_info'],d,clobber=True)
    del flags
    del imflags
    del ngflags
    del dtot
    
    ###############################
    # ngmix file
    # first have to match to psf and s2n
    ng = fitsio.read(config['ngmix'])
    psf = fitsio.read(config['ngmix_psfs'])
    pids = psf['id']
    ngids = ng['coadd_objects_id']
    pinds,nginds = numpy_util.match(pids,ngids)
    assert np.all(psf['flags'][pinds] == 0)    
    ng['psfrec_t'][nginds] = psf['psfrec_T'][pinds]
    ng['psfrec_e_1'][nginds] = psf['psfrec_e'][pinds,0]
    ng['psfrec_e_2'][nginds] = psf['psfrec_e'][pinds,1]
    ng = ng[nginds]    
    del psf,pids,ngids,pinds,nginds
    
    # now match to new s2n measures
    ngs2n = fitsio.read(config['ngmix_s2n'])
    s2nids = ngs2n['id']
    ngids = ng['coadd_objects_id']
    s2ninds,nginds = numpy_util.match(s2nids,ngids)
    ng = ng[nginds]
    ngs2n = ngs2n[s2ninds]
    del s2ninds,nginds,s2nids,ngids
    
    # make sure everything is in the ngmix file
    q, = np.where(d['ngmix_flags'] == 0)
    flagged_ids = set(d['coadd_objects_id'][q])
    ngids = set(ng['coadd_objects_id'])
    assert flagged_ids <= ngids    
    
    # match to total file
    fids = d['coadd_objects_id']
    ngids = ng['coadd_objects_id']
    finds,nginds = numpy_util.match(fids,ngids)
    assert len(finds) == len(d)
    
    # extra tags
    dlist = ng.dtype.descr
    dlist.append(('exp_log10sb_i','f8'))
    dlist.append(('exp_w','f8'))
    dlist.append(('exp_e_sens_avg','f8'))
    dlist.append(('exp_T_r','f8'))
    dlist.append(('exp_s2n_r','f8'))
    dlist.append(('exp_T_s2n_r','f8'))
    dnew = []
    for dt in dlist:
        if '_t' in dt[0]:
            name = dt[0]
            name = name.replace('_t','_T')
            dnew.append((name,dt[1]))
        elif dt[0] not in ['exp_e_1','exp_e_2','exp_e_cov_1_1','exp_e_cov_1_2','exp_e_cov_2_2']:
            dnew.append(dt)
    dlist = dnew
    for tag in ['exp_e_1','exp_e_2','exp_e_cov_1_1','exp_e_cov_1_2','exp_e_cov_2_1','exp_e_cov_2_2']:
        dlist.append((tag,'f8'))

    ng_final = np.zeros(len(d),dtype=dlist)
    ng_final['coadd_objects_id'] = d['coadd_objects_id']
    for tag in ng.dtype.names:
        if tag != 'coadd_objects_id':
            if '_t' in tag:
                final_tag = tag.replace('_t','_T')
            else:
                final_tag = tag
            ng_final[final_tag][finds] = ng[tag][nginds]
        else:
            assert np.array_equal(ng_final[tag][finds],ng[tag][nginds])
    del ng
    
    ng_final['exp_log10sb_i'][finds] = np.log10(np.abs(ng_final['exp_flux_i'][finds]/ng_final['exp_T'][finds]))
    ng_final['exp_e_sens_avg'] = (ng_final['exp_e_sens_1'] + ng_final['exp_e_sens_2'])/2.0
    
    # s2n measures
    for tag in ['exp_T_r','exp_s2n_r','exp_T_s2n_r']:
        ng_final[tag][finds] = ngs2n[tag][nginds]
    del ngs2n
    
    # swap e1 signs
    for tag in ['exp_e_1','exp_e_cov_1_2','psfrec_e_1']:
        ng_final[tag] *= -1.0
    
    # fianl tag comp that depends on psf sign
    ng_final['exp_e_cov_2_1'] = ng_final['exp_e_cov_1_2']
    ng_final['exp_w'] = 1.0/(2.0*0.22*0.22 + ng_final['exp_e_cov_1_1'] + ng_final['exp_e_cov_2_2'])
    
    # unblind
    for tag in ['exp_e_1','exp_e_2']:
        ng_final[tag] = ng_final[tag]/unblind.get_factor()
    
    fitsio.write(config['flatcats_ngmix'],ng_final,clobber=True)
    del ng_final
    
    #####################
    # im3shape file
    # find stuff in catalog
    im = fitsio.read(config['im3shape'])
    
    # make sure everything is in the im3shape file
    q, = np.where(d['im3shape_flags'] == 0)
    flagged_ids = set(d['coadd_objects_id'][q])
    imids = set(im['coadd_objects_id'])
    assert flagged_ids <= imids
    
    fids = d['coadd_objects_id']
    imids = im['coadd_objects_id']
    finds,iminds = numpy_util.match(fids,imids)
    
    # fill data
    dlist = im.dtype.descr
    im_final = np.zeros(len(d),dtype=dlist)
    im_final['coadd_objects_id'] = d['coadd_objects_id']
    for tag in im.dtype.names:
        if tag != 'coadd_objects_id':
            im_final[tag][finds] = im[tag][iminds]
        else:
            assert np.array_equal(im_final[tag][finds],im[tag][iminds])
        
    # clip the weights
    np.clip(im_final['w'], 0.0, 0.24**(-2.0), im_final['w'])
            
    # unblind
    for tag in ['e1','e2']:
        im_final[tag] = im_final[tag]/unblind.get_factor()
    
    fitsio.write(config['flatcats_im3shape'],im_final,clobber=True)
    
if __name__ == "__main__":
    import yaml
    with open(sys.argv[1],'r') as fp:
        config = yaml.load(fp)        
    main(config)

