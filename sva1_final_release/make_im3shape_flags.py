#!/usr/bin/env python
import os
import sys
import numpy as np
import fitsio
import numpy_util

def make_im3shape_flags(config):
    """
    code to make im3shape flags

    People to bug if things are weird
    =================================

    Matthew Becker
    im3shape team
    
    Input Data
    ==========
    
    im3shape v97 - use to set the flags
    
    im3shape_flags
    ==============
    
    bit 2^bit   name              notes
    ----------------------------------------------------------------------
    2   4       NOTORBAD_MEAS     object is not in gold or does not meet required flags below
    
    Required Cuts
    --------------
    Any source which does not meet these cuts is not well measured in im3shape
    
    error_flag == 0
    info_flag == 0
    snr > 15    
    mean_rgpp_rp > 1.2
    finite(psf_e1,psf_e1,nbc_m,nbc_c1,nbc_c2)
    """

    # first read in the data
    gold = fitsio.read(config['gold'],lower=True,columns=['coadd_objects_id'])
    im = fitsio.read(config['im3shape'])
    
    dflag = np.zeros(len(gold),dtype=[('im3shape_flags','i8')])
    
    #####################
    # do im3shape cuts
    tag = 'im3shape_flags'
    
    # find stuff in catalog
    gids = gold['coadd_objects_id']
    imids = im['coadd_objects_id']
    ginds,iminds = numpy_util.match(gids,imids)
    
    # flag stuff not in cat
    sginds = set(ginds)
    totginds = set(range(len(gold)))
    diffginds = totginds - sginds
    q = np.fromiter(diffginds,dtype=int,count=len(diffginds))
    if len(q) > 0:
        dflag[tag][q] |= 4
        
    # required cuts
    q, = np.where(~((im['error_flag'][iminds] == 0) & 
                    (im['info_flag'][iminds] == 0) & 
                    (im['snr'][iminds] > 15.0) &
                    (im['mean_rgpp_rp'][iminds] > 1.2) &
                    (np.isfinite(im['mean_psf_e1_sky'][iminds])) &
                    (np.isfinite(im['mean_psf_e2_sky'][iminds])) & 
                    (np.isfinite(im['nbc_m'][iminds])) &
                    (np.isfinite(im['nbc_c1'][iminds])) &
                    (np.isfinite(im['nbc_c2'][iminds]))
                    ))
    if len(q) > 0:
        dflag[tag][ginds[q]] |= 4
        
    #write to disk
    fitsio.write(config['im3shape_flags'],dflag,clobber=True)
    
if __name__ == "__main__":
    import yaml
    with open(sys.argv[1],'r') as fp:
        config = yaml.load(fp)
    make_im3shape_flags(config)

