#!/usr/bin/env python
import os
import sys
import numpy as np
import fitsio
import numpy_util

def make_ngmix_flags(config):
    """
    code to make ngmix flags

    People to bug if things are weird
    =================================

    Matthew Becker
    Erin Sheldon    
    Mike Jarvis
    
    Input Data
    ==========
    
    ngmix11 - use to set the flags
    
    See these pages
    
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Ngmix011
    https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/A_Modest_Proposal_for_Preliminary_StarGalaxy_Separation
   
    ngmix_flags
    ===========

    We use three levels of flagging, selecting sets of ngmix 
    sources which have various potential levels of additive and 
    multiplicative biases.
    

    bit 2^bit   name              notes
    ----------------------------------------------------------------------
    2   4       NOTORBAD_MEAS     object is not in ngmix or does not meet required flags belo
    
    Required Cuts
    --------------
    Any source which does not meet these cuts is not well measured in ngmix
    
    flags == 0
    exp_flags == 0
    0.4 < exp_arate < 0.6
    exp_e_sens_1 > 0.0
    exp_e_sens_2 > 0.0
    round_flags == 0
    exp_round_flags == 0
    exp_s2n_r > 15
    exp_T_r/psfrec_T > 0.15 
    
    """        
    
    # first read in the data
    gold = fitsio.read(config['gold'],lower=True,columns=['coadd_objects_id'])
    ng = fitsio.read(config['ngmix'])
    psf = fitsio.read(config['ngmix_psfs'])
    ngs2n = fitsio.read(config['ngmix_s2n'])
    
    # match to psf
    pids = psf['id']
    ngids = ng['coadd_objects_id']
    pinds,nginds = numpy_util.match(pids,ngids)
    assert np.all(psf['flags'][pinds] == 0)    
    ng = ng[nginds]
    psf = psf[pinds]
    del pids,pinds
    
    # match to s/n
    s2nids = ngs2n['id']
    ngids = ng['coadd_objects_id']
    s2ninds,nginds = numpy_util.match(s2nids,ngids)
    ng = ng[nginds]
    ngs2n = ngs2n[s2ninds]
    del s2ninds,s2nids
    
    # now do output, matching to gold
    dflag = np.zeros(len(gold),dtype=[('ngmix_flags','i8')])
    
    tag = 'ngmix_flags'
    gids = gold['coadd_objects_id']
    ngids = ng['coadd_objects_id']
    ginds,nginds = numpy_util.match(gids,ngids)
    
    # flag stuff not in ngmix011
    sginds = set(ginds)
    totginds = set(range(len(gold)))
    diffginds = totginds - sginds
    q = np.fromiter(diffginds,dtype=int,count=len(diffginds))
    if len(q) > 0:
        dflag[tag][q] |= 4
        
    # bad measurements
    q, = np.where(~((ng['flags'][nginds] == 0) &
                    (ng['exp_flags'][nginds] == 0) &
                    (ng['exp_arate'][nginds] > 0.4) &
                    (ng['exp_arate'][nginds] < 0.6) &
                    (ng['exp_e_sens_2'][nginds] > 0.0) &
                    (ng['exp_e_sens_1'][nginds] > 0.0) &
                    (ngs2n['round_flags'][nginds] == 0) & 
                    (ngs2n['exp_round_flags'][nginds] == 0) &
                    (ngs2n['exp_s2n_r'][nginds] > 15.0) & 
                    (ngs2n['exp_T_r'][nginds]/psf['psfrec_T'][nginds] > 0.15)                    
                    ))
    if len(q) > 0:
        dflag[tag][ginds[q]] |= 4
    
    #write to disk
    fitsio.write(config['ngmix_flags'],dflag,clobber=True)
    
if __name__ == "__main__":
    import yaml
    with open(sys.argv[1],'r') as fp:
        config = yaml.load(fp)
    make_ngmix_flags(config)

