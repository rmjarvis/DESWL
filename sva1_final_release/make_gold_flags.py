#!/usr/bin/env python
import os
import sys
import numpy as np
import fitsio
import numpy_util

def make_gold_flags(config):
    """
    code to make gold flags for DES WL analysis

    People to bug if things are weird
    =================================
    
    Matthew Becker
    Eli Rykoff
    Erin Sheldon
    
    Input Data
    ==========
    We use the following sources of data.
    
    1) SVA1 gold v1.0.4 with the additional star masking, 
    bad regions mask, crazy colors, and large g- vs i-band 
    offsets. See this page
    
    https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/SVA1_Gold_Catalog_v10#Gold-104
    
    2) ngmix11 - used to cut low surface brightness junk, 
    additonal stars and to set the ngmix flags of course
    
    See these pages
    
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Ngmix011
    https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/A_Modest_Proposal_for_Preliminary_StarGalaxy_Separation
    

    Flag File Description
    ======================
    The flag file consists of three columns which are
    
    sva1_gold_flags
    ---------------
    
    These are flags related to SVA1 image processing and galaxy selection. We 
    have attempted to preserve as much information as possible. The field is 
    created by concatenating the SExtractor flags (i-band), SVA1 gold flags, 
    and additional galaxy selection criteria from ngmix11
    
    bit 2^bit   name              notes
    ----------------------------------------------------------------------
    0   1       SX_KINDOF_BLEND    bit 0 from the i-band sexctractor flags
    1   2       SX_DEF_BLEND       bit 1 from the i-band sexctractor flags
    2   4       MODEST_STAR        SVA1 Gold MODEST_CLASS == 2
    3   8       MODEST_JUNK        SVA1 Gold MODEST_CLASS == 0
    4   16      CRAZY_COLORS       SVA1 Gold crazy mag auto colors
    5   32      BADMASK04          in region with lots of JUNK_DPOS objects (removes ~4% of area)
    6   64      NEAR_2MASS         near a 2MASS star
    7   128     JUNK_DPOS          large offset in g and i band windowed positions
    8   256     NOT_IN_NGMIX011    object did not get measured by ngmix011
    9   512     NGMIX011_STAR      object is a star according to ngmix011 modest classifier 
    10  1024    LOW_SB_NGMIX011    object has low surface brightness in ngmix011
    11  2048    BAD_MEAS_NGMIX011  object does not satisfy good measurement cuts in ngmix
    
    Some cuts from ngmix011 are applied to the full gold sample. These are detailed below 
    and are in flag bits 8 through 12.
    
    We use a modest-style star galaxy separation such that objects 
    with 
    
        exp_t + exp_t_err > 0.02
    
    are galaxies. We also cut objects with 

        exp_mag_i + 3.5*log10(exp_flux_i/exp_t) < 28.0

    where

        exp_flux_i = 10.0^((exp_mag_i-30.0)/(-2.5))

    Note that for the ngmix cuts we require that for ngmix

        flags == 0
        exp_flags == 0
        0.4 < exp_arate < 0.6

    so that object has a good measurement at all. 

    sva1_spte_flags
    ---------------
    
    We flag all galaxies with ra,dec not in 

        50 < ra < 100
        dec < -35.0

    as outside of the SPTE region.

    sva1_gold_mag_flags
    -------------------

    We flag all objects that do not have griz mags.
    
    bit 2^bit   name              notes
    ----------------------------------------------------------------------
    0   1       MISSING_GRIZ_MAGS object does not have all griz mags
    
    """

    # first read in the data
    gold = fitsio.read(config['gold'],lower=True,columns=['coadd_objects_id','modest_class','flags_i','ra','dec'])
    auto = fitsio.read(config['auto'],lower=True)
    badflag = fitsio.read(config['badflag'],lower=True)
    ng = fitsio.read(config['ngmix'])
    
    # make output data
    dflag = np.zeros(len(gold),dtype=[('coadd_objects_id','i8'),('ra','f8'),('dec','f8'),
                                      ('sva1_gold_flags','i8'),('sva1_spte_flags','i8'),
                                      ('sva1_gold_mag_flags','i8'),
                                      ('mag_auto_g','f8'),('mag_auto_r','f8'),
                                      ('mag_auto_i','f8'),('mag_auto_z','f8')])
    for tag in ['coadd_objects_id','ra','dec']:
        dflag[tag] = gold[tag]
    for tag in ['mag_auto_g','mag_auto_r','mag_auto_i','mag_auto_z']:
        dflag[tag] = auto[tag]
    
    # set sva1_gold_flags
    tag = 'sva1_gold_flags'
    dflag[tag] = gold['flags_i']
    
    q, = np.where(gold['modest_class'] == 2)
    if len(q) > 0:
        dflag[tag][q] |= 4
        
    q, = np.where(gold['modest_class'] == 0)
    if len(q) > 0:
        dflag[tag][q] |= 8
    
    # crazy colors
    # (g-r) < -1, (g-r) > 4, (i-z) < -1, (i-z) > 4
    q, = np.where((auto['mag_auto_g'] - auto['mag_auto_r'] < -1.0) |
                  (auto['mag_auto_g'] - auto['mag_auto_r'] > 4.0) |
                  (auto['mag_auto_i'] - auto['mag_auto_z'] < -1.0) |
                  (auto['mag_auto_i'] - auto['mag_auto_z'] > 4.0))
    if len(q) > 0:
        dflag[tag][q] |= 16
         
    # badflags masking
    q, = np.where((badflag['badflag'] & 2) != 0)
    if len(q) > 0:
        dflag[tag][q] |= 32
        
    q, = np.where((badflag['badflag'] & 4) != 0)
    if len(q) > 0:
        dflag[tag][q] |= 64
        
    q, = np.where((badflag['badflag'] & 8) != 0)
    if len(q) > 0:
        dflag[tag][q] |= 128
        
    del badflag
        
    # ngmix cuts
    # first have to match
    gids = gold['coadd_objects_id']
    ngids = ng['coadd_objects_id']
    ginds,nginds = numpy_util.match(gids,ngids)
    
    # flag stuff not in ngmix011
    sginds = set(ginds)
    totginds = set(range(len(gold)))
    diffginds = totginds - sginds
    q = np.fromiter(diffginds,dtype=int,count=len(diffginds))
    if len(q) > 0:
        dflag[tag][q] |= 256
                   
    # bad measurements
    q, = np.where(~((ng['flags'][nginds] == 0) &
                    (ng['exp_flags'][nginds] == 0) &
                    (ng['exp_arate'][nginds] > 0.4) & 
                    (ng['exp_arate'][nginds] < 0.6)))
    if len(q) > 0:
        dflag[tag][ginds[q]] |= 2048

    # low SB
    sb = ng['exp_mag_i'] + 3.5*np.log10(np.abs(ng['exp_flux_i']/ng['exp_t']))
    q, = np.where(~((ng['flags'][nginds] == 0) &
                    (ng['exp_flags'][nginds] == 0) &
                    (ng['exp_arate'][nginds] > 0.4) & 
                    (ng['exp_arate'][nginds] < 0.6) & 
                    (sb[nginds] > 28.0)))
    if len(q) > 0:
        dflag[tag][ginds[q]] |= 1024
        
    # ngmix star
    nsig = 1.0
    sgc = ng['exp_t'] + nsig*ng['exp_t_err'] 
    q, = np.where(~((ng['flags'][nginds] == 0) &
                    (ng['exp_flags'][nginds] == 0) &
                    (ng['exp_arate'][nginds] > 0.4) & 
                    (ng['exp_arate'][nginds] < 0.6) & 
                    (sgc[nginds] > 0.02)))
    if len(q) > 0:
        dflag[tag][ginds[q]] |= 512
        
    # non-99 mag_auto in all bands
    tag = 'sva1_gold_mag_flags'
    q, = np.where(~((np.abs(auto['mag_auto_g']) < 99.0) &
                    (np.abs(auto['mag_auto_r']) < 99.0) &
                    (np.abs(auto['mag_auto_i']) < 99.0) &
                    (np.abs(auto['mag_auto_z']) < 99.0)))
    if len(q) > 0:
        dflag[tag][q] |= 1
    del auto
    
    #set sva1_spte_flags
    q, = np.where(~((gold['ra'] > 50.0) & (gold['ra'] < 100.0) & (gold['dec'] < -35.0)))
    if len(q) > 0:
        dflag['sva1_spte_flags'][q] != 1
        
    #write to disk
    fitsio.write(config['gold_flags'],dflag,clobber=True)
    
if __name__ == "__main__":
    import yaml
    with open(sys.argv[1],'r') as fp:
        config = yaml.load(fp)
    make_gold_flags(config)

