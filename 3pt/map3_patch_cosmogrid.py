import treecorr
import time
import numpy as np
import pickle as pickle
import yaml
import os
import sys
import healpy as hp
import argparse
# File name: run_with_mpi.py
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def load_cosmogrid_catalog(fname, zbin, randomize_radec=False):
    """
    zbins: redshift bin, 1,2,3,4

    tocat: whether convert to treecorr catalog
    """
    # data file
    data = np.load(fname)
    inds = np.load('/global/cfs/cdirs/des/ssunao/3pcf/measurement/Y3Maps_NSIDE1024_Cov_V0/Indices.npy')

    # ra, dec
    nside = 1024
    pix  = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, pix)
    ra = phi
    dec = np.pi/2 - theta
    ra, dec = ra[inds], dec[inds]
    unit = 'rad'

    # shear
    g1 = data[zbin-1, 0]
    g2 =-data[zbin-1, 1]
    
    if randomize_radec:
        ## new feature to test the coriance inaccuracy due to the resolution.
        # We add random disturb in ra, dec as large as size of resolution (healpix cell)
        rng = np.random.default_rng(0)
        hp_resol = hp.nside2resol(nside) # radian
        print('Adding random disturb in ra and dec... : {} radian'.format(hp_resol))
        r   = hp_resol * rng.uniform(0, 1, ra.size)**0.5 # P(r) \propto r in polar coordinate
        phi = rng.uniform(0,2*np.pi, ra.size)
        ra += r * np.cos(phi)
        dec+= r * np.sin(phi)
    
    # weight
    wt = np.load('/global/cfs/cdirs/des/ssunao/3pcf/measurement/weight_map/weiht_map_bin%d.npy'%zbin)[inds]

    return ra, dec, g1, g2, wt, unit


def measure_3pcf_map3_auto_patch(zbin, run, ncut, nreal, method='jackknife', npatch=100, randomize_radec=False):
    # read the data
    fname = '/global/cfs/cdirs/des/ssunao/3pcf/measurement/Y3Maps_NSIDE1024_Cov_20240627/ShearMaps_DESY3_run_{:04d}_ncut{:d}_nreal{:d}.npy'
    fname = fname.format(run, ncut, nreal)
    ra, dec, g1, g2, wt, unit = load_cosmogrid_catalog(fname, zbin, randomize_radec=randomize_radec)
    
    # convert it to treecorr catalog
    print('Making TreeCorr catalog...')
    cat = treecorr.Catalog(ra=ra, dec=dec, 
                           g1=g1, g2=g2, 
                           ra_units=unit, dec_units=unit, 
                           w = wt,
                           npatch=npatch)
    
    # make processor
    print('Processing...')
    maxn = 100
    for_ggg = dict(min_sep=0.5, max_sep=80, 
                   sep_units='arcmin', nbins=20, 
                   max_n=maxn, bin_type='LogMultipole', verbose=2)
    three_pt_corr = treecorr.GGGCorrelation(for_ggg)
    
    # measure 3pcf
    print('Measuring 3PCF...')
    three_pt_corr.process(cat, comm=comm)
    
    # jackknife covariance of map3
    print('{} map3 covariance...'.format(method))
    r2 = np.array([7,14,25,40])
    func = lambda corr: corr.toSAS(phi_bin_size=0.05).calculateMap3(R=r2)[0]
    covMap3 = three_pt_corr.estimate_cov(method=method, func=func, comm=comm)
    print(covMap3)
    
    # save covariance
    if rank == 0:
        method_short = {'jackknife':'JK', 'sample':'SAMP', 'bootstrap': 'BOOTSTRAP'}[method]
        if randomize_radec: 
            method_short = method_short+'_RANDRADEC'
        filename = 'covariance/covMap3_{}_CosmoGrid_zbin{}_run_{:04d}_ncut{:d}_nreal{:d}_npatch{}.dat'
        filename = filename.format(method_short, zbin, run, ncut, nreal, npatch)
        print('Writing to {}...'.format(filename))
        np.savetxt(filename, covMap3)
    

def main():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Run 3-point correlation function measurement with given parameters.")

    # Add arguments
    parser.add_argument("--zbin", type=int, required=True, help="Redshift bin number")
    parser.add_argument("--run", type=int, default=0, help="Run number (default: 0)")
    parser.add_argument("--ncut", type=int, default=0, help="Cut number (default: 0)")
    parser.add_argument("--nreal", type=int, default=0, help="Number of realizations (default: 0)")
    parser.add_argument("--method", type=str, default='jackknife', help="Method to use (default: 'jackknife')")
    parser.add_argument("--npatch", type=int, default=100, help="Number of patches (default: 100)")
    parser.add_argument("--randomize_radec", help="switch to randomize ra, dec", action='store_true')
    

    # Parse arguments
    args = parser.parse_args()

    # Call your function with parsed arguments
    measure_3pcf_map3_auto_patch(
        zbin=args.zbin,
        run=args.run,
        ncut=args.ncut,
        nreal=args.nreal,
        method=args.method,
        npatch=args.npatch,
        randomize_radec=args.randomize_radec
    )

if __name__ == '__main__':
    main()
