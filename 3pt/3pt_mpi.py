import treecorr
import glob, os
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

data_dir = '/global/homes/s/seccolf/shear3pt/SHARED_DATA'

zbin1_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin1', '*.fits')))
zbin2_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin2', '*.fits')))
zbin3_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin3', '*.fits')))
zbin4_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin4', '*.fits')))

params_cat = dict(ra_col='ra', dec_col='dec', ra_units='rad', dec_units='rad',
                  w_col='w', g1_col='g1', g2_col='g2', flip_g1=True)

cats_1 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin1_files)]
cats_2 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin2_files)]
cats_3 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin3_files)]
cats_4 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin4_files)]

rmin = 1
rmax = 25
nr = 10

narrow = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
              min_u=0.0, max_u=1, nubins=20,
              min_v=0.0, max_v=0.1, nvbins=1, verbose=2)

wide = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
            min_u=0.9, max_u=1, nubins=1,
            min_v=0.0, max_v=0.8, nvbins=20, verbose=2)

wider = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
             min_u=0.9, max_u=1, nubins=1,
             min_v=0.8, max_v=0.95, nvbins=20, verbose=2)

widest = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
              min_u=0.9, max_u=1, nubins=1,
              min_v=0.95, max_v=1.0, nvbins=20, verbose=2)

if 1:
    ggg1 = treecorr.GGGCorrelation(narrow, var_method='bootstrap')
    ggg1.process(cats_4, comm=comm)
    if rank == 0:
        ggg1.write('narrow.hdf', write_patch_results=True)

if 1:
    ggg2 = treecorr.GGGCorrelation(wide, var_method='bootstrap')
    ggg2.process(cats_4, comm=comm)
    if rank == 0:
        ggg2.write('wide.hdf', write_patch_results=True)

if 1:
    ggg3 = treecorr.GGGCorrelation(wider, var_method='bootstrap')
    ggg3.process(cats_4, comm=comm)
    if rank == 0:
        ggg3.write('wider.hdf', write_patch_results=True)

if 1:
    ggg4 = treecorr.GGGCorrelation(widest, var_method='bootstrap')
    ggg4.process(cats_4, comm=comm)
    if rank == 0:
        ggg4.write('widest.hdf', write_patch_results=True)
