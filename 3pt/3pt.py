import treecorr
import glob, os

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

wide = dict(min_sep=1, max_sep=5, sep_units='arcmin', nbins=5,
            min_u=0.9, max_u=1, nubins=1,
            min_v=0.0, max_v=1.0, nvbins=50, verbose=2)

narrow = dict(min_sep=1, max_sep=5, sep_units='arcmin', nbins=5,
              min_u=0.0, max_u=1, nubins=20,
              min_v=0.0, max_v=0.1, nvbins=1, verbose=2)

ggg1 = treecorr.GGGCorrelation(wide, var_method='bootstrap')
ggg1.process(cats_4)
ggg1.write('wide.hdf', write_patch_results=True)

ggg2 = treecorr.GGGCorrelation(narrow, var_method='bootstrap')
ggg2.process(cats_4)
ggg2.write('narrow.hdf', write_patch_results=True)
