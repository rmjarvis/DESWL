import treecorr
import glob, os

data_dir = '/global/homes/s/seccolf/shear3pt/SHARED_DATA'

zbin1_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin1', '*.fits')))
zbin2_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin2', '*.fits')))
zbin3_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin3', '*.fits')))
zbin4_files = sorted(glob.glob(os.path.join(data_dir, 'T17', 'patches_zbin4', '*.fits')))
mice_files = sorted(glob.glob(os.path.join(data_dir, 'MICE', '*.fits')))

params_cat = dict(ra_col='ra', dec_col='dec', ra_units='rad', dec_units='rad',
                  w_col='w', g1_col='g1', g2_col='g2', flip_g1=True)

cats_1 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin1_files)]
cats_2 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin2_files)]
cats_3 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin3_files)]
cats_4 = [treecorr.Catalog(f, params_cat, patch=p) for p,f in enumerate(zbin4_files)]
cats_m = [treecorr.Catalog(f, params_cat, patch=p, w_col='0') for p,f in enumerate(mice_files)]

if 0:
    params_2pt = dict(min_sep=1, max_sep=200, sep_units='arcmin', nbins=10)

    gg = treecorr.GGCorrelation(params_2pt, var_method='bootstrap')

    print('zbin1: GG')
    gg.process(cats_2)
    print(gg.xip)
    print(gg.xim)
    print(gg.varxip**0.5)
    print(gg.varxim**0.5)

    print('zbin4: GG')
    gg.process(cats_4)
    print(gg.xip)
    print(gg.xim)
    print(gg.varxip**0.5)
    print(gg.varxim**0.5)

    print('zbin1-4: NG')
    ng = treecorr.NGCorrelation(params_2pt, var_method='bootstrap')
    ng.process(cats_1, cats_4)
    print(ng.xi)
    print(ng.varxi**0.5)

equilateral = dict(min_sep=1, max_sep=5, sep_units='arcmin', nbins=1,
                   min_u=0.9, max_u=1, nubins=1,
                   min_v=0, max_v=0.1, nvbins=1, verbose=2)
isoceles = dict(min_sep=1, max_sep=5, sep_units='arcmin', nbins=5,
                min_u=0.9, max_u=1, nubins=1,
                min_v=0.7, max_v=1.0, nvbins=10, verbose=2)
params_3pt = isoceles
#params_3pt = equilateral
ggg = treecorr.GGGCorrelation(params_3pt, var_method='bootstrap')
ggg.process(cats_4)
print(f'ntri = {ggg.ntri}')
print(f'w = {ggg.weight}')
print(f'meand1 = {ggg.meand1}')
print(f'meand2 = {ggg.meand2}')
print(f'meand3 = {ggg.meand3}')
print(f'meanu = {ggg.meanu}')
print(f'meanv = {ggg.meanv}')
print(f'gam0 = {ggg.gam0}')
print(f'gam1 = {ggg.gam1}')
print(f'gam2 = {ggg.gam2}')
print(f'gam3 = {ggg.gam3}')
print(f'sigma_gam0 = {ggg.vargam0**0.5}')
print(f'sigma_gam1 = {ggg.vargam1**0.5}')
print(f'sigma_gam2 = {ggg.vargam2**0.5}')
print(f'sigma_gam3 = {ggg.vargam3**0.5}')

ggg.write('ggg.hdf', write_patch_results=True)
