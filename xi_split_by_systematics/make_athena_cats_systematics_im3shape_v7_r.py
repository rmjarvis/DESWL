import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pylab
import os
import sys
import pyfits
import glob
sys.path.append('/Users/drgk/DES/SV_tests/DESWL/xi_split_by_systematics')
import homogenise_nz


# Use the im3shape_v6 fits files
# Read in
# Cut with modest proposal
# cut with any other flags
# make maps

nside=2048 # set up empty map
#R=hp.Rotator(coord='cg')
R=hp.Rotator(coord='cg')
gamma1map=np.zeros((12*nside**2))
gamma2map=np.zeros((12*nside**2))
countmap1=np.zeros((12*nside**2))


ra_vals =[]
dec_vals =[]
e1_vals =[]
e2_vals =[]
c1_vals =[]
c2_vals =[]
m_vals =[]
w_vals =[]
error_flag =[]
info_flag =[]
coadd_objects_id =[]

# Read in all fits files
x = glob.glob("/Users/drgk/DES/data/im3shape_cats_v7_sept2014/*r.fits.gz")

k=1;
print 'Reading in fits files'
for ifits in x:
	print str(100*(k+1)/float(len(x)))+' percent complete'

	fits_file = pyfits.open(ifits)

	ra_tmp = fits_file[1].data.field('ra')
	dec_tmp = fits_file[1].data.field('dec')
	e1_tmp = fits_file[1].data.field('e1')
	e2_tmp = fits_file[1].data.field('e2')
	c1_tmp = fits_file[1].data.field('nbc_c1')
	c2_tmp = fits_file[1].data.field('nbc_c2')
	m_tmp  = fits_file[1].data.field('nbc_m')
	w_tmp  = fits_file[1].data.field('w')
	error_flag_tmp = fits_file[1].data.field('error_flag')
	info_flag_tmp = fits_file[1].data.field('info_flag')
	coadd_objects_id_tmp = fits_file[1].data.field('coadd_objects_id')

	ra_vals.extend(ra_tmp)
	dec_vals.extend(dec_tmp)
	e1_vals.extend(e1_tmp)
	e2_vals.extend(e2_tmp)
	c1_vals.extend(c1_tmp)
	c2_vals.extend(c2_tmp)
	m_vals.extend(m_tmp)
	w_vals.extend(w_tmp)
	error_flag.extend(error_flag_tmp)
	info_flag.extend(info_flag_tmp)
	coadd_objects_id.extend(coadd_objects_id_tmp)


	k+=1

ra_vals = np.asanyarray(ra_vals)
dec_vals = np.asanyarray(dec_vals)
e1_vals = np.asanyarray(e1_vals)
e2_vals = np.asanyarray(e2_vals)
c1_vals = np.asanyarray(c1_vals)
c2_vals = np.asanyarray(c2_vals)
m_vals = np.asanyarray(m_vals)
w_vals = np.asanyarray(w_vals)
error_flag = np.asanyarray(error_flag)
info_flag = np.asanyarray(info_flag)
coadd_objects_id = np.asanyarray(coadd_objects_id)
N_fromcat = len(ra_vals)

## CUT ON MODEST PROPOSAL
# (FLAGS_I <=3) AND NOT ( ( (CLASS_STAR_I > 0.3) AND (MAG_AUTO_I < 18.0) ) OR ((SPREAD_MODEL_I + 3*SPREADERR_MODEL_I) < 0.003) OR ((MAG_PSF_I > 30.0) AND (MAG_AUTO_I < 21.0))))
print 'Star/Galaxy Separation: \n'
x = pyfits.open('/Users/drgk/DES/data/im3shape_cats_v7_sept2014/im3shape_v7_r_modest_matched.fits')
coadd_objects_id_modestmatched = x[1].data.field('coadd_objects_id')

x = np.in1d(coadd_objects_id,coadd_objects_id_modestmatched)
coadd_objects_id_modest = coadd_objects_id[x]
ra_modest = ra_vals[x]
dec_modest = dec_vals[x]
e1_modest = e1_vals[x]
e2_modest = e2_vals[x]
c1_modest = c1_vals[x]
c2_modest = c2_vals[x]
m_modest = m_vals[x]
w_modest = w_vals[x]
error_flag_modest = error_flag[x]
info_flag_modest = info_flag[x]
N_modest = len(ra_modest)

## CUT ON OTHER FLAGS
print 'Cut on Flags: \n'

ra_modest_errorflag = ra_modest[np.where(error_flag_modest==0)]
dec_modest_errorflag = dec_modest[np.where(error_flag_modest==0)]
e1_modest_errorflag = e1_modest[np.where(error_flag_modest==0)]
e2_modest_errorflag = e2_modest[np.where(error_flag_modest==0)]
c1_modest_errorflag = c1_modest[np.where(error_flag_modest==0)]
c2_modest_errorflag = c2_modest[np.where(error_flag_modest==0)]
m_modest_errorflag  = m_modest[np.where(error_flag_modest==0)]
w_modest_errorflag  = w_modest[np.where(error_flag_modest==0)]
coadd_objects_id_modest_errorflag = coadd_objects_id_modest[np.where(error_flag_modest==0)]
info_flag_modest_errorflag = info_flag_modest[np.where(error_flag_modest==0)]
N_modest_errorflag = len(ra_modest_errorflag)

ra_modest_errorflag_infoflag = ra_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
dec_modest_errorflag_infoflag = dec_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
e1_modest_errorflag_infoflag = e1_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
e2_modest_errorflag_infoflag = e2_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
c1_modest_errorflag_infoflag = c1_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
c2_modest_errorflag_infoflag = c2_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
m_modest_errorflag_infoflag = m_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
w_modest_errorflag_infoflag = w_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
coadd_objects_id_modest_errorflag_infoflag = coadd_objects_id_modest_errorflag[np.where(info_flag_modest_errorflag==0)]
N_modest_errorflag_infoflag = len(ra_modest_errorflag_infoflag)

# RESTRICT TO SPT-E
include_spte = np.where((ra_modest_errorflag_infoflag>56)&(ra_modest_errorflag_infoflag<94)&(dec_modest_errorflag_infoflag<-42))
include_sptw = np.where((ra_modest_errorflag_infoflag>340)&(dec_modest_errorflag_infoflag<-51))

ra_modest_errorflag_infoflag_spte = ra_modest_errorflag_infoflag[include_spte]
dec_modest_errorflag_infoflag_spte = dec_modest_errorflag_infoflag[include_spte]
e1_modest_errorflag_infoflag_spte = e1_modest_errorflag_infoflag[include_spte]
e2_modest_errorflag_infoflag_spte = e2_modest_errorflag_infoflag[include_spte]
c1_modest_errorflag_infoflag_spte = c1_modest_errorflag_infoflag[include_spte]
c2_modest_errorflag_infoflag_spte = c2_modest_errorflag_infoflag[include_spte]
m_modest_errorflag_infoflag_spte  = m_modest_errorflag_infoflag[include_spte]
w_modest_errorflag_infoflag_spte  = w_modest_errorflag_infoflag[include_spte]
coadd_objects_id_modest_errorflag_infoflag_spte = coadd_objects_id_modest_errorflag_infoflag[include_spte]
N_modest_errorflag_infoflag_spte = len(ra_modest_errorflag_infoflag_spte)

## Match to redshifts
x = pyfits.open('/Users/drgk/DES/data/photoz/gold_neuralnetwork/sva1_gold_1.0_catalog_desdmphotoz.fits.gz')
coadd_objects_id_photoz = x[1].data.field('coadd_objects_id')
z_photoz = x[1].data.field('zp')
z_modest_errorflag_infoflag_spte = np.zeros(np.shape(ra_modest_errorflag_infoflag_spte))
N = len(ra_modest_errorflag_infoflag_spte)
print 'Match to redshifts: \n'

index_dict = {}
index_dict = {key: indexval for (indexval, key) in enumerate(coadd_objects_id_photoz)}
indices_array = np.array([index_dict[single_coadd_object_id] for single_coadd_object_id in coadd_objects_id_modest_errorflag_infoflag_spte])
z_modest_errorflag_infoflag_spte = z_photoz[indices_array]

# cut to 0.3 < z < 1.3
include_zcut = np.where((z_modest_errorflag_infoflag_spte>=0.3)&(z_modest_errorflag_infoflag_spte<1.3))
ra_modest_errorflag_infoflag_spte_zcut = ra_modest_errorflag_infoflag_spte[include_zcut]
dec_modest_errorflag_infoflag_spte_zcut = dec_modest_errorflag_infoflag_spte[include_zcut]
e1_modest_errorflag_infoflag_spte_zcut = e1_modest_errorflag_infoflag_spte[include_zcut]
e2_modest_errorflag_infoflag_spte_zcut = e2_modest_errorflag_infoflag_spte[include_zcut]
c1_modest_errorflag_infoflag_spte_zcut = c1_modest_errorflag_infoflag_spte[include_zcut]
c2_modest_errorflag_infoflag_spte_zcut = c2_modest_errorflag_infoflag_spte[include_zcut]
m_modest_errorflag_infoflag_spte_zcut = m_modest_errorflag_infoflag_spte[include_zcut]
w_modest_errorflag_infoflag_spte_zcut = w_modest_errorflag_infoflag_spte[include_zcut]
z_modest_errorflag_infoflag_spte_zcut = z_modest_errorflag_infoflag_spte[include_zcut]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_full_zmin03_zmax13.dat',X,delimiter =' ')

# convert ra/dec into pixels
degr2rad=np.pi/180.;

theta_modest_errorflag_infoflag_spte_zcut =-degr2rad*dec_modest_errorflag_infoflag_spte_zcut +np.pi/2.
phi_modest_errorflag_infoflag_spte_zcut =degr2rad*ra_modest_errorflag_infoflag_spte_zcut 
pixs=hp.ang2pix(4096, theta_modest_errorflag_infoflag_spte_zcut, phi_modest_errorflag_infoflag_spte_zcut) # turn theta/phi from degrees into radians, nest = False

# ## Full Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_03_z_13.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_03_z_13.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut)))+' percent complete: shear Full'

	ra = ra_modest_errorflag_infoflag_spte_zcut[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut[i] - c1_modest_errorflag_infoflag_spte_zcut[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut[i] - c2_modest_errorflag_infoflag_spte_zcut[i]
	m = m_modest_errorflag_infoflag_spte_zcut[i]
	w = w_modest_errorflag_infoflag_spte_zcut[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

######################## AIRMASS ###########################
# read in airmass r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS__mean.fits.gz')

airmass_pixels = x[1].data.field('pixel')
airmass_r        = x[1].data.field('signal')
airmass_r_median = np.median(airmass_r)

map_airmass_mean = np.zeros(12*4096**2)
map_airmass_mean[airmass_pixels] = airmass_r

full_pix_lower = np.where((map_airmass_mean<airmass_r_median)&(map_airmass_mean>0))
full_pix_upper = np.where(map_airmass_mean>=airmass_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_airmass_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_airmass_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_airmass_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_airmass_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_airmass_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_airmass_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_airmass_lower  = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_airmass_lower  = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_airmass_lower  = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_airmass_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_airmass_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_airmass_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_airmass_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_airmass_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_airmass_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_airmass_upper  = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_airmass_upper  = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_airmass_upper  = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_airmass_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_airmass_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_airmass_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_airmass_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_airmass_upper,w_modest_errorflag_infoflag_spte_zcut_airmass_upper) , (z_modest_errorflag_infoflag_spte_zcut_airmass_lower,w_modest_errorflag_infoflag_spte_zcut_airmass_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper Airmass Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_airmass_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_airmass_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_airmass_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_airmass_upper)))+' percent complete: shear Airmass Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_airmass_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_airmass_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_airmass_upper[i]*nz_weight_bin1[i]
	
	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower Airmass Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_airmass_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_airmass_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_airmass_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_airmass_lower)))+' percent complete: shear Airmass Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_airmass_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_airmass_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_airmass_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

######################## EXPOSURE TIME ###########################
# read in EXPTIME r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz')

exptime_pixels = x[1].data.field('pixel')
exptime_r        = x[1].data.field('signal')
exptime_r_median = np.median(exptime_r)

map_exptime_mean = np.zeros(12*4096**2)
map_exptime_mean[exptime_pixels] = exptime_r

full_pix_lower = np.where((map_exptime_mean<exptime_r_median)&(map_exptime_mean>0))
full_pix_upper = np.where(map_exptime_mean>=exptime_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_exptime_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_exptime_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_exptime_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_exptime_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_exptime_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_exptime_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_exptime_lower = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_exptime_lower = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_exptime_lower = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_exptime_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_exptime_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_exptime_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_exptime_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_exptime_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_exptime_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_exptime_upper = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_exptime_upper = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_exptime_upper = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_exptime_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_exptime_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_exptime_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_exptime_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_exptime_upper,w_modest_errorflag_infoflag_spte_zcut_exptime_upper) , (z_modest_errorflag_infoflag_spte_zcut_exptime_lower,w_modest_errorflag_infoflag_spte_zcut_exptime_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper exptime Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_exptime_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_exptime_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_exptime_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_exptime_upper)))+' percent complete: shear exptime Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_exptime_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_exptime_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_exptime_upper[i]*nz_weight_bin1[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower exptime Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_exptime_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_exptime_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_exptime_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_exptime_lower)))+' percent complete: shear exptime Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_exptime_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_exptime_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_exptime_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

######################## FWHM ###########################
# read in FWHM r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM__mean.fits.gz')

fwhm_pixels = x[1].data.field('pixel')
fwhm_r        = x[1].data.field('signal')
fwhm_r_median = np.median(fwhm_r)

map_fwhm_mean = np.zeros(12*4096**2)
map_fwhm_mean[fwhm_pixels] = fwhm_r

full_pix_lower = np.where((map_fwhm_mean<fwhm_r_median)&(map_fwhm_mean>0))
full_pix_upper = np.where(map_fwhm_mean>=fwhm_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_fwhm_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_fwhm_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_fwhm_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_fwhm_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_fwhm_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_fwhm_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_fwhm_lower = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_fwhm_lower = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_fwhm_lower = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_fwhm_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_fwhm_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_fwhm_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_fwhm_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_fwhm_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_fwhm_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_fwhm_upper = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_fwhm_upper = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_fwhm_upper = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_fwhm_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_fwhm_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_fwhm_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_fwhm_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_fwhm_upper,w_modest_errorflag_infoflag_spte_zcut_fwhm_upper) , (z_modest_errorflag_infoflag_spte_zcut_fwhm_lower,w_modest_errorflag_infoflag_spte_zcut_fwhm_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper fwhm Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_fwhm_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_fwhm_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_fwhm_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_fwhm_upper)))+' percent complete: shear fwhm Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_fwhm_upper[i]*nz_weight_bin1[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower fwhm Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_fwhm_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_fwhm_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_fwhm_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_fwhm_lower)))+' percent complete: shear fwhm Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_fwhm_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()


######################## Magnitude Limit ###########################
# read in maglimit r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit__.fits.gz')

maglimit_pixels = x[1].data.field('pixel')
maglimit_r        = x[1].data.field('signal')
maglimit_r_median = np.median(maglimit_r)

map_maglimit_mean = np.zeros(12*4096**2)
map_maglimit_mean[maglimit_pixels] = maglimit_r

full_pix_lower = np.where((map_maglimit_mean<maglimit_r_median)&(map_maglimit_mean>0))
full_pix_upper = np.where(map_maglimit_mean>=maglimit_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_maglimit_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_maglimit_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_maglimit_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_maglimit_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_maglimit_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_maglimit_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_maglimit_lower = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_maglimit_lower = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_maglimit_lower = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_maglimit_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_maglimit_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_maglimit_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_maglimit_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_maglimit_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_maglimit_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_maglimit_upper = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_maglimit_upper = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_maglimit_upper = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_maglimit_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_maglimit_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_maglimit_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_maglimit_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_maglimit_upper,w_modest_errorflag_infoflag_spte_zcut_maglimit_upper) , (z_modest_errorflag_infoflag_spte_zcut_maglimit_lower,w_modest_errorflag_infoflag_spte_zcut_maglimit_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper maglimit Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_maglimit_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_maglimit_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_maglimit_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_maglimit_upper)))+' percent complete: shear maglimit Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_maglimit_upper[i]*nz_weight_bin1[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower maglimit Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_maglimit_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_maglimit_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_maglimit_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_maglimit_lower)))+' percent complete: shear maglimit Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_maglimit_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()


######################## Sky Brightness ###########################
# read in skybrite r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE__mean.fits.gz')

skybrite_pixels = x[1].data.field('pixel')
skybrite_r        = x[1].data.field('signal')
skybrite_r_median = np.median(skybrite_r)

map_skybrite_mean = np.zeros(12*4096**2)
map_skybrite_mean[skybrite_pixels] = skybrite_r

full_pix_lower = np.where((map_skybrite_mean<skybrite_r_median)&(map_skybrite_mean>0))
full_pix_upper = np.where(map_skybrite_mean>=skybrite_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_skybrite_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_skybrite_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_skybrite_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_skybrite_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_skybrite_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_skybrite_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_skybrite_lower = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_skybrite_lower = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_skybrite_lower = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_skybrite_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_skybrite_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_skybrite_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_skybrite_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_skybrite_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_skybrite_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_skybrite_upper = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_skybrite_upper = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_skybrite_upper = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_skybrite_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_skybrite_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_skybrite_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_skybrite_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_skybrite_upper,w_modest_errorflag_infoflag_spte_zcut_skybrite_upper) , (z_modest_errorflag_infoflag_spte_zcut_skybrite_lower,w_modest_errorflag_infoflag_spte_zcut_skybrite_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper skybrite Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skybrite_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skybrite_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_skybrite_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_skybrite_upper)))+' percent complete: shear skybrite Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_skybrite_upper[i]*nz_weight_bin1[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower skybrite Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skybrite_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skybrite_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_skybrite_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_skybrite_lower)))+' percent complete: shear skybrite Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_skybrite_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

######################## Sky Sigma ###########################
# read in skysigma r

x = pyfits.open('/Users/drgk/DES/data/systematics/nside4096_oversamp4/SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA__mean.fits.gz')

skysigma_pixels = x[1].data.field('pixel')
skysigma_r        = x[1].data.field('signal')
skysigma_r_median = np.median(skysigma_r)

map_skysigma_mean = np.zeros(12*4096**2)
map_skysigma_mean[skysigma_pixels] = skysigma_r

full_pix_lower = np.where((map_skysigma_mean<skysigma_r_median)&(map_skysigma_mean>0))
full_pix_upper = np.where(map_skysigma_mean>=skysigma_r_median)

include_upper = np.in1d(pixs,full_pix_upper)
include_lower = np.in1d(pixs,full_pix_lower)
pixs_lower = pixs[include_lower]
ra_modest_errorflag_infoflag_spte_zcut_skysigma_lower = ra_modest_errorflag_infoflag_spte_zcut[include_lower]
dec_modest_errorflag_infoflag_spte_zcut_skysigma_lower = dec_modest_errorflag_infoflag_spte_zcut[include_lower]
e1_modest_errorflag_infoflag_spte_zcut_skysigma_lower = e1_modest_errorflag_infoflag_spte_zcut[include_lower]
e2_modest_errorflag_infoflag_spte_zcut_skysigma_lower = e2_modest_errorflag_infoflag_spte_zcut[include_lower]
c1_modest_errorflag_infoflag_spte_zcut_skysigma_lower = c1_modest_errorflag_infoflag_spte_zcut[include_lower]
c2_modest_errorflag_infoflag_spte_zcut_skysigma_lower = c2_modest_errorflag_infoflag_spte_zcut[include_lower]
m_modest_errorflag_infoflag_spte_zcut_skysigma_lower = m_modest_errorflag_infoflag_spte_zcut[include_lower]
w_modest_errorflag_infoflag_spte_zcut_skysigma_lower = w_modest_errorflag_infoflag_spte_zcut[include_lower]
z_modest_errorflag_infoflag_spte_zcut_skysigma_lower = z_modest_errorflag_infoflag_spte_zcut[include_lower]
pixs_upper = pixs[include_upper]
ra_modest_errorflag_infoflag_spte_zcut_skysigma_upper = ra_modest_errorflag_infoflag_spte_zcut[include_upper]
dec_modest_errorflag_infoflag_spte_zcut_skysigma_upper = dec_modest_errorflag_infoflag_spte_zcut[include_upper]
e1_modest_errorflag_infoflag_spte_zcut_skysigma_upper = e1_modest_errorflag_infoflag_spte_zcut[include_upper]
e2_modest_errorflag_infoflag_spte_zcut_skysigma_upper = e2_modest_errorflag_infoflag_spte_zcut[include_upper]
c1_modest_errorflag_infoflag_spte_zcut_skysigma_upper = c1_modest_errorflag_infoflag_spte_zcut[include_upper]
c2_modest_errorflag_infoflag_spte_zcut_skysigma_upper = c2_modest_errorflag_infoflag_spte_zcut[include_upper]
m_modest_errorflag_infoflag_spte_zcut_skysigma_upper = m_modest_errorflag_infoflag_spte_zcut[include_upper]
w_modest_errorflag_infoflag_spte_zcut_skysigma_upper = w_modest_errorflag_infoflag_spte_zcut[include_upper]
z_modest_errorflag_infoflag_spte_zcut_skysigma_upper = z_modest_errorflag_infoflag_spte_zcut[include_upper]

a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_skysigma_lower,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_skysigma_lower_zmin03_zmax13.dat',X,delimiter =' ')
a = plt.hist(z_modest_errorflag_infoflag_spte_zcut_skysigma_upper,200); X = np.zeros([200,2]); binsx = a[1]; bins = (binsx[0:-1]+binsx[1:])/2.0; X[:,0] = bins; X[:,1] = a[0]/np.sum(a[0]); np.savetxt('/Users/drgk/DES/data/systematics/nz_v7_r_skysigma_upper_zmin03_zmax13.dat',X,delimiter =' ')

list_split = [ (z_modest_errorflag_infoflag_spte_zcut,w_modest_errorflag_infoflag_spte_zcut), (z_modest_errorflag_infoflag_spte_zcut_skysigma_upper,w_modest_errorflag_infoflag_spte_zcut_skysigma_upper) , (z_modest_errorflag_infoflag_spte_zcut_skysigma_lower,w_modest_errorflag_infoflag_spte_zcut_skysigma_lower) ]

list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3,photoz_nbins=50)
nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# ## Upper skysigma Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skysigma_r_upper.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skysigma_r_upper.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_skysigma_upper)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_skysigma_upper)))+' percent complete: shear skysigma Upper'

	ra = ra_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i] - c1_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i] - c2_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]
	m = m_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]
	w = w_modest_errorflag_infoflag_spte_zcut_skysigma_upper[i]*nz_weight_bin1[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()

# ## Lower skysigma Catalogue
f = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_shears_skysigma_r_lower.dat','w')
g = open('/Users/drgk/DES/SV_tests/athena_cats/im3shape_v7_r_m_skysigma_r_lower.dat','w')

for i in range(0,len(ra_modest_errorflag_infoflag_spte_zcut_skysigma_lower)):
	print str(100*(i+1)/float(len(ra_modest_errorflag_infoflag_spte_zcut_skysigma_lower)))+' percent complete: shear skysigma Lower'

	ra = ra_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]
	dec = dec_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]
	e1 = e1_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i] - c1_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]
	e2 = e2_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i] - c2_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]
	m = m_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]
	w = w_modest_errorflag_infoflag_spte_zcut_skysigma_lower[i]*nz_weight_bin2[i]

	f.write(str(ra)+'   '+str(dec)+'   '+str(e1)+'   '+str(e2)+'   '+str(w)+'\n')
	g.write(str(ra)+'   '+str(dec)+'   '+str(1.0+m)+'   '+str(w)+'\n')

f.close()
g.close()


sys.exit()


### Plot histograms

plt.figure()
plt.subplot(3,2,1)
plt.hist(airmass_r,200)
plt.vlines(airmass_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.ylim([0,80000])
plt.xlabel('Airmass')
plt.subplot(3,2,2)
plt.hist(exptime_r,200)
plt.vlines(exptime_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.ylim([0,150000])
plt.xlabel('Exposure Time')
plt.subplot(3,2,3)
plt.hist(fwhm_r,200)
plt.vlines(fwhm_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.xlabel('FWHM')
plt.ylim([0,50000])
plt.subplot(3,2,4)
plt.hist(maglimit_r,200)
plt.vlines(maglimit_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.xlabel('Magnitude Limit')
plt.ylim([0,100000])
plt.subplot(3,2,5)
plt.hist(skybrite_r,200)
plt.vlines(skybrite_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.xlabel('Sky Brightness')
plt.ylim([0,700000])
plt.subplot(3,2,6)
plt.hist(skysigma_r,200)
plt.vlines(skysigma_r_median,0,800000,color='red',label='median')
plt.legend(loc=0)
plt.xlabel('Sky Sigma')
plt.ylim([0,300000])
plt.show()



