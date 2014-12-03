xi_split_by_systematics
TODO: add a functionality to use fill p(z) for each galaxy instead of redshift point estimate

homogenise_nz.py module usage:

# import module
import homogenise_nz

# example - get some data
des_catalog = pyfits.getdata(some_catalog)

# make the selection by systematics into bins
# let's say bin 0 is a bin taking all galaxies without splitting
select_systematic_bin0 = ... # some selection taking all galaxies 
select_systematic_bin1 = ... # some selection function to split into bins
select_systematic_bin2 = ... # some selection function to split into bins

# now get the redshift point estimate column - for example ZB, and statistical weights
z_bin0 = des_catalog[select_systematic_bin0_z]['ZB']
w_bin0 = des_catalog[select_systematic_bin0_z]['w']
z_bin1 = des_catalog[select_systematic_bin1_z]['ZB']
w_bin1 = des_catalog[select_systematic_bin1_z]['w']
z_bin2 = des_catalog[select_systematic_bin2_z]['ZB']
w_bin2 = des_catalog[select_systematic_bin2_z]['w']

# how create tuples of redshift and weights
list_split = [ (z_bin1,w_bin0), (z_bin1,w_bin1) , (z_bin2,w_bin2) ]

# run the code to get nz_weights
# if you want to use a target nz distribution to be that of one of the bins (let's say a bin without splits), use target_nz_index = 0 (in our case 0 was the index of the bin without splits)
# if you want to use mean of bins as a target n(z), then use target_nz_index=-1
list_nz_weights = homogenise_nz.get_weights(list_split,target_nz_index=0,photoz_min=0.2,photoz_max=1.2,photoz_nbins=50)

nz_weight_bin0 = list_nz_weights[0] # these should be all one (or very close) as we used this bin as a target n(z)
nz_weight_bin1 = list_nz_weights[1]
nz_weight_bin2 = list_nz_weights[2]

# proceed with calculating xis 





