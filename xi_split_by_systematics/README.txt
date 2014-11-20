xi_split_by_systematics
TODO: add a functionality to use fill p(z) for each galaxy instead of redshift point estimate

homogenise_nz.py module usage:

# import module
import homogenise_nz

# example - get some data
des_catalog = pyfits.getdata(some_catalog)

# make the selection by systematics into bins
select_systematic_bin1 = ... # some selection function to split into bins
select_systematic_bin2 = ... # some selection function to split into bins

# now get the redshift point estimate column - for example ZB, and statistical weights
z_bin1 = des_catalog[select_systematic_bin1_z]['ZB']
w_bin1 = des_catalog[select_systematic_bin1_z]['w']

z_bin2 = des_catalog[select_systematic_bin2_z]['ZB']
w_bin2 = des_catalog[select_systematic_bin2_z]['w']

# how create tuples of redshift and weights
list_split = [ (z_bin1,w_bin1) , (z_bin2,w_bin2) ]

list_nz_weights = homogenise_nz.get_weights(list_split)

nz_weight_bin1 = list_nz_weights[0]
nz_weight_bin2 = list_nz_weights[1]

# proceed with calculating xis 





