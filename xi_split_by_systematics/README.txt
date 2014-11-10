xi_split_by_systematics

homogenise_nz.py module usage:

# import module
import homogenise_nz

# example - get some data
des_catalog = pyfits.getdata(some_catalog)

# make the selection by systematics into bins
select_systematic_bin1 = ... # some selection function to split into bins
select_systematic_bin2 = ... # some selection function to split into bins

# now get the redshift point estimate column - for example ZB
z_bin1 = des_catalog[select_systematic_bin1_z]['ZB']
z_bin2 = des_catalog[select_systematic_bin1_z]['ZB']

list_weights = homogenise_nz.get_weights([z_bin1,z_bin2])

weight_bin1 = list_weights[0]
weight_bin2 = list_weights[1]

# proceed with calculating xis 


