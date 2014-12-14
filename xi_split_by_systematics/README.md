# xi_split_by_systematics

homogenise_nz.py module usage has two modes:

## using point estimators for point z estimators

```
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
```

## using point estimators for full p(z)

There are some differences with respect to point-estimator-reweighting:
- you need to have ```DES_pdf_stacker.py``` installed, with branch ```speedup```
get it from https://bitbucket.org/amaraa/wl-photoz, if you can't access it, ask Adam Amara for it.
- you need to have the big h5 files from Christopher Bonnett.
- statistical weights are included at the stage of pulling redshifts from the h5 file, not later.
- there is a parameter that needs to be tuned: ```sigma_regularisation```, with default of ```1e-5```. This parameter is a prior on weights. We would like to make the weights as close to 1 as possible. If if ```sigma_regularisation``` is small, then the prior is weak, and large scatter in weights away from 1 is allowed.
If ```sigma_regularisation``` is big, then the weights are close to one, but at the expense of imperfect redshift homogenisation. Do check the plots produced by ```get_weights_fullPZ``` to assess the quality of reweighting. This potentially could be automated, but I leave it to future volunteers.
- For im3shape, it takes 2GB of RAM to load the redshift tables, for NGMIX it may to up to 6GB. If that's becomming a problem, let me know and I will downsample the redshift vector to save memory.

Example usage:

```
	
	# import module
	import homogenise_nz
	import DES_pdf_stacker
 
	# load data somehow
    res = load_data() ...

    # get a list for coadd ids
    list_snr_bins_cats = []

    # loop over bins, full bin first
    for isb,snr_bin in enumerate(config['snr_bins']):

    	# select galaxies in that bin
        select = (res['SNR'] > snr_bin[0]) & (res['SNR'] < snr_bin[1]) 
        res_bin = res[select]
      
      	# get the redshifts array 
        pdf_array, z_values = DES_pdf_stacker.return_pdf_array(res_bin['coadd_id'],config['name_pzcode'],config['filename_photoz_h5'],weight=res_bin['w'])
        
        # store
        list_snr_bins_cats.append(pdf_array)

    # get the redshift weights
    label = 'snr.%s' % config['method']
    list_weights = homogenise_nz.get_weights_fullPZ(list_snr_bins_cats,z_values=z_values,target_nz_index=0,label=label,plots=True,sigma_regularisation=1e-5)

```



