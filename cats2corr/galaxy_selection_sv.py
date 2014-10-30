import desdb
import scipy.stats
import numpy as np

#Some choices about the redshift range and catalog
#to use
min_redshift = 0.3
max_redshift = 1.3
photoz_catalog = 'joezuntz.tpz_v1'
z_field = 'z_mean'
nbin = 3

#Connect to the DESDM database
connection = desdb.Connection()

#Get all the IDs of objects in the Gold catalog
query = """
SELECT
	gold.coadd_objects_id as id,
	photoz.{z_field} as z
FROM
	des_admin.sva1_gold_basic gold,
	{photoz_catalog} photoz
WHERE
	photoz.coadd_objects_id=gold.coadd_objects_id
	AND
	photoz.{z_field}>{min_redshift}
	AND
	photoz.{z_field}<{max_redshift}
""".format(**locals())

#Run the query and pull out the columns we need
print "Running query"
catalog = connection.quick(query, array=True)
np.save("catalog.npy", catalog)
ids = catalog['id']
z = catalog['z']
n_total = z.size
print "Query complete"

#Get the bin edges by looking at percentiles of the
#redshifts
percentiles = np.linspace(0.,100.,nbin+1)
#for nbin=3: percentiles=[0.0, 33.333, 66.666, 100.0]
bin_edges = scipy.stats.scoreatpercentile(z, percentiles)

#Save these edge values to a text file
edge_header="""
#bin_edges
"""
np.savetxt("bin_edges.txt", bin_edges, header=edge_header)


#Loop through the bins finding the objects in that
#bin and saving their IDs

#This is the header for the ID files (we fill in details in a moment)
bin_header = """
#ids_bin_{i}
#z_min={zmin}
#z_max={zmax}
#n_obj={n_obj}
"""

#Also we want to count up as a little checksum
#that we use all the objects
n_used = 0
for i in xrange(nbin):
	#Get the bin edges for this bin
	zmin = bin_edges[i]
	zmax = bin_edges[i+1]
	print "Bin {0}: z = ({1} = {2})".format(i,zmin,zmax)
	#Find the z values in that bin
	bin_indices = (z>=zmin) & (z<zmax)
	#Find the corresponding coadd_object_ids
	bin_ids = ids[bin_indices]
	#Keep cound of everything
	n_obj = bin_ids.size
	n_used += n_obj
	#Save the bin ids to a text file
	np.savetxt("bin_{0}.txt".format(i), bin_ids, header=bin_header.format(**locals()))

#Just a little check that we use all the objects
if n_used!=n_total:
	print "DID NOT USE ALL GALAXIES - BIN EDGES?"
	print n_used, n_total

