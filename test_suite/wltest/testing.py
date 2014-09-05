from . catalogs import table_from_fits, intersection
from glob import glob

def xtest_read():
	filenames = glob("/Users/jaz/working/wl-tests/v7-testbed/*.fits")
	table = table_from_fits(filenames)
	assert len(table)>0

def xtest_intersect():
	filenames = glob("/Users/jaz/working/wl-tests/v7-testbed/*.fits")
	table = table_from_fits(filenames)
	n = len(table)
	t1, t2 = intersection(table,table)
	assert len(t1) == n
	assert len(t2) == n

	t1, t2 = intersection(table[0:200], table[100:300])
	assert len(t1)==100
	assert len(t2)==100
	assert t1['coadd_objects_id'][0] == table['coadd_objects_id'][100]



