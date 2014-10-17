"""

Methods for loading in astropy tables from:
 - DESDM
 - FITS files

"""
import astropy.table
import numpy as np
import glob
import os

USELESS_COLUMNS = [
	'covmat',
	'levmar',
	'fails',
	'time',
]

class Catalog(astropy.table.Table):
	@classmethod
	def from_multiple_fits(cls, filenames, cat_name, quiet=True):
		tables = []
		for filename in filenames:
			print " - ", filename
			cat = astropy.table.Table.read(filename, format='fits')
			removals = []
			for name in cat.colnames:
				for useless in USELESS_COLUMNS:
					if name.startswith(useless):
						removals.append(name)
						break
			cat.remove_columns(removals)
			tables.append(cat)
		print
		if len(tables)>1:
			cat = astropy.table.vstack(tables, join_type='exact', metadata_conflicts='silent')
		else:
			cat = tables[0]
		cat = cls(cat)
		cat.name = cat_name
		return cat

	@classmethod
	def from_directory(cls, dirname):
		print "Loading from directory: ", dirname
		filenames = glob.glob(dirname+"/*.fits") + glob.glob(dirname+"/*.fits.gz")
		cat_name=dirname.strip(os.path.sep).split(os.path.sep)[-1]
		cat = cls.from_multiple_fits(filenames, cat_name)
		return cat

	@classmethod
	def from_desdb(cls,query,save_as=None):
		"""Get catalog by querying the desdm database. Need desdb installed.
		Arguments: - query: which is an sql query as a string
				   - save_as: filename to save catalog returned by query as .npy file"""
		try:
			import desdb
		except Exception as e:
			print 'while trying to import desdb, got the following exception:'
			print e
			return 1
		conn=desdb.Connection()
		catalog=conn.quick(query,array=True)
		if save_as is not None:
			np.save(save_as,catalog)
		print catalog
		return astropy.table.Table(catalog)


	def intersection(self, cat2, field='coadd_objects_id'):
		ids1 = set(self[field])
		ids2 = set(cat2[field])
		ids_both = ids1.intersection(ids2)

		indices = in_table(self[field], ids_both)
		c1 = self[indices]

		indices = in_table(cat2[field], ids_both)
		c2 = cat2[indices]

		return c1, c2

	def __getitem__(self, *args, **kwargs):
		cat = super(Catalog, self).__getitem__(*args, **kwargs)
		if isinstance(cat, Catalog):
			cat.name = self.name
		return cat



def in_table(col, s):
	#This is not fast 
	indices = []
	for i in xrange(len(col)):
		if col[i] in s:
			indices.append(i)
	indices = np.array(indices)
	return indices

