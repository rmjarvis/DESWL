"""

Methods for loading in astropy tables from:
 - DESDM
 - FITS files

"""
import astropy.table
import numpy as np
import glob

class Catalog(astropy.table.Table):
	@classmethod
	def from_multiple_fits(cls, filenames, name, quiet=True):
		tables = []
		for filename in filenames:
			print " - ", filename
			tables.append(astropy.table.Table.read(filename, format='fits'))
		print
		if len(tables)>1:
			cat = astropy.table.vstack(tables, join_type='exact', metadata_conflicts='silent')
		else:
			cat = tables[0]
		cat = cls(cat)
		cat.name = name
		return cat

	@classmethod
	def from_directory(cls, dirname):
		print "Loading from directory: ", dirname
		filenames = glob.glob(dirname+"/*.fits") + glob.glob(dirname+"/*.fits.gz")
		cat = cls.from_multiple_fits(filenames, dirname)
		return cat


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

