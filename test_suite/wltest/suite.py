from .suite_base import Suite
from . import single_tests
from . import pair_tests

class YellowSuite(Suite):
	#Take the intersection of the catalogs -
	#we will be comparing object-by-object
	intersect=True
	classes = [
		#Tests where we get a stat from each
		#catalog, rather than a comparison stat
		single_tests.MeanE,
		single_tests.EWithSNR,
		single_tests.EWithPSF11,
		single_tests.EWithPSF21,
		single_tests.EWithPSF22,
		single_tests.EWithPSF12,
		
		#Tests where the stat is a comparison
		#of two catalogs
		pair_tests.EDifferenceTest,
	]

	def select(self, cat):
		yellow = cat['error_flag']==0
		return cat[yellow]
