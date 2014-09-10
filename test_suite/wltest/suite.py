from .suite_base import Suite, TICK, CROSS
from . import single_tests
from . import pair_tests
import numpy as np

class GreenSuite(Suite):
	#Take the intersection of the catalogs -
	#we will be comparing object-by-object
	intersect=True
	classes = [
		#Tests where we get a stat from each
		#catalog, rather than a comparison stat
		single_tests.MeanE,
		single_tests.HighSNRMeanE,
		single_tests.EWithSNR,
		single_tests.EWithPSF11,
		single_tests.EWithPSF21,
		single_tests.EWithPSF22,
		single_tests.EWithPSF12,
		
		#Tests where the stat is a comparison
		#of two catalogs
		pair_tests.EDifferenceTest,
		pair_tests.PSFDifferenceTest11,
	]

	def select(self, cat):
		green = (cat['info_flag']==0)#&(cat['info_flag'])
		return cat[green]


class Extrapolation(object):
	extrapolations = [
		(pair_tests.EDifferenceTest, single_tests.MeanE),
		(pair_tests.PSFDifferenceTest11, single_tests.EWithPSF11)
	]

	def find_suite_test(self, suite, test):
		for i,t in enumerate(suite.classes):
			if t is test:
				return i
		raise ValueError("Test %s not found in %s"%(test, suite))

	def run(self, testbed_suite, full_suite):
		print "Extrapolating from test-bed to full:"

		for (pair_type, single_type) in self.extrapolations:
			i1 = self.find_suite_test(testbed_suite, pair_type)
			i2 = self.find_suite_test(full_suite, single_type)
			pair_test = testbed_suite.tests[i1]
			single_test = full_suite.tests[i2]

			delta = pair_test.statistic
			base = single_test.statistic
			target = single_test.statistic_target

			extrapolated = base-delta

			if np.all(abs(extrapolated)<=target):
				result = TICK
			else: 
				result = CROSS

			print "    %s = %s = (%s - %s) (target %s)  %s" % (
				single_test.name, extrapolated, base, delta, target, result)
			print 



