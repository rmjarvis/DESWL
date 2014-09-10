from .catalogs import Catalog
from .suite import GreenSuite, Extrapolation
import sys

def suite_compare():
	dirname1 = 	sys.argv[1]
	dirname2 = 	sys.argv[2]

	cat1 = Catalog.from_directory(dirname1)
	cat2 = Catalog.from_directory(dirname2)
	suite = GreenSuite()
	print
	print "Running tests"
	print
	filenames = suite.run(cat1, cat2)
	print
	if filenames:
		print "Made files:"
		for filename in filenames:
			print " - ", filename

def suite_extrapolate():
	testbed_old = sys.argv[1]
	testbed_new = sys.argv[2]
	full_old = sys.argv[3]

	testbed_old = Catalog.from_directory(testbed_old)
	testbed_new = Catalog.from_directory(testbed_new)
	full_old = Catalog.from_directory(full_old)
	
	testbed_suite = GreenSuite()
	full_suite = GreenSuite()

	filenames = testbed_suite.run(testbed_old, testbed_new)
	filenames = full_suite.run(full_old)

	print
	print

	extrapolation_suite = Extrapolation()
	extrapolation_suite.run(testbed_suite, full_suite)






if __name__ == '__main__':
	suite_extrapolate()
	# suite_compare()