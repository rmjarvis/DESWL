from .catalogs import Catalog
from .suite import YellowSuite
import sys

def suite_compare():
	dirname1 = 	sys.argv[1]
	dirname2 = 	sys.argv[2]

	cat1 = Catalog.from_directory(dirname1)
	cat2 = Catalog.from_directory(dirname2)
	suite = YellowSuite()
	print
	print "Running tests"
	print
	filenames = suite.run(cat1, cat2)
	print
	if filenames:
		print "Made files:"
		for filename in filenames:
			print " - ", filename


if __name__ == '__main__':
	suite_compare()