#coding: utf-8
from . import test_base
from . import single_tests
from . import pair_tests

from glob import glob
from .catalogs import Catalog

from . import termcolor

TICK = termcolor.colored(u"PASS" ,color='green')
CROSS = termcolor.colored(u"FAIL" ,color='red')
CHANGE = termcolor.colored(u"DIFFERENT" ,color='blue')
NO_CHANGE = termcolor.colored(u"SAME" ,color='green')

class Suite(object):
	classes = [
	]
	intersect=True


	def __init__(self):
		options = {}
		self.tests = [cls(**options) for cls in self.classes]

	def select(self, cat):
		return cat

	def run_single_test(self, test, cat):
		if cat is None: return
		passed = test(cat)
		if passed: print u" - %s  %s  %s  (%s<%s)"% (
			cat.name, test.name, TICK, test.statistic, test.statistic_target)
		else: print u" - %s  %s  %s  (%s>%s)"% (
			cat.name, test.name, CROSS, test.statistic, test.statistic_target)

	def run_pair_test(self, test, cat1, cat2):
		if cat2 is None: return
		passed = test(cat1,cat2)
		if passed: print u" - %s %s (%s<%s)"% (
			test.name, NO_CHANGE, test.statistic, test.statistic_target)
		else: print u" - %s  %s  (%s>%s)"% (
			test.name, CHANGE, test.statistic, test.statistic_target)



	def run(self, cat1, cat2=None):
		cat1 = self.select(cat1)
		if cat2 is not None:
			print "Testing %s vs %s" % (cat1.name, cat2.name)
			cat2 = self.select(cat2)
			if self.intersect:
				print " - Intersecting catalogs"
				print
				cat1, cat2 = cat1.intersection(cat2)
		else:
			print "Testing %s alone" % cat1.name
		filenames = []
		for test in self.tests:
			if isinstance(test, test_base.SingleCatalogTest):
				self.run_single_test(test, cat1)
				self.run_single_test(test, cat2)
			elif isinstance(test, test_base.PairCatalogTest):
				self.run_pair_test(test, cat1,cat2)
			print
			filenames.extend(test.figures.keys())
			test.save_figures()

		return filenames
