from .test_base import PairCatalogTest
import numpy as np
from . import lazy_pylab as pylab

class EDifferenceTest(PairCatalogTest):
	name = "e_difference"
	statistic_target = 0.001
	def run(self, cat1, cat2):
		d1 = (cat1['e1'] - cat2['e1']).mean()
		d2 = (cat1['e2'] - cat2['e2']).mean()
		filename = self.filename("e1_difference")
		fig = self.figure(filename)
		pylab.plot(cat1['e1'], cat2['e1']-cat1['e1'],',')
		pylab.xlabel("$e^A_1$")
		pylab.ylabel("$e^B_1 - e^A_1$")

		filename = self.filename("e2_difference")
		fig = self.figure(filename)
		pylab.plot(cat1['e2'], cat2['e2']-cat1['e2'],',')
		pylab.xlabel("$e^A_2$")
		pylab.ylabel("$e^B_2 - e^A_2$")

		return np.array([d1,d2])