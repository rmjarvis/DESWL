#coding: utf-8
from .test_base import PairCatalogTest
from .binnings import BinnedTrendMethods
import numpy as np
from . import lazy_pylab as pylab

class EDifferenceTest(PairCatalogTest):
	name = "e_difference"
	statistic_target = 0.001
	def run(self, cat1, cat2):
		d1 = (cat1['e1'] - cat2['e1']).mean()
		d2 = (cat1['e2'] - cat2['e2']).mean()
		d1_error = (cat1['e1'] - cat2['e1']).std()/len(cat1['e1'])**0.5
		d2_error = (cat1['e2'] - cat2['e2']).std()/len(cat2['e2'])**0.5

		print "(difference values d1 = %e ± %e)" % (d1, d1_error)
		print "(difference values d2 = %e ± %e)" % (d2, d2_error)
		print
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

class PSFDifferenceTest11(PairCatalogTest, BinnedTrendMethods):
	name = "psf_m_difference_11"
	statistic_target = [0.001, 0.01]
	def run(self, cat1, cat2):
		x_axis = "mean_psf_e1_sky"
		y_axis = "e1"
		x1 = cat1[x_axis]
		y1 = cat1[y_axis]
		x2 = cat2[x_axis]
		y2 = cat2[y_axis]
		n = 10
		p1 = BinnedTrendMethods.binned_mean_equal_count_fit(x1, y1, n)
		p2 = BinnedTrendMethods.binned_mean_equal_count_fit(x2, y2, n)

		m1 = p1[0]
		m2 = p2[0]
		c1 = p1[1]
		c2 = p2[1]
		print "(difference values m11_A = %e)" % (m1, )
		print "(difference values m11_B = %e)" % (m2, )
		print "(difference values c11_A = %e)" % (c1, )
		print "(difference values c11_B = %e)" % (c2, )
		print

		dm = m2-m1
		dc = c2-c1

		return np.array([dm,dc])

	
