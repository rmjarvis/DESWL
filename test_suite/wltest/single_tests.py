from .test_base import SingleCatalogTest
from .binnings import BinnedTrendMethods
import numpy as np
from . import lazy_pylab as pylab

class MeanE(SingleCatalogTest):
    name = "mean_e"
    statistic_names = ["mean_e1", "mean_e2"]
    statistic_target = 0.001

    def run(self, cat):
        stat = np.array([cat['e1'].mean(), cat['e2'].mean()])
        return stat

class HighSNRMeanE(SingleCatalogTest):
    name = "mean_e_high_snr"
    statistic_names = ["mean_e1", "mean_e2"]    
    statistic_target = 0.001

    def run(self, cat):
        high = cat['snr']>100
        stat = np.array([cat['e1'][high].mean(), cat['e2'][high].mean()])
        return stat


class BinnedTrend(SingleCatalogTest, BinnedTrendMethods):
    x_axes = []
    y_axes = []
    xlog=False
    ylog=False
    n = 10

    def run(self, cat):
        outputs = []
        for x_axis in self.x_axes:
            for y_axis in self.y_axes:
                filename = self.filename("%s_vs_%s"%(y_axis, x_axis))
                self.figure(filename)
                x = cat[x_axis]
                y = cat[y_axis]
                p = self.binned_mean_equal_count_plot(x, y, self.n, label=cat.name)
                pylab.xlabel(x_axis)
                pylab.ylabel(y_axis)
                X = [x.min(), x.max()]
                pylab.plot(X, np.polyval(p,X),label=cat.name+" fit")
                outputs.extend(p)
                pylab.legend(loc='lower right')
        return np.array(outputs)


class EWithSNR(BinnedTrend):
    name = "e_with_snr"
    x_axes = ["snr"]
    y_axes = ["e1", "e2"]
    statistic_names = ["c1", "m1","c2","m2"]    
    n = 10

    #Make the x-axis logarithmic after plotting
    def run(self, cat):
        r = super(EWithSNR,self).run(cat)
        for x_axis in self.x_axes:
            for y_axis in self.y_axes:
                filename = self.filename("%s_vs_%s"%(y_axis, x_axis))
                self.figure(filename)
                pylab.xscale("log")
        return r

    statistic_target = np.array([
        0.01, 0.01, 0.01, 0.01 ])
        #c1    m1    c2     m2



class EWithPSF(BinnedTrend):
    n = 10
    statistic_target = np.array([
        0.01, #c1_1
        0.01, #m1_1
    ])


class EWithPSF11(EWithPSF):
    name = "e1_with_psf1"
    x_axes = ["mean_psf_e1_sky"]
    y_axes = ["e1"]

class EWithPSF12(EWithPSF):
    name = "e2_with_psf1"
    x_axes = ["mean_psf_e1_sky"]
    y_axes = ["e2"]

class EWithPSF21(EWithPSF):
    name = "e1_with_psf2"
    x_axes = ["mean_psf_e2_sky"]
    y_axes = ["e1"]

class EWithPSF22(EWithPSF):
    name = "e2_with_psf2"
    x_axes = ["mean_psf_e2_sky"]
    y_axes = ["e2"]



