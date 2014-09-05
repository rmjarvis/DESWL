from .test_base import SingleCatalogTest
import numpy as np
from . import lazy_pylab as pylab

class MeanE(SingleCatalogTest):
    name = "mean_e"
    statistic_target = 0.001

    def run(self, cat):
        stat = np.array([cat['e1'].mean(), cat['e2'].mean()])
        return stat

class BinnedTrend(SingleCatalogTest):
    @staticmethod
    def binned_means(y, bin, n):
        y_mean = []
        y_std = []
        for i in xrange(n):
            w = (bin==i)
            y_bin = y[w]
            y_mean.append(y_bin.mean())
            y_std.append(y_bin.std() / y_bin.size**0.5)
        y_mean = np.array(y_mean)
        y_std = np.array(y_std)
        return y_mean, y_std

    @staticmethod
    def find_equal_count_bins(x,n):
        xs = np.sort(x)
        m = len(xs)
        r = np.linspace(0.0,1.0,n+1) * (m-1)
        return xs[r.astype(int)]


    @staticmethod
    def binned_mean_equal_count_plot(x, y, n, **plot_args):
        bins = BinnedTrend.find_equal_count_bins(x, n)
        bin = np.digitize(x, bins) - 1
        y_mean, y_std = BinnedTrend.binned_means(y, bin, n)
        x_mid = (bins[1:]+bins[:-1])/2.
        pylab.errorbar(x_mid, y_mean, y_std, fmt='.', **plot_args)
        return np.polyfit(x_mid, y_mean, 1)

    def binned_mean_equal_width_plot(x, y, n, **plot_args):
        bins = np.linspace(x.min(), x.max(), n+1)
        bin = np.digitize(x, bins) - 1
        y_mean, y_std = BinnedTrend.binned_means(y, bin, n)
        x_mid = (bins[1:]+bins[:-1])/2.
        pylab.errorbar(x_mid, y_mean, y_std, fmt='.', **plot_args)
        return np.polyfit(x_mid, y_mean, 1)  

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
        return np.abs(np.array(outputs))


class EWithSNR(BinnedTrend):
    name = "e_with_snr"
    x_axes = ["snr"]
    y_axes = ["e1", "e2"]
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
