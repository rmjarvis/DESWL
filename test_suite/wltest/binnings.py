import numpy as np
from . import lazy_pylab as pylab

class BinnedTrendMethods(object):
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
    def binned_mean_equal_count_fit(x, y, n, **plot_args):
        bins = BinnedTrendMethods.find_equal_count_bins(x, n)
        bin = np.digitize(x, bins) - 1
        y_mean, y_std = BinnedTrendMethods.binned_means(y, bin, n)
        x_mid = (bins[1:]+bins[:-1])/2.
        return np.polyfit(x_mid, y_mean, 1)


    @staticmethod
    def binned_mean_equal_count_plot(x, y, n, **plot_args):
        bins = BinnedTrendMethods.find_equal_count_bins(x, n)
        bin = np.digitize(x, bins) - 1
        y_mean, y_std = BinnedTrendMethods.binned_means(y, bin, n)
        x_mid = (bins[1:]+bins[:-1])/2.
        pylab.errorbar(x_mid, y_mean, y_std, fmt='.', **plot_args)
        return np.polyfit(x_mid, y_mean, 1), x_mid, y_mean, y_std 

    @staticmethod
    def binned_mean_equal_width_plot(x, y, n, **plot_args):
        bins = np.linspace(x.min(), x.max(), n+1)
        bin = np.digitize(x, bins) - 1
        y_mean, y_std = BinnedTrendMethods.binned_means(y, bin, n)
        x_mid = (bins[1:]+bins[:-1])/2.
        pylab.errorbar(x_mid, y_mean, y_std, fmt='.', **plot_args)
        return np.polyfit(x_mid, y_mean, 1)  
