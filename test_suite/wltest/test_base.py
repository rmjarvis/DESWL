#coding: utf-8
import numpy as np
from . import lazy_pylab as pylab

class BaseCatalogTest(object):
    statistic_target = np.nan
    name = "Unknown"

    def __init__(self, file_type='png', output_dir='.', prefix=''):
        self.target = self.statistic_target
        self.statistic = np.nan
        self.plot = None

        self.file_type = file_type
        self.output_dir = output_dir
        self.prefix = prefix
        if self.prefix: self.prefix=self.prefix + "_"
        self.figures = {}

    def passed(self):
        return np.all(self.statistic < self.target)

    def filename(self, base, *bases):
        if bases:
            base = base + "_" + ("_".join(bases))
        return "{0}/{1}{2}.{3}".format(self.output_dir, self.prefix, base, self.file_type)

    def figure(self, name):
        #we want to be able to plot multiple chains on the same
        #figure at some point.  So when we make figures we should
        #call this function
        fig = self.figures.get(name)
        if fig is None:
            fig = pylab.figure()
            self.figures[name] = fig
        else:
            pylab.figure(fig.number)
        return fig

    def save_figures(self):
        for filename, figure in self.figures.items():
            pylab.figure(figure.number)
            pylab.savefig(filename)
            pylab.close()
        self.figures = {}       

class SingleCatalogTest(BaseCatalogTest):
    def __call__(self, cat):
        self.statistic = self.run(cat)
        return self.passed()

    def run(self, cat):
        # maybe make a plot - use self.filename and self.figure
        # compute and return statistic
        return np.nan


class PairCatalogTest(BaseCatalogTest):
    def __call__(self, cat1, cat2):
        self.statistic = self.run(cat1, cat2)
        return self.passed()

    def run(self, cat1, cat2):
        # maybe make a plot - use self.filename and self.figure
        # compute and return statistic
        return np.nan
