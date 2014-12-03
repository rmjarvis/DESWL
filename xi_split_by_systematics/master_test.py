import numpy as np
import numpy.ma as ma
import treecorr
import sys
import pyfits
import os
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
#plt.style.use('SVA1StyleSheet.mplstyle')
#plt.minorticks_on()
#plt.tight_layout()
import pylab
import numpy.linalg as linalg
import argparse
import scipy.interpolate as interp
import scipy.integrate as scint
import numpy.random as npr
import pandas as pd
import misc_funcs

noval=999999


cat=misc_funcs.CatalogMethods.cat_num('im3shapev7')
i3=misc_funcs.CatalogStore(cat,True)

cat=misc_funcs.CatalogMethods.cat_num('ngmix010')
ng=misc_funcs.CatalogStore(cat,True)

misc_funcs.SVA1.Fig1_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig2_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig3_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig4_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig5_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig1_ShearPaper(i3,ng)

#misc_funcs.SVA1.Fig1_ShearPaper(i3,ng)


