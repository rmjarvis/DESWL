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
import fitsio as fio
import kmeans_radec

noval=999999

cat=misc_funcs.CatalogMethods.cat_num('im3shapev7')
i3=misc_funcs.CatalogStore(cat,True)

cat=misc_funcs.CatalogMethods.cat_num('ngmix010')
ng=misc_funcs.CatalogStore(cat,True)

i3.use_zrw=True
i3.ztyp='DESDM'
ng.use_zrw=True
ng.ztyp='DESDM'

fits=fio.FITS('sva1_gold_1.0_catalog_desdmphotoz.fits')
tmp=fits[-1].read()

mask1=np.in1d(ng.coadd,tmp['COADD_OBJECTS_ID'],assume_unique=True)
mask2=np.in1d(tmp['COADD_OBJECTS_ID'],ng.coadd,assume_unique=True)
tmp2=tmp[mask2]
sort=np.argsort(tmp2['COADD_OBJECTS_ID'])[np.argsort(np.argsort(ng.coadd[mask1]))]
tmp2=tmp2[sort]
ng.zpeak=tmp2['ZP']

mask1=np.in1d(i3.coadd,tmp['COADD_OBJECTS_ID'],assume_unique=True)
mask2=np.in1d(tmp['COADD_OBJECTS_ID'],i3.coadd,assume_unique=True)
tmp2=tmp[mask2]
sort=np.argsort(tmp2['COADD_OBJECTS_ID'])[np.argsort(np.argsort(i3.coadd[mask1]))]
tmp2=tmp2[sort]
i3.zpeak=tmp2['ZP']

misc_funcs.SVA1.Fig1_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig2_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig3_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig4_2ptPaper(i3,ng)
misc_funcs.SVA1.Fig5_2ptPaper(i3,ng)

misc_funcs.SVA1.Fig1_ShearPaper(i3,ng)
misc_funcs.SVA1.Fig2_ShearPaper(i3,ng)




