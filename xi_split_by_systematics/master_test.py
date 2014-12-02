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


sys.exit()

misc_funcs.split_systematics.split_gals_lin_along_(i3.psf2,i3,'psf2')


#=======================================================
# ignore below this point #
#=======================================================

theta,ggp,ggm,ggperr,ggmerr=misc_funcs.xi_2pt_shear_test_methods.xi_2pt(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,bs,wt,regs,num_reg)

print ' '
print ' '
print '-----------------------'
print 'b/s - '+str(bs)+', w - '+str(wt)
print '-----------------------'
mean_e1,std_e1,rms_e1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e1,m,c1,s1,w,cat,bs,wt)
mean_e2,std_e2,rms_e2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e2,m,c2,s2,w,cat,bs,wt)
print 'mean e'
print 'e1',mean_e1
print 'e2',mean_e2
print 'sigma e'
print 'e1',std_e1
print 'e2',std_e2
print 'rms e'
print 'e1',rms_e1
print 'e2',rms_e2
mean_psf1,std_psf1,rms_psf1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(psf1,m,c1,s1,w,cat,False,False)
mean_psf2,std_psf2,rms_psf2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(psf2,m,c2,s2,w,cat,False,False)
print 'mean psf'
print 'psf1',mean_psf1
print 'psf2',mean_psf2
print 'sigma psf'
print 'psf1',std_psf1
print 'psf2',std_psf2
print 'rms psf'
print 'psf1',rms_psf1
print 'psf2',rms_psf2
if misc_funcs.CatalogMethods.cat_name(cat)=='im3shapev7':
  mean_m,std_m,rms_m=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(m,m,c2,s1,w,cat,False,False)
  mean_c1,std_c1,rms_c1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(c1,m,c1,s1,w,cat,False,False)
  mean_c2,std_c2,rms_c2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(c2,m,c2,s2,w,cat,False,False)
  mean_s1=0
  mean_s2=0
  std_s1=0
  std_s2=0
  rms_s1=0
  rms_s2=0
  print 'mean bias params'
  print 'm',mean_m
  print 'c1',mean_c1
  print 'c2',mean_c2
  print 'sigma bias params'
  print 'm',std_m
  print 'c1',std_c1
  print 'c2',std_c2
  print 'rms bias params'
  print 'm',rms_m
  print 'c1',rms_c1
  print 'c2',rms_c2
if misc_funcs.CatalogMethods.cat_name(cat)=='ngmix009':
  mean_s1,std_s1,rms_s1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(s1,m,c1,s1,w,cat,False,False)
  mean_s2,std_s2,rms_s2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(s2,m,c2,s2,w,cat,False,False)
  mean_m=0
  mean_c1=0
  mean_c2=0
  std_m=0
  std_c1=0
  std_c2=0
  rms_m=0
  rms_c1=0
  rms_c2=0
  print 'mean sensitivity'
  print 's1',mean_s1
  print 's2',mean_s2
  print 'sigma sensitivity'
  print 's1',std_s1
  print 's2',std_s2
  print 'rms sensitivity'
  print 's1',rms_s1
  print 's2',rms_s2     
mean_w,std_w,rms_w=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(w,m,c2,s1,w,cat,False,False)
print 'mean weight'
print 'w',mean_w
print 'sigma weight'
print 'w',std_w
print 'rms weight'
print 'w',rms_w
misc_funcs.split_systematics.split_gals_lin_along_(w,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'weight')
misc_funcs.split_systematics.split_gals_lin_along_(psf1,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'psf1')
misc_funcs.split_systematics.split_gals_2pt_along_(psf1,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'psf1')
misc_funcs.split_systematics.split_gals_lin_along_(psf2,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'psf2')
misc_funcs.split_systematics.split_gals_2pt_along_(psf2,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'psf2')
misc_funcs.split_systematics.split_gals_lin_along_(fwhm,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'psf size')
misc_funcs.split_systematics.split_gals_2pt_along_(fwhm,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'fwhm')
print ' '
theta,ggp,ggm,ggperr,ggmerr=misc_funcs.xi_2pt_shear_test_methods.xi_2pt(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,bs,wt,regs,num_reg)
gpp,gpm,gpperr,ppp,ppm,ppperr,alpha,alphaerr,leakage,leakageerr=misc_funcs.xi_2pt_shear_test_methods.xi_2pt_psf(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,num_reg)
label='b/s - '+str(bs)+', w - '+str(wt)
misc_funcs.plotting_methods.fig_open_xi_2pt(theta,ggp,gpp,ppp,ggperr,gpperr,ppperr,alpha,alphaerr,leakage,leakageerr,bins,cat,bs,wt,label)
#misc_funcs.xi_2pt_shear_test_methods.save_xi_2pt(theta,ggp,ggm,ggperr,gpp,gpm,gpperr,ppp,ppm,ppperr,alpha,alphaerr,bins,cat,bs,wt)
mean_g1,std_g1,rms_g1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e1-alpha[len(alpha)/2]*psf1,m,c1,s1,w,cat,bs,wt)
mean_g2,std_g2,rms_g2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e2-alpha[len(alpha)/2]*psf2,m,c2,s2,w,cat,bs,wt)
print 'mean g'
print 'g1',mean_g1
print 'g2',mean_g2
print 'sigma g'
print 'g1',std_g1
print 'g2',std_g2
print 'rms g'
print 'g1',rms_g1
print 'g2',rms_g2

misc_funcs.plotting_methods.fig_close_xi_2pt(bins,cat)


    # misc_funcs.split_systematics.split_gals_lin_along_(ebv,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'ebv')
    # misc_funcs.split_systematics.split_gals_lin_along_(ebverr,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'ebverr')
    # misc_funcs.split_systematics.split_gals_lin_along_(ra,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'ra')
    # misc_funcs.split_systematics.split_gals_lin_along_(dec,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'dec')
    # misc_funcs.split_systematics.split_gals_lin_along_(exptime,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'exptime')
    # misc_funcs.split_systematics.split_gals_lin_along_(maglimit,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'maglimit')
    # misc_funcs.split_systematics.split_gals_lin_along_(skysigma,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'skysigma')
    # misc_funcs.split_systematics.split_gals_lin_along_(skybrite,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'skybrite')
    # misc_funcs.split_systematics.split_gals_lin_along_(airmass,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'airmass')
    # misc_funcs.split_systematics.split_gals_lin_along_(fwhm,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'fwhm')
    # misc_funcs.split_systematics.split_gals_lin_along_(skysigmacoadd,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'skysigmacoadd')
    # misc_funcs.split_systematics.split_gals_lin_along_(skybritecoadd,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'skybritecoadd')
    # misc_funcs.split_systematics.split_gals_lin_along_(airmasscoadd,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'airmasscoadd')
    # misc_funcs.split_systematics.split_gals_lin_along_(fwhmcoadd,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,num_reg,'fwhmcoadd')
misc_funcs.split_systematics.split_gals_2pt_along_(snr,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'SNR')
misc_funcs.split_systematics.split_gals_2pt_along_(ebv,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'ebv')
misc_funcs.split_systematics.split_gals_2pt_along_(ebverr,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'ebverr')
misc_funcs.split_systematics.split_gals_2pt_along_(ra,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'ra')
misc_funcs.split_systematics.split_gals_2pt_along_(dec,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'dec')
misc_funcs.split_systematics.split_gals_2pt_along_(skybrite,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'skybrite')
misc_funcs.split_systematics.split_gals_2pt_along_(exptime,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'exptime')
misc_funcs.split_systematics.split_gals_2pt_along_(maglimit,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'maglimit')
misc_funcs.split_systematics.split_gals_2pt_along_(skysigma,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'skysigma')
misc_funcs.split_systematics.split_gals_2pt_along_(airmass,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'airmass')
misc_funcs.split_systematics.split_gals_2pt_along_(fwhm2,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'fwhm2')
misc_funcs.split_systematics.split_gals_2pt_along_(skysigmacoadd,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'skysigmacoadd')
misc_funcs.split_systematics.split_gals_2pt_along_(skybritecoadd,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'skybritecoadd')
misc_funcs.split_systematics.split_gals_2pt_along_(airmasscoadd,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'airmasscoadd')
misc_funcs.split_systematics.split_gals_2pt_along_(fwhmcoadd,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'fwhmcoadd')
misc_funcs.split_systematics.split_gals_2pt_along_(psf1,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'psf1')
misc_funcs.split_systematics.split_gals_2pt_along_(psf2,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'psf2')
misc_funcs.split_systematics.split_gals_2pt_along_(fwhm,zpeak,ra,dec,e1,e2,nbin,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,num_reg,'fwhm')


    # theta,ggp,ggm,ggperr,ggmerr=misc_funcs.xi_2pt_shear_test_methods.xi_2pt(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,bs,wt,regs,num_reg)
    # gpp,gpm,gpperr,ppp,ppm,ppperr,alpha,alphaerr,leakage,leakageerr=misc_funcs.xi_2pt_shear_test_methods.xi_2pt_psf(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,num_reg)
    # label='b/s - '+str(bs)+', w - '+str(wt)
    # misc_funcs.plotting_methods.fig_open_xi_2pt(theta,ggp,gpp,ppp,ggperr,gpperr,ppperr,alpha,alphaerr,leakage,leakageerr,bins,cat,bs,wt,label)
    # mean_g1,std_g1,rms_g1=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e1-c1-alpha[len(alpha)/2]*psf1,m,0,s1,w,cat,bs,wt)
    # mean_g2,std_g2,rms_g2=misc_funcs.lin_shear_test_methods.calc_mean_std_rms_e(e2-c2-alpha[len(alpha)/2]*psf2,m,0,s2,w,cat,bs,wt)
    # print 'mean g'
    # print 'g1',mean_g1
    # print 'g2',mean_g2
    # print 'sigma g'
    # print 'g1',std_g1
    # print 'g2',std_g2
    # print 'rms g'
    # print 'g1',rms_g1
    # print 'g2',rms_g2
    # #misc_funcs.xi_2pt_shear_test_methods.save_xi_2pt(theta,ggp,ggm,ggperr,gpp,gpm,gpperr,ppp,ppm,ppperr,alpha,alphaerr,bins,cat,bs,wt)

misc_funcs.plotting_methods.fig_close_xi_2pt(bins,cat)

