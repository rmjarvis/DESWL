import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
#plt.style.use('SVA1StyleSheet.mplstyle')
import pylab
import treecorr
import healpy as hp
import glob
import pyfits as py
import astropy.table as apt
import multiprocessing as mp
import multiprocessing.sharedctypes as mpc
import homogenise_nz as hnz
import fitsio as fio
import numpy.lib.recfunctions as nlr
import time
import ctypes

  
class CatalogStore(object):

  def __init__(self,coadd,ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,tbins,sbins,lbins,sep,slop,use_jk,bs,wt,regs,num_reg):
    self.coadd=coadd
    self.ra=ra
    self.dec=dec
    self.e1=e1
    self.e2=e2
    self.m=m
    self.c1=c1
    self.c2=c2
    self.s1=s1
    self.s2=s2
    self.w=w
    self.regs=regs
    self.cat=cat
    self.tbins=tbins
    self.sbins=sbins
    self.lbins=lbins
    self.sep=sep
    self.slop=slop
    self.use_jk=use_jk
    self.bs=bs
    self.wt=wt
    self.num_reg=num_reg


class lin_shear_test_methods(object):

  @staticmethod
  def get_lin_e_w_norm(e,m,c,s,w,cat,bs,wt):

    if CatalogMethods.cat_name(cat)=='im3shapev7':

      if bs & wt:
        mw=w
        me=(e-c)
        norm=w*(m+1)
      elif bs and not wt:
        mw=np.ones((len(e)))
        me=e-c
        norm=m+1
      elif wt and not bs:
        mw=w
        me=e
        norm=w
      else:
        mw=np.ones((len(e)))
        me=e
        norm=np.ones((len(e)))

    elif CatalogMethods.cat_name(cat)=='ngmix009':

      if bs & wt:
        mw=w
        me=e
        norm=w*s
      elif bs and not wt:
        mw=np.ones((len(e)))
        me=e
        norm=s
      elif wt and not bs:
        mw=w
        me=e
        norm=w
      else:
        mw=np.ones((len(e)))
        me=e
        norm=np.ones((len(e)))

    return me,mw,norm

  @staticmethod
  def calc_mean_std_rms_e(e,m,c,s,w,cat,bs,wt):

    me,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(e,m,c,s,w,cat,bs,wt)

    mean=np.sum(mw*me)/np.sum(norm)
    std=(np.sum(mw*(me-mean)**2)/np.sum(norm))
    rms=np.sqrt(np.sum((mw*me)**2)/np.sum(mw**2))

    return mean,std,rms

  @staticmethod
  def find_bin_edges(x,nbins):

    xs=np.sort(x)
    m=len(xs)
    r=np.linspace(0.,1.,nbins+1.)*(m-1)

    return xs[r.astype(int)]

  @staticmethod
  def binned_means(y,bin,nbins,m,c,s,w,cat,bs,wt):

    y_mean=[]
    y_std=[]

    for i in xrange(nbins):
      mask=(bin==i)
      mean,std,rms=lin_shear_test_methods.calc_mean_std_rms_e(y[mask],m[mask],c[mask],s[mask],w[mask],cat,bs,wt)
      y_mean.append(mean)
      y_std.append(std/np.sqrt(len(y[mask])))

    y_mean=np.array(y_mean)
    y_std=np.array(y_std)

    return y_mean,y_std

  @staticmethod
  def bin_means(y,x,nbins,m,c,s,w,cat,bs,wt):

    edge=lin_shear_test_methods.find_bin_edges(x,nbins)
    xbin=np.digitize(x,edge)-1

    y_mean,y_std=lin_shear_test_methods.binned_means(y,xbin,nbins,m,c,s,w,cat,bs,wt)
    x_mean,x_std=lin_shear_test_methods.binned_means(x,xbin,nbins,m,c,s,w,cat,False,False)

    return x_mean,x_std,y_mean,y_std

  @staticmethod
  def ngmix_weight_calc(cov11,cov22,cov12):
    sn=0.16

    w=1./(2.*sn**2.+cov11+cov22+2.*cov12)
    print w[w<0]

    return w

class xi_2pt_shear_test_methods(object):

  @staticmethod
  def xi_2pt(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,bs,wt,regs,reg_num):


    me1,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(e1,m,c1,np.sqrt(s1*s2),w,cat,bs,wt)
    me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(e2,m,c2,np.sqrt(s1*s2),w,cat,bs,wt)

    gg = treecorr.GGCorrelation(nbins=bins, min_sep=sep[0], max_sep=sep[1], sep_units='arcmin',binslop=slop,verbose=0)
    kk = treecorr.KKCorrelation(nbins=bins, min_sep=sep[0], max_sep=sep[1], sep_units='arcmin',binslop=slop,verbose=0)

    cate=treecorr.Catalog(g1=mw*me1, g2=mw*me2, ra=ra, dec=dec, ra_units='deg', dec_units='deg')
    catm=treecorr.Catalog(k=norm, ra=ra, dec=dec, ra_units='deg', dec_units='deg')

    kk.process(catm,catm)
    kkp = np.copy(kk.xi)

    gg.process(cate)
    ggp = np.copy(gg.xip)/kkp
    ggm = np.copy(gg.xim)/kkp
    ggperr = ggmerr = np.sqrt(np.copy(gg.varxi))/kkp

    theta = np.exp(gg.meanlogr)
    print 'xip',ggp

    if use_jk:
      ggperr,ggmerr=jackknife_methods.xi_2pt_err0(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,bs,wt,regs,reg_num)

    return theta,ggp,ggm,ggperr,ggmerr

  @staticmethod
  def xi_2pt_psf(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,use_jk,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,reg_num):

    me1,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(e1,m,c1,np.sqrt(s1*s2),w,cat,bs,wt)
    me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(e2,m,c2,np.sqrt(s1*s2),w,cat,bs,wt)

    gg = treecorr.GGCorrelation(nbins=bins, min_sep=sep[0], max_sep=sep[1], sep_units='arcmin',binslop=slop,verbose=0)
    nk = treecorr.NKCorrelation(nbins=bins, min_sep=sep[0], max_sep=sep[1], sep_units='arcmin',binslop=slop,verbose=0)

    cate=treecorr.Catalog(g1=mw*me1, g2=mw*me2, ra=ra, dec=dec, ra_units='deg', dec_units='deg')
    catpsf=treecorr.Catalog(g1=psf1, g2=psf2, ra=ra, dec=dec, ra_units='deg', dec_units='deg')
    catm=treecorr.Catalog(k=norm, ra=ra, dec=dec, ra_units='deg', dec_units='deg')
    catn=treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg')

    nk.process(catn,catm)
    nkp = np.copy(nk.xi)

    gg.process(cate,catpsf)
    gpp = np.copy(gg.xip)/nkp
    gpm = np.copy(gg.xim)/nkp
    gpperr = np.sqrt(np.copy(gg.varxi))/nkp

    gg.process(catpsf)
    ppp=np.copy(gg.xip)
    ppm=np.copy(gg.xim)
    ppperr=np.sqrt(np.copy(gg.varxi))

    alpha=xi_2pt_shear_test_methods.calc_alpha(gpp,ppp,mean_e1,mean_e2,mean_psf1,mean_psf2)
    alphaerr=0.
    print 'alpha = ',alpha[bins/2]

    leakage=xi_2pt_shear_test_methods.calc_psf_leakage(alpha,ppp,mean_e1,mean_e2,mean_psf1,mean_psf2)
    leakage/=ggp
    leakageerr=0.
    print 'psf leakage = ',leakage[bins/2]

    theta = np.exp(gg.meanlogr)

    if use_jk:
      gpperr,gpmerr,ppperr,ppmerr,alphaerr,leakageerr=jackknife_methods.xi_2pt_psf_err(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,reg_num)

    return gpp,gpm,gpperr,ppp,ppm,ppperr,alpha,alphaerr,leakage,leakageerr

  @staticmethod
  def calc_alpha(gp,pp,e1,e2,psf1,psf2):

    alpha=(gp-e1*psf1-e2*psf2)/(pp-psf1**2-psf2**2)

    return alpha

  @staticmethod
  def calc_psf_leakage(alpha,pp,e1,e2,psf1,psf2):

    leakage=alpha**2.*(pp-2.*(psf1**2.+psf2**2.))+2.*alpha*(e1*psf1+e2*psf2)

    return leakage

  @staticmethod
  def save_xi_2pt(theta,ggp,ggm,ggv,gpp,gpm,gpv,ppp,ppm,ppv,alpha,alphaerr,nbins,cat,bs,wt):

    tmp=np.column_stack((theta,ggp,ggm,ggv,gpp,gpm,gpv,ppp,ppm,ppv,alpha,alphaerr))

    np.savetxt(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_2pt_xi.dat',tmp)
    np.save(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_2pt_xi.npy',tmp)

    return

  @staticmethod
  def save_xi_2pt0(tmp,nbins,cat,bs,wt,label):

    np.savetxt(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_2pt_'+label+'.dat',tmp)
    np.save(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_2pt_'+label+'.npy',tmp)

    return

class split_systematics(object):

  @staticmethod
  def init_systematics_maps(ra,dec):

    dir='/share/des/sv/systematics_maps/'

    ebv=split_systematics.load_sys_map_to_array(ra,dec,dir+'Planck_EBV_2048r_Q.fits',False,2048,True)
    ebverr=split_systematics.load_sys_map_to_array(ra,dec,dir+'Planck_EBVerr_2048r_Q.fits',False,2048,True)
    exptime=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    maglimit=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    skysigmacoadd=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA__mean.fits.gz',False,4096,False)
    skybritecoadd=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE__mean.fits.gz',False,4096,False)
    airmasscoadd=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS__mean.fits.gz',False,4096,False)
    fwhmcoadd=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM__mean.fits.gz',False,4096,False)

    return ebv,ebverr,exptime,maglimit,skysigmacoadd,skysigma,skybritecoadd,skybrite,airmasscoadd,airmass,fwhmcoadd,fwhm


  @staticmethod
  def load_sys_map_to_array(ra,dec,sys_file,nested,nside,map):

    print sys_file
    if map:
      sys = hp.read_map(sys_file)
      array = sys[hp.ang2pix(nside, np.pi/2.-np.radians(dec),np.radians(ra), nest=nested)]
    else:
      sys=py.open(sys_file)
      sys=sys[1].data
      sysmap = np.zeros(12*4096**2)
      sysmap[sys['pixel']] = sys['signal']
      array = sysmap[hp.ang2pix(nside, np.pi/2.-np.radians(dec),np.radians(ra), nest=nested)]

    return array

  @staticmethod
  def split_gals_lin_along_(array,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,use_jk,regs,reg_num,label):

    arr1,tmp,me1,e1err=lin_shear_test_methods.bin_means(e1,array,nbins,m,c1,s1,w,cat,bs,wt)
    arr1,tmp,me2,e2err=lin_shear_test_methods.bin_means(e2,array,nbins,m,c2,s2,w,cat,bs,wt)

    slp1,tmp=np.polyfit(arr1, me1, 1)
    slp2,tmp=np.polyfit(arr1, me2, 1)

    if use_jk:
      edge=lin_shear_test_methods.find_bin_edges(array,nbins)
      bin=np.digitize(array,edge)-1
      e1err[i]=jackknife_methods.lin_err0(e1,array,m,c1,s1,w,nbins,cat,bs,wt,regs,reg_num)
      e2err[i]=jackknife_methods.lin_err0(e2,array,m,c2,s2,w,nbins,cat,bs,wt,regs,reg_num)

    plotting_methods.fig_create_e12vs(me1,me2,e1err,e2err,arr1,nbins,cat,bs,wt,label)
    print 'slope of e1 vs '+label,slp1
    print 'slope of e2 vs '+label,slp2

    return

  @staticmethod
  def split_gals_2pt_along_(array,zpeak,ra,dec,e1,e2,nbins,m,c1,c2,s1,s2,w,cat,bs,wt,bins,sep,slop,use_jk,regs,reg_num,label):

    xip=np.zeros((nbins+1,bins))
    xim=np.zeros((nbins+1,bins))
    var=np.zeros((nbins+1,bins))
    xiperr=np.zeros((nbins+1,bins))
    ximerr=np.zeros((nbins+1,bins))
    edge=lin_shear_test_methods.find_bin_edges(array,nbins)
    print edge
    bin=np.digitize(array,edge)-1
    for i in xrange(nbins):
      mask=(bin==i)
      print 'mean splits',np.mean(array[mask])
      list_split=[(zpeak,w),(zpeak[mask],w[mask])]
      list_nz_weights = hnz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3)
      w2=list_nz_weights[1]
      print len(w2),len(w[mask])
      theta,xip[i,:],xim[i,:],xiperr[i,:],tmp=xi_2pt_shear_test_methods.xi_2pt(ra[mask],dec[mask],e1[mask],e2[mask],m[mask],c1[mask],c2[mask],s1[mask],s2[mask],w[mask]*w2,cat,bins,sep,slop,False,bs,wt,regs[mask],reg_num)
      if use_jk:
        xiperr[i,:],ximerr[i,:]=jackknife_methods.xi_2pt_err0(ra[mask],dec[mask],e1[mask],e2[mask],m[mask],c1[mask],c2[mask],s1[mask],s2[mask],w[mask],cat,bins,sep,slop,bs,wt,regs[mask],reg_num)
      plotting_methods.fig_open_xi_2pt0(theta*(1+.02*(i+1)),xip[i,:],xiperr[i,:],nbins,cat,bs,wt,'bin '+str(i+1))

    theta,xip[nbins,:],xim[nbins,:],xiperr[nbins,:],tmp=xi_2pt_shear_test_methods.xi_2pt(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,False,bs,wt,regs,reg_num)
    if use_jk:
      xiperr[nbins,:],ximerr[nbins,:]=jackknife_methods.xi_2pt_err0(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,bs,wt,regs,reg_num)
    plotting_methods.fig_open_xi_2pt0(theta,xip[nbins,:],xiperr[nbins,:],nbins,cat,bs,wt,'full')
    plotting_methods.fig_close_xi_2pt0(nbins,cat,bs,wt,label)

    #tmp=np.column_stack((theta,xip,xiperr))
    #xi_2pt_shear_test_methods.save_xi_2pt0(tmp,nbins,cat,bs,wt,label)

    return 

def xi_2pt_err_slave0(nsi):

  ns=nsi[0]
  i=nsi[1]
  mask=(ns.regs==i)
  print mask
  xip=0
  xim=0
  print 'test'
  if np.sum(mask)>0:
    print ns
    print 'test'
    tmp,xip,xim,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ns.ra[mask],ns.dec[mask],ns.e1[mask],ns.e2[mask],ns.m[mask],ns.c1[mask],ns.c2[mask],ns.s1[mask],ns.s2[mask],ns.w[mask],ns.cat,ns.bins,ns.sep,ns.slop,False,ns.bs,ns.wt,ns.regs,ns.num_reg)

  return xip,xim



class jackknife_methods(object):


  @staticmethod
  def load_jk_regs(ra,dec,healno):

    regfile=np.genfromtxt('pix-list_im3shape-v7_ieflags_base-no_nside'+str(healno)+'-nest.txt')

    regs=np.zeros((len(ra)))
    for i in xrange(len(regfile)):
      print 'jack reg',i
      mask=np.in1d(hp.ang2pix(healno, np.pi/2.-np.radians(dec),np.radians(ra),nest=True),regfile[i][:])
      regs[mask]=i
    
    return len(regfile),regs

  @staticmethod
  def jk_cov(array,nbins,num_reg):

    cov=np.zeros((nbins,nbins))

    for i in xrange(nbins):
      print np.mean(array[:,i]),array[:,i]-np.mean(array[:,i])
      for j in xrange(nbins):
        cov[i,j]=np.sum((array[:,i]-np.mean(array[:,i]))*(array[:,j]-np.mean(array[:,j])))/num_reg*(num_reg-1.)

    return cov

  @staticmethod
  def lin_err0(e,array,m,c,s,w,nbins,cat,bs,wt,regs,num_reg):

    e1=np.zeros((num_reg,nbins))

    for i in xrange(reg_num):
      mask=(regs!=i)
      arr1,tmp,e1[i,:],tmp=lin_shear_test_methods.bin_means(e[mask],array[mask],nbins,m[mask],c[mask],s[mask],w[mask],cat,bs,wt)

    err=np.sqrt(np.diagonal(jackknife_methods.jk_cov(e1,nbins,num_reg)))

    return err


  @staticmethod
  def xi_2pt_err0(ra,dec,e1,e2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,bs,wt,regs,num_reg):


    num_reg=20

    xip=np.zeros((num_reg,bins))
    xim=np.zeros((num_reg,bins))
    tmptime=time.time()

    for i in xrange(num_reg):
      print 'jack_iter', i, time.time()-tmptime
      mask=(regs!=i)
      tmp,xip[i,:],xim[i,:],tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ra[mask],dec[mask],e1[mask],e2[mask],m[mask],c1[mask],c2[mask],s1[mask],s2[mask],w[mask],cat,bins,sep,slop,False,bs,wt,regs,num_reg)
      print 'jack xip',xip[i,:]

    tmp=jackknife_methods.jk_cov(xip,bins,num_reg)
    print 'cov before sqrt',tmp

    xiperr=np.sqrt(np.diagonal(tmp))
    ximerr=0#np.sqrt(np.diagonal(jackknife_methods.jk_cov(xim,bins,num_reg)))

    return xiperr,ximerr

  @staticmethod
  def xi_2pt_psf_err(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,num_reg):

    gpp=np.zeros((num_reg,bins))
    gpm=np.zeros((num_reg,bins))
    ppp=np.zeros((num_reg,bins))
    ppm=np.zeros((num_reg,bins))
    alpha=np.zeros((num_reg,bins))
    leakage=np.zeros((num_reg,bins))

    for i in xrange(reg_num):
      mask=(regs!=i)
      gpp[i,:],gpm[i,:],tmp,ppp[i,:],ppm[i,:],tmp,alpha[i,:],tmp,leakage[i,:],tmp=xi_2pt_shear_test_methods.xi_2pt_psf(ra[mask],dec[mask],e1[mask],e2[mask],psf1[mask],psf2[mask],m[mask],c1[mask],c2[mask],s1[mask],s2[mask],w[mask],cat,bins,sep,slop,False,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,num_reg)

    gpperr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(gpp,bins,num_reg)))
    gpmerr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(gpm,bins,num_reg)))
    ppperr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(ppp,bins,num_reg)))
    ppmerr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(ppm,bins,num_reg)))
    alphaerr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(alpha,bins,num_reg)))
    leakageerr=np.sqrt(np.diagonal(jackknife_methods.jk_cov(leakage,bins,num_reg)))

    return gpperr,gpmerr,ppperr,ppmerr,alphaerr,leakageerr

def task_2pt_function(arg_set):
  tmp,xipi,ximi,tmp,tmp = xi_2pt_shear_test_methods.xi_2pt(*arg_set)
  print 'in parallel loop'
  return xipi,ximi

class plotting_methods(object):

  @staticmethod
  def imshow_symlog(my_matrix, logthresh=3):
    plt.imshow( my_matrix,interpolation='nearest',origin='lower',vmin=np.min(my_matrix), vmax=np.max(my_matrix),norm=matplotlib.colors.SymLogNorm(10**-logthresh) )
    maxlog=int(np.ceil( np.log10(np.max(my_matrix)) ))
    minlog=int(np.ceil( np.log10(-np.min(my_matrix)) ))
    #generate logarithmic ticks 
    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]+[0.0]+[(10**x) for x in xrange(-logthresh,maxlog+1)] )
    plt.colorbar(ticks=tick_locations)
    return 

  @staticmethod
  def fig_create_e12vs(e1,e2,e1err,e2err,x,nbins,cat,bs,wt,label):

    plt.figure(1)
    plt.errorbar(x,e1, yerr=e1err, xerr=None,label='<e1>')
    plt.errorbar(x,e2, yerr=e2err, xerr=None,label='<e2>')
    plt.legend()
    plt.ylabel('<e>')
    plt.xlabel(label)
    plt.ylim((-5e-3,5e-3))
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_eVS'+label+'.png', bbox_inches='tight')
    plt.close(1)

    return


  @staticmethod
  def fig_create_xi_2pt(theta,ggp,gpp,ppp,ggv,gpv,ppv,alpha,nbins,cat,bs,wt):


    plt.figure(1)
    plt.errorbar(theta,ggp, yerr=ggv, xerr=None,label='')
    plt.ylabel(r'$\xi_+(gg)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('symlog')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_xi+gg.png', bbox_inches='tight')
    plt.close(1)

    plt.figure(1)
    plt.errorbar(theta,gpp, yerr=gpv, xerr=None,label='')
    plt.ylabel(r'$\xi_+(gp)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('symlog')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_xi+gp.png', bbox_inches='tight')
    plt.close(1)

    plt.figure(1)
    plt.errorbar(theta,ppp, yerr=ppv, xerr=None,label='')
    plt.ylabel(r'$\xi_+(pp)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('symlog')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_xi+pp.png', bbox_inches='tight')
    plt.close(1)

    plt.figure(1)
    plt.errorbar(theta,alpha, yerr=None, xerr=None,label='')
    plt.ylabel(r'$\alpha$')
    plt.xlabel(r'$\theta$')
    plt.yscale('symlog')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_alpha.png', bbox_inches='tight')
    plt.close(1)

    return

  @staticmethod
  def fig_open_xi_2pt0(theta,xi,err,nbins,cat,bs,wt,label):

    plt.figure(10)
    plt.errorbar(theta,xi, yerr=err, xerr=None,marker='o',linestyle='',label=label)

    return

  @staticmethod
  def fig_close_xi_2pt0(nbins,cat,bs,wt,label):

    plt.figure(10)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\xi_+$')
    plt.xlabel(r'$\theta (arcmin)$')
    plt.xlim((1,200))
    plt.ylim((5e-7,1e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(bs)+'_weight-'+str(wt)+'_xi_'+label+'.png', bbox_inches='tight')
    plt.close(10)

    return



  @staticmethod
  def fig_open_xi_2pt(theta,ggp,gpp,ppp,ggv,gpv,ppv,alpha,alphaerr,leakage,leakageerr,nbins,cat,bs,wt,label):

    plt.figure(5)
    plt.errorbar(theta,ggp, yerr=ggv, xerr=None,label=label)

    plt.figure(6)
    plt.errorbar(theta,np.abs(gpp), yerr=gpv, xerr=None,label=label)

    plt.figure(7)
    plt.errorbar(theta,ppp, yerr=ppv, xerr=None,label=label)

    plt.figure(8)
    plt.errorbar(theta,alpha, yerr=alphaerr, xerr=None,label=label)

    plt.figure(9)
    plt.errorbar(theta,leakage, yerr=leakageerr, xerr=None,label=label)

    return

  @staticmethod
  def fig_close_xi_2pt(nbins,cat):

    plt.figure(5)
    plt.ylabel(r'$\xi_+(gg)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_xi+gg.png', bbox_inches='tight')
    plt.close(5)

    plt.figure(6)
    plt.ylabel(r'$\xi_+(gp)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,1e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_xi+gp.png', bbox_inches='tight')
    plt.close(6)

    plt.figure(7)
    plt.ylabel(r'$\xi_+(pp)$')
    plt.xlabel(r'$\theta$')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-4,3e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_xi+pp.png', bbox_inches='tight')
    plt.close(7)

    plt.figure(8)
    plt.ylabel(r'$\alpha$')
    plt.xlabel(r'$\theta$')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((-.1,.4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_alpha.png', bbox_inches='tight')
    plt.close(8)

    plt.figure(9)
    plt.ylabel(r'$Frac PSF Leakage$')
    plt.xlabel(r'$\theta$')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((-.1,.4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat)+'_nbins-'+str(nbins)+'_leakage.png', bbox_inches='tight')
    plt.close(9)

    return





class ColError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class CatalogMethods(object):


  @staticmethod
  def default_im3shapev7_cuts():
    noval=999999

    cuts=CatalogMethods.add_cut(np.array([]),'error_flag',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'info_flag',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'modest_class',noval,1,noval)
    cuts=CatalogMethods.add_cut(cuts,'gold_match',noval,1,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',56.,noval,94.)
    cuts=CatalogMethods.add_cut(cuts,'dec',noval,noval,-42.)

    return cuts

  @staticmethod
  def default_ngmix009_cuts():
    noval=999999

    cuts=CatalogMethods.add_cut(np.array([]),'exp_flags',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'flags_i',noval,noval,4)
    cuts=CatalogMethods.add_cut(cuts,'modest_class',noval,1,noval)
    cuts=CatalogMethods.add_cut(cuts,'exp_arate',0.4,noval,0.6)
    cuts=CatalogMethods.add_cut(cuts,'exp_s2n_w',10,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'exp_t_s2n',3,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',56.,noval,94.)
    cuts=CatalogMethods.add_cut(cuts,'dec',noval,noval,-42.)

    return cuts

  @staticmethod
  def get_cat_cols(dir,cols,cuts,full,maxiter):
    noval=999999
    len0=0

    for ifile,file in enumerate(glob.glob(dir+'*.fits')):
      if ifile>maxiter:
        break
      print 'file',ifile,file
      try:
        fits=fio.FITS(file)
        #tmparray=apt.Table.read(file)
      except IOError:
        print 'error loading fits file: ',file
        return

      tmparray=fits[-1].read()
      len0+=len(tmparray)

      colex,colist=CatalogMethods.col_exists(cols,fits[-1].get_colnames())#tmparray.dtype.names)
      if colex<1:
        raise ColError('columns '+colist+' do not exist in file: '+file)
        print 'there are column name(s) not in file: ',file
        return

      colex,colist=CatalogMethods.col_exists(cuts['col'],fits[-1].get_colnames())#tmparray.dtype.names)
      if colex<1:
        print 'there are cut column name(s) not in file: ',file
        raise ColError('cut columns '+colist+' do not exist in file: '+file)
        return

      mask=np.array([])
      for icut,cut in enumerate(cuts): 
        mask=CatalogMethods.cuts_on_col(mask,tmparray,cut['col'],cut['min'],cut['eq'],cut['max'])

      tmparray=fits[-1].read(columns=cols)

      #tmparray.keep_columns(cols)

      if ifile==0:
        array=tmparray[mask]
      else:
        array=np.append(array,tmparray[mask],axis=0)

    zpeak,zmean,zmask=CatalogMethods.get_pz(array['coadd_objects_id'])
    zpeak0=np.zeros((len(zpeak),),dtype=[('zpeak','>f8')])
    zmean0=np.zeros((len(zmean),),dtype=[('zmean','>f8')])
    zpeak0['zpeak']=zpeak
    zmean0['zmean']=zmean

    len1=len(array)
    array=array[zmask]
    array=nlr.merge_arrays([array,zpeak0,zmean0],flatten=True)

    print 'selected '+str(len(array))+' galaxies from catalog out of '+str(len0)+' ('+str(len1)+' before pz cut)'

    cols=np.append(cols,'zpeak')
    cols=np.append(cols,'zmean')

    array=array[(array['zpeak']>0.3)&(array['zpeak']<1.3)]

    if full:
      return np.copy(array)
    else:
      return [(array[col]) for col in cols]

  @staticmethod
  def get_cat_cols_ng(dir,cols,cuts,full,maxiter):
    noval=999999
    len0=0

    for ifile,file in enumerate(glob.glob(dir+'*.fits.gz')):
      print 'test'
      if ifile>maxiter:
        print 'test break found'
        break
      print 'file',ifile,file
      try:
        fits=fio.FITS(file)
      except IOError:
        print 'error loading fits file: ',file
        return

      colex,colist=CatalogMethods.col_exists(cols,fits[-1].get_colnames())
      if colex<1:
        raise ColError('columns '+colist+' do not exist in file: '+file)
        print 'there are column name(s) not in file: ',file
        return

      colex,colist=CatalogMethods.col_exists(cuts['col'],fits[-1].get_colnames())
      if colex<1:
        print 'there are cut column name(s) not in file: ',file
        raise ColError('cut columns '+colist+' do not exist in file: '+file)
        return

      tmparray=fits[-1].read()
      len0+=len(tmparray)

      mask=np.array([])
      for icut,cut in enumerate(cuts): 
        mask=CatalogMethods.cuts_on_col(mask,tmparray,cut['col'],cut['min'],cut['eq'],cut['max'])

      tmparray=fits[-1].read(columns=cols)

      if ifile==0:
        array=tmparray[mask]
      else:
        array=np.append(array,tmparray[mask],axis=0)

    print 'selected '+str(len(array))+' galaxies from catalog out of '+str(len0)

    if full:
      return np.copy(array)
    else:
      return [(array[col]) for col in cols]

  @staticmethod
  def col_exists(cols,colnames):

    colist=''
    exists=np.in1d(cols,colnames)
    for i,val in enumerate(exists):
      if val==0:
        colist+=' '+cols[i]


    return np.sum(exists)/len(cols),colist

  @staticmethod
  def cuts_on_col(mask,array,col,valmin,valeq,valmax):
    noval=999999

    if mask.size==0:
      mask=np.ones((len(array[col])), dtype=bool)

    if (valmin==noval) & (valmax==noval):
      if valeq==noval:
        print 'warning, no range or equality set in cut on column '+col
      else:
        mask=mask & (array[col]==valeq)
    elif (valmin!=noval) & (valmax!=noval):
      if valeq!=noval:
        print 'cannot have both equality and range cut on column '+col
      else:
        mask=mask & (valmin<array[col]) & (array[col]<valmax)
    elif (valmin!=noval):
      mask=mask & (valmin<array[col])
    else:
      mask=mask & (array[col]<valmax)
    return mask

  @staticmethod
  def col_view(array, cols):

    dtype2 = np.dtype({name:array.dtype.fields[name] for name in cols})
    return np.ndarray(array.shape, dtype2, array, 0, array.strides)

  @staticmethod
  def add_cut(cuts,col,min,eq,max):
    
    if cuts.size==0:
      cuts=np.zeros((1), dtype=[('col',np.str,20),('min',np.float64),('eq',np.float64),('max',np.float64)])
      cuts[0]['col']=col
      cuts[0]['min']=min
      cuts[0]['eq']=eq
      cuts[0]['max']=max
    else:
      cuts0=np.zeros((1), dtype=[('col',np.str,20),('min',np.float64),('eq',np.float64),('max',np.float64)])
      cuts0[0]['col']=col
      cuts0[0]['min']=min
      cuts0[0]['eq']=eq
      cuts0[0]['max']=max
      cuts=np.append(cuts,cuts0,axis=0)

    return cuts

  @staticmethod
  def match_cat_mask(cat1,cat2):

    mask1=np.in1d(cat1,cat2,assume_unique=True)
    mask2=np.in1d(cat2,cat1,assume_unique=True)

    return mask1,mask2

  @staticmethod
  def match_cat(cat1,cat2,both):

    if both:
      mask1,mask2=match_cat_mask(cat1,cat2)
      return cat1[mask1],cat2[mask2]
    else:
      mask1=match_cat_mask(cat1,cat2,both)
      return cat1[mask1]

    return

  @staticmethod
  def cat_name(cat):

    if cat==0:
      return 'im3shapev7'
    else:
      return 'ngmix010'

    return

  @staticmethod
  def cat_num(cat):

    if cat=='im3shapev7':
      return 0
    elif cat=='ngmix010':
      return 1

    return

  @staticmethod
  def get_pz(coadd):

    dir='/share/des/sv/photoz/tpz/'

    for ifile,file in enumerate(glob.glob(dir+'*.fits')):
      print 'file',ifile,file
      try:
        fits=fio.FITS(file)
      except IOError:
        print 'error loading fits file: ',file
        #return
      if ifile==0:
        zcoadd=fits[-1].read(columns=['coadd_objects_id'])['coadd_objects_id']
        zpeak=fits[-1].read(columns=['z_peak'])['z_peak']
        zmean=fits[-1].read(columns=['z_mean'])['z_mean']
      else:
        zcoadd=np.append(zcoadd,fits[-1].read(columns=['coadd_objects_id'])['coadd_objects_id'],axis=0)
        zmean=np.append(zmean,fits[-1].read(columns=['z_mean'])['z_mean'],axis=0)
        zpeak=np.append(zpeak,fits[-1].read(columns=['z_peak'])['z_peak'],axis=0)

    mask1,mask2=CatalogMethods.match_cat_mask(coadd,zcoadd)

    zcoadd=zcoadd[mask2]
    zmean=zmean[mask2]
    zpeak=zpeak[mask2]

    sort=np.argsort(zcoadd)[np.argsort(np.argsort(coadd[mask1]))]

    return zpeak[sort],zmean[sort],mask1


  @staticmethod
  def get_new_nbcw(coadd):

    file='/share/des/sv/im3shape/v7.2/id_nbc_w.fits'
    fits=fio.FITS(file)
    coadd2=fits[-1].read(columns=['coadd_id'])['coadd_id']
    c12=fits[-1].read(columns=['nbc_c1'])['nbc_c1']
    c22=fits[-1].read(columns=['nbc_c2'])['nbc_c2']
    m2=fits[-1].read(columns=['nbc_m'])['nbc_m']
    w2=fits[-1].read(columns=['w'])['w']

    mask1,mask2=CatalogMethods.match_cat_mask(coadd2,coadd)

    coadd2=coadd2[mask1]
    c12=c12[mask1]
    c22=c22[mask1]
    m2=m2[mask1]
    w2=w2[mask1]

    sort=np.argsort(coadd2)[np.argsort(np.argsort(coadd[mask2]))]

    coadd2=coadd2[sort]
    c12=c12[sort]
    c22=c22[sort]
    m2=m2[sort]
    w2=w2[sort]

    return c12,c22,m2,w2


