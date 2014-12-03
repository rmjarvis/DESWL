import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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
import sys
import ctypes
import numpy.linalg as linalg
  
class CatalogStore(object):

  coadd=[]
  ra=[]
  dec=[]
  e1=[]
  e2=[]
  m=[]
  c1=[]
  c2=[]
  s1=[]
  s2=[]
  w=[]
  psf1=[]
  psf2=[]
  psffwhm=[]
  radius=[]
  snr=[]
  regs=[]
  cat=0
  tbins=7
  sbins=2
  lbins=5
  sep=np.array([1,120])
  slop=0.1
  use_jk=False
  bs=True
  wt=True
  num_reg=152
  zpeak=[]
  ebv=[]
  exptime=[]
  maglimit=[]
  skysigma=[]
  skybrite=[]
  airmass=[]
  fwhm=[]
  colour=[]

  def __init__(self,cat,setup):

    noval=999999

    if setup:
      if CatalogMethods.cat_name(cat)=='im3shapev7':
        dir='/share/des/sv/im3shape/v7/main_cats/'
        cuts=CatalogMethods.default_im3shapev7_cuts()
        cols=np.array(['coadd_objects_id','ra','dec','e1','e2','mean_psf_e1_sky','mean_psf_e2_sky','mean_psf_fwhm','nbc_m','nbc_c1','nbc_c2','weight','snr','radius'])
        coadd,ra,dec,e1,e2,psf1,psf2,fwhm,m,c2,c1,w,snr,radius,zpeak=CatalogMethods.get_cat_cols(dir,cols,cuts,False,noval)
        s1=np.ones((len(ra)))
        s2=np.ones((len(ra)))
        c12,c22,m2,w2=CatalogMethods.get_new_nbcw(coadd)
        c11=c1
        c21=c2
        m1=m
        w1=w
        c1=c11
        c2=c21
        m=m1
        w=w1
        c1=c12
        c2=c22
        m=m2
        w=w2
      elif CatalogMethods.cat_name(cat)=='ngmix010':
        dir='/share/des/sv/ngmix/v009/'
        cuts=CatalogMethods.default_ngmix009_cuts()
        cols=np.array(['coadd_objects_id','ra','dec','exp_e_1','exp_e_2','psfrec_e_1','psfrec_e_2','psfrec_t','exp_e_sens_1','exp_e_sens_2','exp_e_cov_1_1','exp_e_cov_2_2','exp_e_cov_1_2','exp_t','exp_s2n_w'])
        coadd,ra,dec,e1,e2,psf1,psf2,fwhm,s1,s2,cov11,cov22,cov12,radius,snr,zpeak=CatalogMethods.get_cat_cols(dir,cols,cuts,False,noval)
        e1=-e1
        psf1=-psf1
        cov12=-cov12
        fwhm=np.sqrt(fwhm/2.)
        radius=np.sqrt(radius/2.)
        w=lin_shear_test_methods.ngmix_weight_calc(cov11,cov22,cov12)
        m=np.ones((len(ra))) #np.sqrt(s1*s2)
        c1=np.zeros((len(ra)))
        c2=np.zeros((len(ra)))

      num_reg,regs=jackknife_methods.load_jk_regs(ra,dec,64)


      lbins=5
      sbins=2
      tbins=7
      slop=0.1
      sep=np.array([1,120])
      use_jk=True
      bs=True
      wt=True

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
      self.psf1=psf1
      self.psf2=psf2
      self.psffwhm=fwhm
      self.radius=radius
      self.snr=snr
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
      self.zpeak=zpeak

      ebv,exptime,maglimit,skysigma,skybrite,airmass,fwhm=split_systematics.init_systematics_maps(ra,dec)

      if CatalogMethods.cat_name(cat)=='im3shapev7':
        self.ebv=ebv
        self.exptime=exptime[1,:]
        self.maglimit=maglimit[1,:]
        self.skysigma=skysigma[1,:]
        self.skybrite=skybrite[1,:]
        self.airmass=airmass[1,:]
        self.fwhm=fwhm[1,:]
      if CatalogMethods.cat_name(cat)=='ngmix010':
        self.ebv=ebv
        self.exptime=np.mean(exptime,axis=0)
        self.maglimit=np.mean(maglimit,axis=0)
        self.skysigma=np.mean(skysigma,axis=0)
        self.skybrite=np.mean(skybrite,axis=0)
        self.airmass=np.mean(airmass,axis=0)
        self.fwhm=np.mean(fwhm,axis=0)

      match=np.genfromtxt('../catcutDES/cat_match_dec.dat', names=True)

      mask1=np.in1d(coadd,match['coadd_objects_id'],assume_unique=True)
      mask2=np.in1d(match['coadd_objects_id'],coadd,assume_unique=True)
      match2=match[mask2]
      sort=np.argsort(match2['coadd_objects_id'])[np.argsort(np.argsort(coadd[mask1]))]
      match2=match2[sort]
      self.colour=match2['mag_auto_r']-match2['mag_auto_g']

class lin_shear_test_methods(object):

  @staticmethod
  def get_lin_e_w_norm(cat,xi):

    if CatalogMethods.cat_name(cat.cat)=='im3shapev7':

      if cat.bs & cat.wt:
        mw=cat.w
        me1=(cat.e1-cat.c1)
        me2=(cat.e2-cat.c2)
        norm=cat.w*(cat.m+1)
      elif cat.bs and not cat.wt:
        mw=np.ones((len(cat.e1)))
        me1=(cat.e1-cat.c1)
        me2=(cat.e2-cat.c2)
        norm=cat.m+1
      elif cat.wt and not cat.bs:
        mw=cat.w
        me1=cat.e1
        me2=cat.e2
        norm=cat.w
      else:
        mw=np.ones((len(cat.e1)))
        me1=cat.e1
        me2=cat.e2
        norm=np.ones((len(cat.e1)))

    elif CatalogMethods.cat_name(cat.cat)=='ngmix010':

      if cat.bs & cat.wt:
        mw=cat.w
        me1=cat.e1
        me2=cat.e2
        if xi:
          norm=cat.w*np.sqrt(cat.s1*cat.s2)
        else:
          norm=cat.w*(cat.s1+cat.s2)/2.
      elif cat.bs and not cat.wt:
        mw=np.ones((len(cat.e1)))
        me1=cat.e1
        me2=cat.e2
        if xi:
          norm=cat.w*np.sqrt(cat.s1*cat.s2)
        else:
          norm=cat.w*(cat.s1+cat.s2)/2.
      elif cat.wt and not cat.bs:
        mw=cat.w
        me1=cat.e1
        me2=cat.e2
        norm=cat.w
      else:
        mw=np.ones((len(cat.e1)))
        me1=cat.e1
        me2=cat.e2
        norm=np.ones((len(cat.e1)))

    return me1,me2,mw,norm

  @staticmethod
  def calc_mean_std_rms_e(cat,mask):

    me1,me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(cat,False)
    mean1=np.sum(mw[mask]*me1[mask])/np.sum(norm[mask])
    mean2=np.sum(mw[mask]*me2[mask])/np.sum(norm[mask])
    std1=(np.sum(mw[mask]*(me1[mask]-mean1)**2)/np.sum(norm[mask]))
    std2=(np.sum(mw[mask]*(me2[mask]-mean2)**2)/np.sum(norm[mask]))
    rms1=np.sqrt(np.sum((mw[mask]*me1[mask])**2)/np.sum(mw[mask]**2))
    rms2=np.sqrt(np.sum((mw[mask]*me2[mask])**2)/np.sum(mw[mask]**2))

    return mean1,mean2,std1,std2,rms1,rms2

  @staticmethod    
  def find_bin_edges(x,nbins):

    xs=np.sort(x)
    m=len(xs)
    r=np.linspace(0.,1.,nbins+1.)*(m-1)

    return xs[r.astype(int)]

  @staticmethod
  def binned_means_e(bin,cat,jkmask):

    y_mean1=[]
    y_std1=[]
    y_mean2=[]
    y_std2=[]

    for i in xrange(cat.lbins):
      mask=(bin==i)&jkmask
      mean1,mean2,std1,std2,rms1,rms2=lin_shear_test_methods.calc_mean_std_rms_e(cat,mask)
      y_mean1.append(mean1)
      y_std1.append(std1/np.sqrt(len(cat.e1[mask])))
      y_mean2.append(mean2)
      y_std2.append(std2/np.sqrt(len(cat.e2[mask])))

    y_mean1=np.array(y_mean1)
    y_std1=np.array(y_std1)
    y_mean2=np.array(y_mean2)
    y_std2=np.array(y_std2)

    return y_mean1,y_std1,y_mean2,y_std2

  @staticmethod
  def bin_means(x,cat,jkmask):

    edge=lin_shear_test_methods.find_bin_edges(x[jkmask],cat.lbins)
    xbin=np.digitize(x,edge)-1

    e1_mean,e1_std,e2_mean,e2_std=lin_shear_test_methods.binned_means_e(xbin,cat,jkmask)

    x_mean=[]
    x_std=[]
    for i in xrange(cat.lbins):
      mask=(xbin==i)&jkmask
      x_mean.append(np.mean(x[mask]))
      x_std.append(np.std(x[mask]))
    x_mean=np.array(x_mean)
    x_std=np.array(x_std)

    return x_mean,x_std,e1_mean,e1_std,e2_mean,e2_std

  @staticmethod
  def ngmix_weight_calc(cov11,cov22,cov12):
    sn=0.16

    w=1./(2.*sn**2.+cov11+cov22+2.*cov12)
    print w[w<0]

    return w

class xi_2pt_shear_test_methods(object):

  @staticmethod
  def xi_2pt(cat,mask,w2):

    me1,me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(cat,True)

    gg = treecorr.GGCorrelation(nbins=cat.tbins, min_sep=cat.sep[0], max_sep=cat.sep[1], sep_units='arcmin',binslop=cat.slop,verbose=0)
    kk = treecorr.KKCorrelation(nbins=cat.tbins, min_sep=cat.sep[0], max_sep=cat.sep[1], sep_units='arcmin',binslop=cat.slop,verbose=0)

    cate=treecorr.Catalog(g1=mw[mask]*me1[mask]*w2, g2=mw[mask]*me2[mask]*w2, ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')
    catm=treecorr.Catalog(k=norm[mask]*w2, ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')

    kk.process(catm,catm)
    kkp = np.copy(kk.xi)

    gg.process(cate)
    ggp = np.copy(gg.xip)/kkp
    ggm = np.copy(gg.xim)/kkp
    ggperr = ggmerr = np.sqrt(np.copy(gg.varxi))/kkp

    theta = np.exp(gg.meanlogr)
    print 'xip',ggp

    if cat.use_jk:
      ggperr,ggmerr=jackknife_methods.xi_2pt_err0(cat,mask,w2)

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

    exptime=np.zeros((5,len(ra)))
    maglimit=np.zeros((5,len(ra)))
    skysigma=np.zeros((5,len(ra)))
    skybrite=np.zeros((5,len(ra)))
    airmass=np.zeros((5,len(ra)))
    fwhm=np.zeros((5,len(ra)))

    ebv=split_systematics.load_sys_map_to_array(ra,dec,dir+'Planck_EBV_2048r_Q.fits',False,2048,True)
    exptime[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    maglimit[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    skysigma[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    airmass[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[4,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_Y_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)

    return ebv,exptime,maglimit,skysigma,skybrite,airmass,fwhm


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
  def split_gals_lin_along_(array,cat,label):

    arr1,tmp,me1,e1err,me2,e2err=lin_shear_test_methods.bin_means(array,cat,np.ones((len(cat.e1))).astype(bool))

    slp1,b1=np.polyfit(arr1, me1, 1)
    slp2,b2=np.polyfit(arr1, me2, 1)

    if cat.use_jk:
      e1err,e2err=jackknife_methods.lin_err0(array,cat)

    #plotting_methods.fig_create_e12vs(me1,me2,e1err,e2err,arr1,cat,label)
    print 'slope of e1 vs '+label,slp1
    print 'slope of e2 vs '+label,slp2

    return arr1,me1,me2,e1err,e2err,slp1,slp2,b1,b2

  @staticmethod
  def split_gals_2pt_along_(array,cat,label):

    jkst=cat.use_jk
    cat.use_jk=False

    xip=np.zeros((cat.sbins+1,cat.tbins))
    xim=np.zeros((cat.sbins+1,cat.tbins))
    var=np.zeros((cat.sbins+1,cat.tbins))
    xiperr=np.zeros((cat.sbins+1,cat.tbins))
    ximerr=np.zeros((cat.sbins+1,cat.tbins))
    edge=lin_shear_test_methods.find_bin_edges(array,cat.sbins)
    print 'bin edges',edge
    bin=np.digitize(array,edge)-1
    for i in xrange(cat.sbins):
      mask=(bin==i)
      print 'mean splits',np.mean(array[mask])
      list_split=[(cat.zpeak,cat.w),(cat.zpeak[mask],cat.w[mask])]
      list_nz_weights = hnz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3)
      w2=list_nz_weights[1]
      if i==0:
        w20=list_nz_weights[1]
        mask0=mask
      if i==(cat.sbins-1):
        w21=list_nz_weights[1]
        mask1=mask
      theta,xip[i,:],xim[i,:],xiperr[i,:],tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask,w2)
      #if cat.use_jk:
      #  xiperr[i,:],tmp,invcovp,tmp=jackknife_methods.xi_2pt_err0(cat,mask,w2)
      plotting_methods.fig_open_xi_2pt0(theta*(1+.02*(i+1)),xip[i,:],xiperr[i,:],cat,'bin '+str(i+1))

    mask=np.ones((len(cat.coadd))).astype(bool)
    theta,xip[cat.sbins,:],tmp,xiperr[cat.sbins,:],tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask,np.ones((len(cat.coadd))))
    dxip1=xip[0,:]-xip[cat.sbins,:]
    dxip2=xip[cat.sbins-1,:]-xip[cat.sbins,:]
    cat.use_jk=jkst
    if cat.use_jk:
      #xiperr[cat.sbins,:],tmp,invcovp,tmp=jackknife_methods.xi_2pt_err0(cat,mask,np.ones((len(cat.coadd))))
      xiperr1,xiperr2,chi2=jackknife_methods.xi_2pt_split_err0(cat,mask0,w20,mask1,w21)
    else:
      xiperr1=xiperr[0,:]
      xiperr2=xiperr[cat.sbins-1,:]
      chi2=0.
    plotting_methods.fig_open_xi_2pt0(theta,xip[cat.sbins,:],xiperr[cat.sbins,:],cat,'full')
    plotting_methods.fig_close_xi_2pt0(cat,label)

    return theta,edge,xip[cat.sbins,:],xiperr[cat.sbins,:],dxip1,dxip2,xiperr1,xiperr2,chi2


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
  def jk_chi2(invcov,xi1,xi2,cat):

    chi2=0.
    for i in xrange(cat.tbins):
      for j in xrange(cat.tbins):
        chi2+=((xi1[i]-xi2[i])*(xi1[j]-xi2[j])*invcov[i,j])

    return chi2

  @staticmethod
  def jk_corr_invcov_lin(cov,cat):

    invcov=(cat.num_reg-cat.lbins-2)/(cat.num_reg-1)*linalg.inv(cov)

    return invcov

  @staticmethod
  def jk_corr_invcov(cov,cat):

    invcov=(cat.num_reg-cat.tbins-2)/(cat.num_reg-1)*linalg.inv(cov)

    return invcov

  @staticmethod
  def jk_cov_lin(array,cat):

    cov=np.zeros((cat.lbins,cat.lbins))

    for i in xrange(cat.lbins):
      for j in xrange(cat.lbins):
        cov[i,j]=np.sum((array[:,i]-np.mean(array[:,i]))*(array[:,j]-np.mean(array[:,j])))/(cat.num_reg-1.)

    return cov

  @staticmethod
  def jk_cov(array,cat):

    cov=np.zeros((cat.tbins,cat.tbins))

    for i in xrange(cat.tbins):
      for j in xrange(cat.tbins):
        cov[i,j]=np.sum((array[:,i]-np.mean(array[:,i]))*(array[:,j]-np.mean(array[:,j])))/(cat.num_reg-1.)

    return cov

  @staticmethod
  def lin_err0(array,cat):

    e1=np.zeros((cat.num_reg,cat.lbins))
    e2=np.zeros((cat.num_reg,cat.lbins))

    for i in xrange(cat.num_reg):
      print 'jk',i
      mask=(cat.regs!=i)
      arr1,tmp,e1[i,:],tmp,e2[i,:],tmp=lin_shear_test_methods.bin_means(array,cat,mask)

    err1=np.sqrt(np.diagonal(jackknife_methods.jk_cov_lin(e1,cat)))
    err2=np.sqrt(np.diagonal(jackknife_methods.jk_cov_lin(e2,cat)))

    return err1,err2

  @staticmethod
  def xi_2pt_err0(cat,mask,w2):

    xip=np.zeros((cat.num_reg,cat.tbins))
    xim=np.zeros((cat.num_reg,cat.tbins))
    tmptime=time.time()

    jkst=cat.use_jk
    cat.use_jk=False
    for i in xrange(cat.num_reg):
      print 'jack_iter', i, time.time()-tmptime
      jkmask=(cat.regs!=i)
      tmp,xip[i,:],xim[i,:],tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask,w2)
      print 'jack xip',xip[i,:]
    cat.use_jk=jkst

    covp=jackknife_methods.jk_cov(xip,cat)
    covm=jackknife_methods.jk_cov(xim,cat)

    #invcovp=jackknife_methods.jk_corr_invcov(covp,cat)
    #invcovm=jackknife_methods.jk_corr_invcov(covm,cat)

    xiperr=np.sqrt(np.diagonal(covp))
    ximerr=np.sqrt(np.diagonal(covm))

    return xiperr,ximerr,invcovp,invcovm


  @staticmethod
  def xi_2pt_split_err0(cat,mask0,w20,mask1,w21):

    dxip=np.zeros((cat.num_reg,cat.tbins))
    dxip1=np.zeros((cat.num_reg,cat.tbins))
    dxip2=np.zeros((cat.num_reg,cat.tbins))
    tmptime=time.time()

    jkst=cat.use_jk
    cat.use_jk=False
    for i in xrange(cat.num_reg):
      print 'jack_iter_split', i, time.time()-tmptime
      jkmask=(cat.regs!=i)
      tmp,xip1,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask0,w20[jkmask[mask0]])
      tmp,xip2,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask1,w21[jkmask[mask1]])
      tmp,xip,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask,np.ones((len(cat.coadd)))[jkmask])
      dxip[i,:]=xip2-xip1
      dxip1[i,:]=xip1-xip
      dxip2[i,:]=xip2-xip
      print 'jack dxip',dxip[i,:]
    tmp,xip01,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask0,w20)
    tmp,xip02,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask1,w21)
    cat.use_jk=jkst

    covp=jackknife_methods.jk_cov(dxip,cat)
    covp1=jackknife_methods.jk_cov(dxip1,cat)
    covp2=jackknife_methods.jk_cov(dxip2,cat)

    if linalg.cond(covp) < 1/sys.float_info.epsilon:
      invcovp=jackknife_methods.jk_corr_invcov(covp,cat)
      chi2=jackknife_methods.jk_chi2(invcovp,xip02,xip01,cat)
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2=jackknife_methods.jk_chi2(1./covp,xip02,xip01,cat)

    if linalg.cond(covp1) < 1/sys.float_info.epsilon:
      invcovp1=jackknife_methods.jk_corr_invcov(covp1,cat)
      xiperr1=1./np.sqrt(np.diagonal(invcovp1))
    else:
      print 'xiperr1 val likely wrong, using approx due to singular covariance matrix'
      xiperr1=np.sqrt(np.diagonal(covp1))

    if linalg.cond(covp2) < 1/sys.float_info.epsilon:
      invcovp2=jackknife_methods.jk_corr_invcov(covp2,cat)
      xiperr2=1./np.sqrt(np.diagonal(invcovp2))
    else:
      print 'xiperr1 val likely wrong, using approx due to singular covariance matrix'
      xiperr2=np.sqrt(np.diagonal(covp2))


    return xiperr1,xiperr2,chi2


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
  def fig_create_e12vs(e1,e2,e1err,e2err,x,cat,label):

    plt.figure(10)
    plt.errorbar(x,e1, yerr=e1err, xerr=None,label='<e1>')
    plt.errorbar(x,e2, yerr=e2err, xerr=None,label='<e2>')
    plt.legend()
    plt.ylabel('<e>')
    plt.xlabel(label)
    plt.ylim((-5e-3,5e-3))
    plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_eVS'+label+'.png', bbox_inches='tight')
    plt.close(10)

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
  def fig_open_xi_2pt0(theta,xi,err,cat,label):

    plt.figure(10)
    plt.errorbar(theta,xi, yerr=err, xerr=None,marker='o',linestyle='',label=label)

    return

  @staticmethod
  def fig_close_xi_2pt0(cat,label):

    plt.figure(10)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\xi_+$')
    plt.xlabel(r'$\theta (arcmin)$')
    plt.xlim((1,200))
    plt.ylim((5e-7,1e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.tbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_xi_'+label+'.png', bbox_inches='tight')
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
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)

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
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)

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
    zpeak0['zpeak']=zpeak

    len1=len(array)
    array=array[zmask]
    array=nlr.merge_arrays([array,zpeak0],flatten=True)

    print 'selected '+str(len(array))+' galaxies from catalog out of '+str(len0)+' ('+str(len1)+' before pz cut)'

    cols=np.append(cols,'zpeak')

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

    zpeak,zmean,zmask=CatalogMethods.get_pz(array['coadd_objects_id'])
    zpeak0=np.zeros((len(zpeak),),dtype=[('zpeak','>f8')])
    zpeak0['zpeak']=zpeak

    len1=len(array)
    array=array[zmask]
    array=nlr.merge_arrays([array,zpeak0],flatten=True)

    print 'selected '+str(len(array))+' galaxies from catalog out of '+str(len0)+' ('+str(len1)+' before pz cut)'

    cols=np.append(cols,'zpeak')

    array=array[(array['zpeak']>0.3)&(array['zpeak']<1.3)]

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

class SVA1(object):

  @staticmethod
  def Fig_2ptPaper_amp(cat,xip,dxip,dxiperr):

    chi2st=999999
    amp=0.
    for a in xrange(-200,200):
      chi2=0.
      for i in xrange(cat.tbins):
        chi2=np.sum((dxip-a*xip/100.)**2./dxiperr**2.)
      if chi2<chi2st:
        chi2st=chi2
        amp=a/100.

    return amp

  @staticmethod
  def Fig_2ptPaper_subplot(cat,fig,r,c,n,array,label):

    plt.figure(fig)
    print ' '
    print ' '
    print '..........'+label+'..........'
    print ' '
    print ' '
    theta,edge,xip,xiperr,dxip1,dxip2,xiperr1,xiperr2,chi2=split_systematics.split_gals_2pt_along_(array,cat,label)
    edge=edge.astype(int)
    A1=SVA1.Fig_2ptPaper_amp(cat,xip,dxip1,xiperr1)
    A2=SVA1.Fig_2ptPaper_amp(cat,xip,dxip2,xiperr2)
    print 'theta',theta
    print 'xip',xip
    print 'xiperr',xiperr
    print 'dxip1',dxip1
    print 'dxiperr1',xiperr1
    print 'dxip2',dxip2
    print 'dxiperr2',xiperr2
    print 'chi',chi2
    print 'A',A1,A2
    ax=plt.subplot(r,c,n)
    plt.errorbar(theta[xip>0],xip[xip>0]*theta[xip>0],yerr=xiperr[xip>0]*theta[xip>0],marker='.',linestyle='',color='k')
    plt.errorbar(theta[xip<0],xip[xip<0]*theta[xip<0],yerr=xiperr[xip<0]*theta[xip<0],marker='.',linestyle='',color='k')
    plt.errorbar(theta,np.abs(xip*A1*theta),marker='',linestyle=':',color='r')
    plt.errorbar(theta[dxip1>0]*(1.2),dxip1[dxip1>0]*theta[dxip1>0],yerr=xiperr1[dxip1>0]*theta[dxip1>0],marker='v',linestyle='',color='r',label='snr=['+str(edge[0])+','+str(edge[1])+']')
    plt.errorbar(theta[dxip1<0]*(1.2),-dxip1[dxip1<0]*theta[dxip1<0],yerr=xiperr1[dxip1<0]*theta[dxip1<0],marker='1',linestyle='',color='r')
    plt.errorbar(theta,np.abs(xip*A2*theta),marker='',linestyle=':',color='b')
    plt.errorbar(theta[dxip2>0]*(1.4),dxip2[dxip2>0]*theta[dxip2>0],yerr=xiperr2[dxip2>0]*theta[dxip2>0],marker='^',linestyle='',color='b',label='snr=['+str(edge[1])+','+str(edge[2])+']')
    plt.errorbar(theta[dxip2<0]*(1.4),-dxip2[dxip2<0]*theta[dxip2<0],yerr=xiperr2[dxip2<0]*theta[dxip2<0],marker='2',linestyle='',color='b')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,200))
    plt.ylim((1e-6,9e-4))
    #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
    if n<5:
      ax.set_xticklabels([])
    else:
      plt.xlabel(r'$\theta$')
    if n%2==0:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\theta\xi_+$')

    return

  @staticmethod
  def Fig_ShearPaper_subplot(cat,fig,r,c,n,array,label):

    plt.figure(fig)
    print ' '
    print ' '
    print '..........'+label+'..........'
    print ' '
    print ' '
    arr1,me1,me2,e1err,e2err,slp1,slp2,b1,b2=split_systematics.split_gals_lin_along_(array,cat,label)
    print 'psf',arr1
    print 'me',me1,me2
    print 'err',e1err,e2err
    print 'slp',slp1,slp2
    print 'b',b1,b2
    ax=plt.subplot(r,c,n)
    if n==3:
      plt.errorbar(arr1,me1,yerr=e1err,marker='.',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
      plt.errorbar(arr1,me2,yerr=e2err,marker='.',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
    else:
      plt.errorbar(arr1,me1,yerr=e1err,marker='.',linestyle='',color='r')
      plt.errorbar(arr1,me2,yerr=e2err,marker='.',linestyle='',color='b')
    plt.errorbar(arr1,slp1*arr1+b1,marker='',linestyle=':',color='r')
    plt.errorbar(arr1,slp2*arr1+b2,marker='',linestyle=':',color='b')
    plt.ylim((-.0025,.0025))
    plt.legend(loc='upper right',ncol=2, frameon=False,prop={'size':12})
    if n<4:
      ax.set_xticklabels([])
    else:
      if label=='psfe1':
        plt.xlabel(r'$PSF\quad e_1$')
      if label=='psfe2':
        plt.xlabel(r'$PSF\quad e_2$')
      if label=='psffwhm':
        plt.xlabel(r'$PSF\quad FWHM$')
    if (n==2)|(n==5):
      ax.set_yticklabels([])
      #plt.xlim((-.0075,.0225))
    elif (n==3)|(n==6):
      ax.set_yticklabels([])
      #plt.xlim((-.0075,.0225))
    else:
      plt.ylabel(r'$\langle e \rangle$')
      #plt.xlim((-.0025,.0175))

    return


  @staticmethod
  def Fig1_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 1'
    print '...........'

    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,1,i3.snr,'snr')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,2,ng.snr,'snr')
    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,3,i3.radius,'radius')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,4,ng.radius,'radius')
    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,5,i3.colour,'colour')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,6,ng.colour,'colour')
    plt.figure(1)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_Fig1.png', bbox_inches='tight')
    plt.close(1)   

    return

  @staticmethod
  def Fig2_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 2'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,1,i3.ra,'ra')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,2,ng.ra,'ra')
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,3,i3.dec,'dec')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,4,ng.dec,'dec')
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,5,i3.ebv,'ebv')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,6,ng.ebv,'ebv')
    plt.figure(2)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_Fig2.png', bbox_inches='tight')
    plt.close(2)   
    
    return

  @staticmethod
  def Fig3_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 3'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,1,i3.exptime,'exptime')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,2,ng.exptime,'exptime')
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,3,i3.maglimit,'maglimit')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,4,ng.maglimit,'maglimit')
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,5,i3.airmass,'airmass')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,6,ng.airmass,'airmass')
    plt.figure(3)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_Fig3.png', bbox_inches='tight')
    plt.close(3)   
    
    return


  @staticmethod
  def Fig4_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 4'
    print '...........'

  
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,1,i3.skysigma,'skysigma')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,2,ng.skysigma,'skysigma')
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,3,i3.skybrite,'skybrite')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,4,ng.skybrite,'skybrite')
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,5,i3.fwhm,'fwhm')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,6,ng.fwhm,'fwhm')
    plt.figure(4)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_Fig4.png', bbox_inches='tight')
    plt.close(4)   
    
    return

  @staticmethod
  def Fig5_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 5'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,1,i3.psf1,'psf1')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,2,ng.psf1,'psf1')
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,3,i3.psf2,'psf2')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,4,ng.psf2,'psf2')
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,5,i3.psffwhm,'psffwhm')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,6,ng.psffwhm,'psffwhm')
    plt.figure(5)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_Fig5.png', bbox_inches='tight')
    plt.close(5)   
    
    return

  @staticmethod
  def Fig1_ShearPaper(i3,ng):


    print '...........'
    print 'Figure 1'
    print '...........'

    SVA1.Fig_ShearPaper_subplot(i3,1,2,3,1,i3.psf1,'psfe1')    
    SVA1.Fig_ShearPaper_subplot(i3,1,2,3,2,i3.psf2,'psfe2')    
    SVA1.Fig_ShearPaper_subplot(i3,1,2,3,3,i3.psffwhm,'psffwhm')    
    SVA1.Fig_ShearPaper_subplot(ng,1,2,3,4,ng.psf1,'psfe1')    
    SVA1.Fig_ShearPaper_subplot(ng,1,2,3,5,ng.psf2,'psfe2')    
    SVA1.Fig_ShearPaper_subplot(ng,1,2,3,6,ng.psffwhm,'psffwhm')    
    plt.figure(1)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.lbins)+'_Fig1.png', bbox_inches='tight')
    plt.close(1)   
    
    return


