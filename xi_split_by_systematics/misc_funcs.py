import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
#plt.style.use('SVA1StyleSheet.mplstyle')
import matplotlib.ticker
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
  sbins=2
  lbins=10
  tbins=9
  sep=np.array([1,400])
  slop=0.15
  num_reg=152
  jkxipst0=np.zeros((num_reg,tbins))
  meanst=np.zeros((sbins))
  edgest=np.zeros((sbins+1))
  astp=np.zeros((4))
  chi2stp=np.zeros((3))
  astm=np.zeros((4))
  chi2stm=np.zeros((3))
  use_jk=True
  bs=True
  wt=True
  use_zrw=True
  ztyp='tpzv1'
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
        dir='/share/des/sv/im3shape/v7.2b/'
        cuts=CatalogMethods.default_im3shapev7_cuts()
        cols=np.array(['coadd_objects_id','ra','dec','e1','e2','mean_psf_e1_sky','mean_psf_e2_sky','mean_psf_fwhm','nbc_m','nbc_c1','nbc_c2','w','snr','mean_rgpp_rp','z_peak'])
        coadd,ra,dec,e1,e2,psf1,psf2,fwhm,m,c1,c2,w,snr,radius,zpeak=CatalogMethods.get_cat_cols_i3(dir,cols,cuts,False,noval)
        s1=np.ones((len(ra)))
        s2=np.ones((len(ra)))
      elif CatalogMethods.cat_name(cat)=='ngmix010':
        dir='/share/des/sv/ngmix/v010/'
        cuts=CatalogMethods.default_ngmix009_cuts()
        cols=np.array(['coadd_objects_id','ra','dec','exp_e_1','exp_e_2','psfrec_e_1','psfrec_e_2','psfrec_t','exp_e_sens_1','exp_e_sens_2','exp_e_cov_1_1','exp_e_cov_2_2','exp_e_cov_1_2','exp_t','exp_s2n_w','z_peak'])
        coadd,ra,dec,e1,e2,psf1,psf2,fwhm,s1,s2,cov11,cov22,cov12,radius,snr,zpeak=CatalogMethods.get_cat_cols_i3(dir,cols,cuts,False,noval)
        e1=-e1
        psf1=-psf1
        cov12=-cov12
        fwhm=np.sqrt(fwhm/2.)
        radius=np.sqrt(radius/2.)
        w=lin_shear_test_methods.ngmix_weight_calc(cov11,cov22,cov12)
        m=np.ones((len(ra))) #np.sqrt(s1*s2)
        c1=np.zeros((len(ra)))
        c2=np.zeros((len(ra)))

      num_reg,regs=jackknife_methods.load_jk_regs(ra,dec,32)


      self.lbins=10
      self.sbins=2
      self.slop=.15
      self.tbins=9
      self.sep=np.array([1,400])
      self.use_jk=True
      self.bs=True
      self.wt=True
      self.use_zrw=True
      self.jkxipst0=np.zeros((self.num_reg,self.tbins))
      self.ztyp='tpzv1'
      self.meanst=np.zeros((self.sbins))
      self.edgest=np.zeros((self.sbins+1))
      self.astp=np.zeros((4))
      self.chi2stp=np.zeros((3))
      self.astm=np.zeros((4))
      self.chi2stm=np.zeros((3))
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

      match=np.genfromtxt('/share/des/sv/des_sva1_gold_v1_0_id_type_absmags_stellarmass.dat', names=['coadd_objects_id','colour','1','2','3','4','5','6'],usecols=['coadd_objects_id','colour'])

      mask1=np.in1d(coadd,match['coadd_objects_id'],assume_unique=True)
      mask2=np.in1d(match['coadd_objects_id'],coadd,assume_unique=True)
      match2=match[mask2]
      sort=np.argsort(match2['coadd_objects_id'])[np.argsort(np.argsort(coadd[mask1]))]
      match2=match2[sort]
      self.colour=match2['colour']

class lin_shear_test_methods(object):

  @staticmethod
  def get_lin_e_w_norm(cat,xi):

    if CatalogMethods.cat_name(cat.cat)=='im3shapev7':

      if cat.bs & cat.wt:
        mw=cat.w
        me1=(cat.e1-cat.c1)
        me2=(cat.e2-cat.c2)
        norm=(cat.m+1)
      elif cat.bs and not cat.wt:
        mw=np.ones((len(cat.e1)))
        me1=(cat.e1-cat.c1)
        me2=(cat.e2-cat.c2)
        norm=cat.m+1
      elif cat.wt and not cat.bs:
        mw=cat.w
        me1=cat.e1
        me2=cat.e2
        norm=np.ones((len(cat.e1)))
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
          norm=np.sqrt(cat.s1*cat.s2)
        else:
          norm=(cat.s1+cat.s2)/2.
      elif cat.bs and not cat.wt:
        mw=np.ones((len(cat.e1)))
        me1=cat.e1
        me2=cat.e2
        if xi:
          norm=np.sqrt(cat.s1*cat.s2)
        else:
          norm=(cat.s1+cat.s2)/2.
      elif cat.wt and not cat.bs:
        mw=cat.w
        me1=cat.e1
        me2=cat.e2
        norm=np.ones((len(cat.e1)))
      else:
        mw=np.ones((len(cat.e1)))
        me1=cat.e1
        me2=cat.e2
        norm=np.ones((len(cat.e1)))

    return me1,me2,mw,norm

  @staticmethod
  def calc_mean_std_rms_e(cat,mask):

    me1,me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(cat,False)
    mean1=np.sum(mw[mask]*me1[mask])/np.sum(mw[mask]*norm[mask])
    mean2=np.sum(mw[mask]*me2[mask])/np.sum(mw[mask]*norm[mask])
    std1=(np.sum(mw[mask]*(me1[mask]-mean1)**2)/np.sum(mw[mask]*norm[mask]))
    std2=(np.sum(mw[mask]*(me2[mask]-mean2)**2)/np.sum(mw[mask]*norm[mask]))
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

    cate=treecorr.Catalog(g1=me1[mask], g2=me2[mask], w=mw[mask]*w2, ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')
    catm=treecorr.Catalog(k=norm[mask], w=mw[mask]*w2, ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')

    kk.process(catm,catm)
    kkp = np.copy(kk.xi)

    gg.process(cate)
    ggp = np.copy(gg.xip)/kkp
    ggm = np.copy(gg.xim)/kkp
    ggperr = ggmerr = np.sqrt(np.copy(gg.varxi))/kkp

    chi2p=chi2m=0.

    theta = np.exp(gg.meanlogr)
    print 'xip',ggp
    print 'xim',ggm

    if cat.use_jk:
      ggperr,ggmerr,chi2p,chi2m=jackknife_methods.xi_2pt_err0(cat,mask,w2)

    return theta,ggp,ggm,ggperr,ggmerr,chi2p,chi2m

  @staticmethod
  def xi_2pt_psf(cat,mask):

    me1,me2,mw,norm=lin_shear_test_methods.get_lin_e_w_norm(cat,True)

    gg = treecorr.GGCorrelation(nbins=cat.tbins, min_sep=cat.sep[0], max_sep=cat.sep[1], sep_units='arcmin',binslop=cat.slop,verbose=0)
    kk = treecorr.KKCorrelation(nbins=cat.tbins, min_sep=cat.sep[0], max_sep=cat.sep[1], sep_units='arcmin',binslop=cat.slop,verbose=0)
    nk = treecorr.NKCorrelation(nbins=cat.tbins, min_sep=cat.sep[0], max_sep=cat.sep[1], sep_units='arcmin',binslop=cat.slop,verbose=0)

    cate=treecorr.Catalog(g1=me1[mask], g2=me2[mask], w=mw[mask], ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')
    catpsf=treecorr.Catalog(g1=cat.psf1[mask], g2=cat.psf2[mask], ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')
    catm=treecorr.Catalog(k=norm[mask], w=mw[mask], ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')
    catn=treecorr.Catalog(ra=cat.ra[mask], dec=cat.dec[mask], ra_units='deg', dec_units='deg')

    nk.process(catn,catm)
    nkp = np.copy(nk.xi)
    kk.process(catm,catm)
    kkp = np.copy(kk.xi)

    gg.process(cate)
    ggp = np.copy(gg.xip)/kkp
    ggm = np.copy(gg.xim)/kkp
    ggperr = ggmerr = np.sqrt(np.copy(gg.varxi))/kkp

    gg.process(cate,catpsf)
    gpp = np.copy(gg.xip)/nkp
    gpm = np.copy(gg.xim)/nkp
    gpperr = gpmerr = np.sqrt(np.copy(gg.varxi))/nkp

    gg.process(catpsf)
    ppp=np.copy(gg.xip)
    ppm=np.copy(gg.xim)
    ppperr=ppmerr=np.sqrt(np.copy(gg.varxi))

    alphap=xi_2pt_shear_test_methods.calc_alpha(gpp,ppp,cat)
    alphaperr=0.
    alpham=xi_2pt_shear_test_methods.calc_alpha(gpm,ppm,cat)
    alphamerr=0.
    print 'alphap = ',alphap[cat.tbins/2]
    print 'alpham = ',alpham[cat.tbins/2]

    leakagep=xi_2pt_shear_test_methods.calc_psf_leakage(alphap,ppp,cat)
    leakagep/=ggp
    leakageperr=0.
    print 'psf leakagep = ',leakagep[cat.tbins/2]
    leakagem=xi_2pt_shear_test_methods.calc_psf_leakage(alpham,ppm,cat)
    leakagem/=ggp
    leakagemerr=0.
    print 'psf leakagem = ',leakagem[cat.tbins/2]

    theta = np.exp(gg.meanlogr)

    if cat.use_jk:
      gpperr,gpmerr,ppperr,ppmerr,alphaerr,leakageerr=jackknife_methods.xi_2pt_psf_err(ra,dec,e1,e2,psf1,psf2,m,c1,c2,s1,s2,w,cat,bins,sep,slop,mean_e1,mean_e2,mean_psf1,mean_psf2,ggp,bs,wt,regs,reg_num)

    return theta,ggp,ggm,ggperr,ggmerr,gpp,gpm,gpperr,gpmerr,ppp,ppm,ppperr,ppmerr,alphap,alphaperr,leakagep,leakageperr,alpham,alphamerr,leakagem,leakagemerr

  @staticmethod
  def calc_alpha(gp,pp,cat):

    e1,e2,tmp,tmp,tmp,tmp=lin_shear_test_methods.calc_mean_std_rms_e(cat,np.ones((len(cat.coadd))).astype(bool))
    psf1=np.mean(cat.psf1)
    psf2=np.mean(cat.psf2)
    alpha=(gp-e1*psf1-e2*psf2)/(pp-psf1**2-psf2**2)

    return alpha

  @staticmethod
  def calc_psf_leakage(alpha,pp,cat):

    e1,e2,tmp,tmp,tmp,tmp=lin_shear_test_methods.calc_mean_std_rms_e(cat,np.ones((len(cat.coadd))).astype(bool))
    psf1=np.mean(cat.psf1)
    psf2=np.mean(cat.psf2)
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

    exptime=np.zeros((4,len(ra)))
    maglimit=np.zeros((4,len(ra)))
    skysigma=np.zeros((4,len(ra)))
    skybrite=np.zeros((4,len(ra)))
    airmass=np.zeros((4,len(ra)))
    fwhm=np.zeros((4,len(ra)))

    ebv=split_systematics.load_sys_map_to_array(ra,dec,dir+'Planck_EBV_2048r_Q.fits',False,2048,True)
    exptime[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    exptime[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_EXPTIME__total.fits.gz',False,4096,False)
    maglimit[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    maglimit[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_maglimit__.fits.gz',False,4096,False)
    skysigma[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skysigma[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    skybrite[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',False,4096,False)
    airmass[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    airmass[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[0,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[1,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[2,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)
    fwhm[3,:]=split_systematics.load_sys_map_to_array(ra,dec,dir+'SVA1_IMAGE_SRC_band_z_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',False,4096,False)

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

    xip=np.zeros((cat.sbins+1,cat.tbins))
    xim=np.zeros((cat.sbins+1,cat.tbins))
    var=np.zeros((cat.sbins+1,cat.tbins))
    xiperr=np.zeros((cat.sbins+1,cat.tbins))
    ximerr=np.zeros((cat.sbins+1,cat.tbins))
    jkxipst=np.zeros((cat.sbins,cat.num_reg,cat.tbins))
    jkximst=np.zeros((cat.sbins,cat.num_reg,cat.tbins))

    if label=='Colour':
      bin=-np.ones((len(cat.coadd)))
      bin[(array==1)|(array==2)]=0
      bin[(array>3)]=1
      edge=np.array([0,3,0])
    else:
      edge=lin_shear_test_methods.find_bin_edges(array,cat.sbins)
      bin=np.digitize(array,edge)-1
    print 'bin edges',edge
    edgemean=np.zeros((cat.sbins))
    for i in xrange(cat.sbins):
      mask=(bin==i)
      edgemean[i]=np.mean(array[mask])
      print 'mean splits',np.mean(array[mask]),len(array[mask])
      if cat.use_zrw:
        list_split=[(cat.zpeak,cat.w),(cat.zpeak[mask],cat.w[mask])]
        list_nz_weights = hnz.get_weights(list_split,target_nz_index=0,photoz_min=0.3,photoz_max=1.3)
        w2=list_nz_weights[1]
        if i==0:
          w20=list_nz_weights[1]
        if i==(cat.sbins-1):
          w21=list_nz_weights[1]
      else:
        w2=np.ones((len(cat.coadd[mask])))
        if i==0:
          w20=np.ones((len(cat.coadd[mask])))
        if i==(cat.sbins-1):
          w21=np.ones((len(cat.coadd[mask])))
      if i==0:
        mask0=mask
      if i==(cat.sbins-1):
        mask1=mask

      theta,xip[i,:],xim[i,:],xiperr[i,:],ximerr[i,:],tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask,w2)
      plotting_methods.fig_open_xi_2pt0(theta*(1+.02*(i+1)),xip[i,:],xiperr[i,:],cat,'bin '+str(i+1))
      if cat.use_jk:
        jkxipst[i,:,:]=cat.jkxipst
        jkximst[i,:,:]=cat.jkximst

    jkst=cat.use_jk
    if cat.use_jk:
      if cat.jkxipst0[0,0]!=0.:
        cat.use_jk=False
    mask=np.ones((len(cat.coadd)))
    theta,xip[cat.sbins,:],xim[cat.sbins,:],xiperr[cat.sbins,:],ximerr[cat.sbins,:],tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask.astype(bool),mask)
    cat.use_jk=jkst
    dxip1=xip[0,:]-xip[cat.sbins,:]
    dxip2=xip[cat.sbins-1,:]-xip[cat.sbins,:]
    dxim1=xim[0,:]-xim[cat.sbins,:]
    dxim2=xim[cat.sbins-1,:]-xim[cat.sbins,:]
    plotting_methods.fig_open_xi_2pt0(theta,xip[cat.sbins,:],xiperr[cat.sbins,:],cat,'full')
    plotting_methods.fig_close_xi_2pt0(cat,label)

    if cat.use_jk:
      covp0=jackknife_methods.jk_cov(jkxipst[cat.sbins-1,:,:]-jkxipst[0,:,:],cat)
      covm0=jackknife_methods.jk_cov(jkximst[cat.sbins-1,:,:]-jkxipst[0,:,:],cat)
      covp1=jackknife_methods.jk_cov(jkxipst[0,:,:]-cat.jkxipst0,cat)
      covm1=jackknife_methods.jk_cov(jkximst[0,:,:]-cat.jkxipst0,cat)
      covp2=jackknife_methods.jk_cov(jkxipst[cat.sbins-1,:,:]-cat.jkxipst0,cat)
      covm2=jackknife_methods.jk_cov(jkximst[cat.sbins-1,:,:]-cat.jkxipst0,cat)

      print 'cov',covp0,covm0,covp1,covm1,covp2,covm2

      if cat.jkxipst0[0,0]!=0.:
        covp=jackknife_methods.jk_cov(cat.jkxipst0,cat)
        if linalg.cond(covp) < 1/sys.float_info.epsilon:
          invcovp=jackknife_methods.jk_corr_invcov(covp,cat)
          xiperr[cat.sbins,:]=1./np.sqrt(np.diagonal(invcovp))
        else:
          print 'chi2 val likely wrong, using approx due to singular covariance matrix'
          xiperr[cat.sbins,:]=np.sqrt(np.diagonal(covp))
        covm=jackknife_methods.jk_cov(cat.jkximst0,cat)
        if linalg.cond(covm) < 1/sys.float_info.epsilon:
          invcovm=jackknife_methods.jk_corr_invcov(covm,cat)
          ximerr[cat.sbins,:]=1./np.sqrt(np.diagonal(invcovm))
        else:
          print 'chi2 val likely wrong, using approx due to singular covariance matrix'
          ximerr[cat.sbins,:]=np.sqrt(np.diagonal(covm))

      if linalg.cond(covp0) < 1/sys.float_info.epsilon:
        invcovp0=jackknife_methods.jk_corr_invcov(covp0,cat)
        print 'invcovp0',invcovp0
        chi2p0=jackknife_methods.jk_chi2(invcovp0,dxip2,dxip1,cat)
        xiperr0=1./np.sqrt(np.diagonal(invcovp0))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2p0=jackknife_methods.jk_chi2(1./covp0,dxip2,dxip1,cat)
        xiperr0=np.sqrt(np.diagonal(covp0))
      if linalg.cond(covm0) < 1/sys.float_info.epsilon:
        invcovm0=jackknife_methods.jk_corr_invcov(covm0,cat)
        print 'invcovm0',invcovm0
        chi2m0=jackknife_methods.jk_chi2(invcovm0,dxim2,dxim1,cat)
        ximerr0=1./np.sqrt(np.diagonal(invcovm0))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2m0=jackknife_methods.jk_chi2(1./covm0,dxim2,dxim1,cat)
        ximerr0=np.sqrt(np.diagonal(covm0))

      if linalg.cond(covp1) < 1/sys.float_info.epsilon:
        invcovp1=jackknife_methods.jk_corr_invcov(covp1,cat)
        print 'invcovp1',invcovp1
        chi2p1=jackknife_methods.jk_chi2(invcovp1,dxip1,np.zeros((len(dxip1))),cat)
        xiperr1=1./np.sqrt(np.diagonal(invcovp1))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2p1=jackknife_methods.jk_chi2(1./covp1,dxip1,np.zeros((len(dxip1))),cat)
        xiperr1=np.sqrt(np.diagonal(covp1))
      if linalg.cond(covm1) < 1/sys.float_info.epsilon:
        invcovm1=jackknife_methods.jk_corr_invcov(covm1,cat)
        print 'invcovm1',invcovm1
        chi2m1=jackknife_methods.jk_chi2(invcovm1,dxim1,np.zeros((len(dxim1))),cat)
        ximerr1=1./np.sqrt(np.diagonal(invcovm1))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2m1=jackknife_methods.jk_chi2(1./covm1,dxim1,np.zeros((len(dxim1))),cat)
        ximerr1=np.sqrt(np.diagonal(covm1))

      if linalg.cond(covp2) < 1/sys.float_info.epsilon:
        invcovp2=jackknife_methods.jk_corr_invcov(covp2,cat)
        print 'invcovp2',invcovp2
        chi2p2=jackknife_methods.jk_chi2(invcovp2,dxip2,np.zeros((len(dxip2))),cat)
        xiperr2=1./np.sqrt(np.diagonal(invcovp2))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2p2=jackknife_methods.jk_chi2(1./covp2,dxip2,np.zeros((len(dxip2))),cat)
        xiperr2=np.sqrt(np.diagonal(covp2))
      if linalg.cond(covm2) < 1/sys.float_info.epsilon:
        invcovm2=jackknife_methods.jk_corr_invcov(covm2,cat)
        print 'invcovm2',invcovm2
        chi2m2=jackknife_methods.jk_chi2(invcovm2,dxim2,np.zeros((len(dxim2))),cat)
        ximerr2=1./np.sqrt(np.diagonal(invcovm2))
      else:
        print 'chi2 val likely wrong, using approx due to singular covariance matrix'
        chi2m2=jackknife_methods.jk_chi2(1./covm2,dxim2,np.zeros((len(dxim2))),cat)
        ximerr2=np.sqrt(np.diagonal(covm2))
    else:
      xiperr0=xiperr[cat.sbins,:]
      xiperr1=xiperr[0,:]
      xiperr2=xiperr[cat.sbins-1,:]
      chi2p0=0.
      chi2p1=0.
      chi2p2=0.
      ximerr0=ximerr[cat.sbins,:]
      ximerr1=ximerr[0,:]
      ximerr2=ximerr[cat.sbins-1,:]
      chi2m0=0.
      chi2m1=0.
      chi2m2=0.

    return theta,edge,edgemean,xip[cat.sbins,:],xiperr[cat.sbins,:],xim[cat.sbins,:],ximerr[cat.sbins,:],dxip2-dxip1,dxip1,dxip2,xiperr0,xiperr1,xiperr2,chi2p0,chi2p1,chi2p2,dxim2-dxim1,dxim1,dxim2,ximerr0,ximerr1,ximerr2,chi2m0,chi2m1,chi2m2

class jackknife_methods(object):


  @staticmethod
  def load_jk_regs(ra,dec,healno):

    regfile=np.genfromtxt('/share/des/sv/pix-list_im3shape-v7_ieflags_base-no_nside'+str(healno)+'-nest.txt')

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

    invcov=linalg.inv(cov)#*(cat.num_reg-cat.lbins-2.)/(cat.num_reg-1.)

    return invcov

  @staticmethod
  def jk_corr_invcov(cov,cat):

    invcov=linalg.inv(cov)#*(cat.num_reg-cat.tbins-2.)/(cat.num_reg-1.)

    return invcov

  @staticmethod
  def jk_cov_lin(array,cat):

    cov=np.zeros((cat.lbins,cat.lbins))

    for i in xrange(cat.lbins):
      for j in xrange(cat.lbins):
        cov[i,j]=np.sum((array[:,i]-np.mean(array[:,i]))*(array[:,j]-np.mean(array[:,j])))*(cat.num_reg-1.)/cat.num_reg

    return cov

  @staticmethod
  def jk_cov(array,cat):

    cov=np.zeros((cat.tbins,cat.tbins))

    for i in xrange(cat.tbins):
      for j in xrange(cat.tbins):
        cov[i,j]=np.sum((array[:,i]-np.mean(array[:,i]))*(array[:,j]-np.mean(array[:,j])))*(cat.num_reg-1.)/cat.num_reg


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
      tmp,xip[i,:],xim[i,:],tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask,w2[jkmask[mask]])
    tmp,xip0,xim0,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask,w2)
    cat.use_jk=jkst

    cat.jkxipst=xip
    cat.jkximst=xim
    if len(cat.coadd[mask])==len(cat.coadd):
      cat.jkxipst0=xip
      cat.jkximst0=xim


    covp=jackknife_methods.jk_cov(xip,cat)
    covm=jackknife_methods.jk_cov(xim,cat)

    print 'cov',covp,covm


    if linalg.cond(covp) < 1/sys.float_info.epsilon:
      invcovp=jackknife_methods.jk_corr_invcov(covp,cat)
      print 'invcovp',invcovp
      chi2p=jackknife_methods.jk_chi2(invcovp,xip0,np.zeros((len(xip0))),cat)
      xiperr=1./np.sqrt(np.diagonal(invcovp))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2p=jackknife_methods.jk_chi2(1./covp,xip0,np.zeros((len(xip0))),cat)
      xiperr=np.sqrt(np.diagonal(covp))
    if linalg.cond(covm) < 1/sys.float_info.epsilon:
      invcovm=jackknife_methods.jk_corr_invcov(covm,cat)
      print 'invcovm',invcovm
      chi2m=jackknife_methods.jk_chi2(invcovm,xim0,np.zeros((len(xip0))),cat)
      ximerr=1./np.sqrt(np.diagonal(invcovm))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2m=jackknife_methods.jk_chi2(1./covp,xim0,np.zeros((len(xip0))),cat)
      ximerr=np.sqrt(np.diagonal(covm))


    return xiperr,ximerr,chi2p,chi2m


  @staticmethod
  def xi_2pt_split_err0(cat,mask0,w20,mask1,w21):

    dxipp=np.zeros((cat.num_reg,cat.tbins))
    dxip1p=np.zeros((cat.num_reg,cat.tbins))
    dxip2p=np.zeros((cat.num_reg,cat.tbins))
    dxipm=np.zeros((cat.num_reg,cat.tbins))
    dxip1m=np.zeros((cat.num_reg,cat.tbins))
    dxip2m=np.zeros((cat.num_reg,cat.tbins))
    tmptime=time.time()

    jkst=cat.use_jk
    cat.use_jk=False
    for i in xrange(cat.num_reg):
      print 'jack_iter_split', i, time.time()-tmptime
      jkmask=(cat.regs!=i)
      tmp,xip1,xim1,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask0,w20[jkmask[mask0]])
      tmp,xip2,xim2,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask&mask1,w21[jkmask[mask1]])
      if cat.xi_jk_st==0:
        tmp,xip,xim,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,jkmask,np.ones((len(cat.coadd)))[jkmask])
      dxipp[i,:]=xip2-xip1
      dxip1p[i,:]=xip1-xip
      dxip2p[i,:]=xip2-xip
      dxipm[i,:]=xim2-xim1
      dxip1m[i,:]=xim1-xim
      dxip2m[i,:]=xim2-xim
      print 'jack dxipp',dxipp[i,:],dxip1p[i,:],dxip2p[i,:]
      print 'jack dxipm',dxipm[i,:],dxip1m[i,:],dxip2m[i,:]
    tmp,xip01,xim01,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask0,w20)
    tmp,xip02,xim02,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask1,w21)
    tmp,xip0,xim0,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,np.ones((len(cat.coadd))).astype(bool),np.ones((len(cat.coadd))))
    cat.use_jk=jkst

    covp=jackknife_methods.jk_cov(dxipp,cat)
    covp1=jackknife_methods.jk_cov(dxip1p,cat)
    covp2=jackknife_methods.jk_cov(dxip2p,cat)
    covm=jackknife_methods.jk_cov(dxipm,cat)
    covm1=jackknife_methods.jk_cov(dxip1m,cat)
    covm2=jackknife_methods.jk_cov(dxip2m,cat)
    if linalg.cond(covp) < 1/sys.float_info.epsilon:
      invcovp=jackknife_methods.jk_corr_invcov(covp,cat)
      chi2p=jackknife_methods.jk_chi2(invcovp,xip02,xip01,cat)
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2p=jackknife_methods.jk_chi2(1./covp,xip02,xip01,cat)
    if linalg.cond(covm) < 1/sys.float_info.epsilon:
      invcovm=jackknife_methods.jk_corr_invcov(covm,cat)
      chi2m=jackknife_methods.jk_chi2(invcovm,xim02,xim01,cat)
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2m=jackknife_methods.jk_chi2(1./covm,xim02,xim01,cat)

    if linalg.cond(covp1) < 1/sys.float_info.epsilon:
      invcovp1=jackknife_methods.jk_corr_invcov(covp1,cat)
      chi2p1=jackknife_methods.jk_chi2(invcovp1,xip01,xip0,cat)
      xiperr1=1./np.sqrt(np.diagonal(invcovp1))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2p1=jackknife_methods.jk_chi2(1./covp1,xip01,xip0,cat)
      xiperr1=np.sqrt(np.diagonal(covp1))
    if linalg.cond(covm1) < 1/sys.float_info.epsilon:
      invcovm1=jackknife_methods.jk_corr_invcov(covm1,cat)
      chi2m1=jackknife_methods.jk_chi2(invcovm1,xim01,xim0,cat)
      ximerr1=1./np.sqrt(np.diagonal(invcovm1))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2m1=jackknife_methods.jk_chi2(1./covm1,xim01,xim0,cat)
      ximerr1=np.sqrt(np.diagonal(covm1))

    if linalg.cond(covp2) < 1/sys.float_info.epsilon:
      invcovp2=jackknife_methods.jk_corr_invcov(covp2,cat)
      chi2p2=jackknife_methods.jk_chi2(invcovp2,xip02,xip0,cat)
      xiperr2=1./np.sqrt(np.diagonal(invcovp2))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2p2=jackknife_methods.jk_chi2(1./covp2,xip02,xip0,cat)
      xiperr2=np.sqrt(np.diagonal(covp2))
    if linalg.cond(covm2) < 1/sys.float_info.epsilon:
      invcovm2=jackknife_methods.jk_corr_invcov(covm2,cat)
      chi2m2=jackknife_methods.jk_chi2(invcovm2,xim02,xim0,cat)
      ximerr2=1./np.sqrt(np.diagonal(invcovm2))
    else:
      print 'chi2 val likely wrong, using approx due to singular covariance matrix'
      chi2m2=jackknife_methods.jk_chi2(1./covm2,xim02,xim0,cat)
      ximerr2=np.sqrt(np.diagonal(covm2))


    if cat.xi_jk_st==0:
      cat.xi_jk_st=1



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
    plt.xlim((1,500))
    plt.ylim((5e-7,1e-4))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
    plt.savefig(CatalogMethods.cat_name(cat.cat)+'_jk-'+str(cat.use_jk)+'_nbins-'+str(cat.tbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_t0-'+str(cat.sep[0])+'_t1-'+str(cat.sep[1])+'_pzreweight-'+str(cat.use_zrw)+'_ztype-'+cat.ztyp+'_xi_'+label+'.png', bbox_inches='tight')
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
    cuts=CatalogMethods.add_cut(cuts,'ra',60.,noval,95.)
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)
    cuts=CatalogMethods.add_cut(cuts,'z_peak',0.3,noval,1.3)
    cuts=CatalogMethods.add_cut(cuts,'mean_psf_fwhm',0.,noval,7.4)

    return cuts

  @staticmethod
  def default_ngmix009_cuts():
    noval=999999

    cuts=CatalogMethods.add_cut(np.array([]),'exp_flags',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'flags_i',noval,noval,4)
    cuts=CatalogMethods.add_cut(cuts,'modest_class',noval,1,noval)
    cuts=CatalogMethods.add_cut(cuts,'exp_arate',0.4,noval,0.6)
    cuts=CatalogMethods.add_cut(cuts,'exp_s2n_w',10,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'exp_t_s2n',4,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',60.,noval,95.)
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)
    #cuts=CatalogMethods.add_cut(cuts,'exp_t',noval,noval,200)

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
  def get_cat_cols_i3(dir,cols,cuts,full,maxiter):
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

      if ifile==0:
        array=tmparray[mask]
      else:
        array=np.append(array,tmparray[mask],axis=0)

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

    print ' '
    print ' '
    print '..........'+label+'..........'
    print ' '
    print ' '
    theta,edge,edgemean,xip,xiperr,xim,ximerr,dxip0,dxip1,dxip2,xiperr0,xiperr1,xiperr2,chi2p0,chi2p1,chi2p2,dxim0,dxim1,dxim2,ximerr0,ximerr1,ximerr2,chi2m0,chi2m1,chi2m2=split_systematics.split_gals_2pt_along_(array,cat,label)
    #edge=edge.astype(int)
    A1b=SVA1.Fig_2ptPaper_amp(cat,xip,dxip1,xiperr1)
    A2b=SVA1.Fig_2ptPaper_amp(cat,xip,dxip2,xiperr2)
    A1bm=SVA1.Fig_2ptPaper_amp(cat,xim,dxim1,ximerr1)
    A2bm=SVA1.Fig_2ptPaper_amp(cat,xim,dxim2,ximerr2)
    
    cat.edgest=np.vstack((cat.edgest,edge))
    cat.meanst=np.vstack((cat.meanst,edgemean))
    cat.astp=np.vstack((cat.astp,np.array([A1a,A2a,A1b,A2b])))
    cat.astm=np.vstack((cat.astm,np.array([A1am,A2am,A1bm,A2bm])))
    cat.chi2stp=np.vstack((cat.chi2stp,np.array([chi2p0,chi2p1,chi2p2])))
    cat.chi2stm=np.vstack((cat.chi2stm,np.array([chi2m0,chi2m1,chi2m2])))

    print 'theta',theta
    print 'xip',xip
    print 'xiperr',xiperr
    print 'dxip1',dxip1
    print 'dxiperr1',xiperr1
    print 'dxip2',dxip2
    print 'dxiperr2',xiperr2
    print 'chi',chi2p1,chi2p2,chi2p0
    print 'Aa',A1a,A2a
    print 'Ab',A1b,A2b
    plt.figure(fig)
    ax=plt.subplot(r,c,n)
    ax.fill_between(theta,-xiperr/xip,xiperr/xip,facecolor='gray',alpha=0.4)
    plt.errorbar(theta,np.zeros((len(theta))),marker='',linestyle='-',color='k')
    plt.errorbar(theta,A1b*np.ones((len(theta))),marker='',linestyle='-',color='r')
    plt.errorbar(theta*(1),dxip1/xip,yerr=xiperr1/xip,marker='v',linestyle='',color='r')
    plt.errorbar(theta,A2b*np.ones((len(theta))),marker='',linestyle='-',color='b')
    plt.errorbar(theta*(1.2),dxip2/xip,yerr=xiperr2/xip,marker='^',linestyle='',color='b')
    t=ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=11,verticalalignment='top')
    t.set_bbox(dict(color='white', alpha=0.5, edgecolor='white'))
    #plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,500))
    plt.ylim((-1.49,1.49))
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
    if n<5:
      ax.set_xticklabels([])
    else:
      plt.xlabel(r'$\theta$ (arcmin)')
    if n%2==0:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\Delta\xi_+/\xi_+$')
    if n==1:
      ax.set_title('im3shape7.2')
    elif n==2:
      ax.set_title('ngmix010')

    plt.figure(fig+1)
    ax=plt.subplot(r,c,n)
    ax.fill_between(theta,-np.abs(ximerr/xim),np.abs(ximerr/xim),facecolor='gray',alpha=0.4)
    plt.errorbar(theta,np.zeros((len(theta))),marker='',linestyle='-',color='k')
    plt.errorbar(theta,A1bm*np.ones((len(theta))),marker='',linestyle='-',color='r')
    plt.errorbar(theta*(1),dxim1/xim,yerr=ximerr1/xim,marker='v',linestyle='',color='r')
    plt.errorbar(theta,A2bm*np.ones((len(theta))),marker='',linestyle='-',color='b')
    plt.errorbar(theta*(1.2),dxim2/xim,yerr=ximerr2/xim,marker='^',linestyle='',color='b')
    t=ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=11,verticalalignment='top')
    t.set_bbox(dict(color='white', alpha=0.5, edgecolor='white'))    
    #plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,500))
    plt.ylim((-1.49,1.49))
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
    if n<5:
      ax.set_xticklabels([])
    else:
      plt.xlabel(r'$\theta$ (arcmin)')
    if n%2==0:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\Delta\xi_{-}/\xi_{-}$')
    if n==1:
      ax.set_title('im3shape7.2')
    elif n==2:
      ax.set_title('ngmix010')
    return


  @staticmethod
  def Fig_2ptPaper_subplotold(cat,fig,r,c,n,array,label):

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
    #plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1,500))
    plt.ylim((1e-6,9e-4))
    #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
    if n<5:
      ax.set_xticklabels([])
    else:
      plt.xlabel(r'$\theta$ (arcmin)')
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
    lbinst=cat.lbins
    cat.lbins=100
    arr1,me1,me2,e1err,e2err,slp1,slp2,b1,b2=split_systematics.split_gals_lin_along_(array,cat,label)
    cat.lbins=1
    tmp,mm1,mm2,tmp,tmp,tmp,tmp,tmp,tmp=split_systematics.split_gals_lin_along_(array,cat,label)
    cat.lbins=lbinst
    arr1,me1,me2,e1err,e2err,tmp,tmp,tmp,tmp=split_systematics.split_gals_lin_along_(array,cat,label)
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
    plt.errorbar(arr1,slp1*arr1+b1,marker='',linestyle='-',color='r')
    plt.errorbar(arr1,slp2*arr1+b2,marker='',linestyle='-',color='b')
    plt.errorbar(arr1,mm1*np.ones((len(arr1))),marker='',linestyle=':',color='r')
    plt.errorbar(arr1,mm2*np.ones((len(arr1))),marker='',linestyle=':',color='b')
    plt.legend(loc='upper right',ncol=2, frameon=False,prop={'size':12})
    if n<4:
      ax.set_xticklabels([])
      plt.ylim((-.00175,.00175))
    else:
      plt.ylim((-.0005,.00175))
      if label=='psfe1':
        plt.xlabel(r'PSF $e_1$')
        plt.xlim((-.019,.0325))
      if label=='psfe2':
        plt.xlabel(r'PSF $e_2$')
        plt.xlim((-.019,.029))
      if label=='psffwhm':
        plt.xlabel(r'PSF FWHM')
        plt.xlim((.51,.7))
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
  def Fig_ShearPaper_subplot2(cat,ng,i3,fig,n,label):

    #if cat.cat==ng.cat:
    #  mask1=np.in1d(ng.coadd,i3.coadd,assume_unique=True)
    #  mask2=np.in1d(i3.coadd,ng.coadd,assume_unique=True)
    #  tmp=i3.w[mask2]
    #  sort=np.argsort(tmp)[np.argsort(np.argsort(ng.coadd[mask1]))]
    #  mask=mask1
    #  cat.w[mask1]=tmp[sort]
    #else:
    #  mask=np.ones((len(cat.coadd))).astype(bool)
    mask=np.ones((len(cat.coadd))).astype(bool)

    print ' '
    print ' '
    print '..........psf-shear '+label+'..........'
    print ' '
    print ' '
    theta,ggp,ggm,ggperr,ggmerr,gpp,gpm,gpperr,gpmerr,ppp,ppm,ppperr,ppmerr,alphap,alphaperr,leakagep,leakageperr,alpham,alphamerr,leakagem,leakagemerr=xi_2pt_shear_test_methods.xi_2pt_psf(cat,mask)
    print 'ggp',ggp
    print 'ggperr',ggperr
    print 'gpp',gpp
    print 'gpperr',gpperr
    print 'ppp',ppp
    print 'ppperr',ppperr
    print 'alpha',alphap,alpham
    print 'leakage',leakagep,leakagem

    plt.figure(fig)

    ax=plt.subplot(3,2,1+n)
    if n==0:
      ax.set_title('im3shapev7.2')
    else:
      ax.set_title('ngmix009')
    plt.errorbar(theta,ggp,yerr=ggperr,marker='.',linestyle='',color='k',label='gg')
    plt.errorbar(theta,np.abs(gpp),yerr=gpperr,marker='.',linestyle='',color='r',label='gp')
    plt.errorbar(theta,ppp,yerr=ppperr,marker='.',linestyle='',color='b',label='pp')
    plt.legend(loc='upper right',ncol=3, frameon=False,prop={'size':12})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((1,500))
    plt.ylim(5e-7,6e-4)
    ax.set_xticklabels([])
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\xi_+$')

    ax=plt.subplot(3,2,3+n)
    plt.errorbar(theta,alphap,marker='',linestyle='-',color='k')
    ax.fill_between(theta,alphap-alphaperr,alphap+alphaperr,facecolor='gray',alpha=0.5)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim((1,500))
    plt.ylim(-.1,.4)
    ax.set_xticklabels([])
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\alpha$')

    ax=plt.subplot(3,2,5+n)
    plt.errorbar(theta,leakagep,yerr=leakageperr,marker='',linestyle='-',color='k')
    ax.fill_between(theta,leakagep-leakageperr,leakagep+leakageperr,facecolor='gray',alpha=0.5)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim((1,500))
    plt.ylim(-.1,.5)
    plt.xlabel(r'$\theta$ (arcmin)')
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'Leakage into $\xi_+$')

    plt.figure(fig+1)

    ax=plt.subplot(3,2,1+n)
    plt.errorbar(theta,ggm,yerr=ggmerr,marker='.',linestyle='',color='k',label='gg')
    plt.errorbar(theta,gpm,yerr=gpmerr,marker='.',linestyle='',color='r',label='gp')
    plt.errorbar(theta,ppm,yerr=ppmerr,marker='.',linestyle='',color='b',label='pp')
    plt.legend(loc='upper right',ncol=3, frameon=False,prop={'size':12})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((1,500))
    plt.ylim(5e-7,6e-4)
    ax.set_xticklabels([])
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\xi_-$')

    ax=plt.subplot(3,2,3+n)
    plt.errorbar(theta,alpham,marker='',linestyle='-',color='k')
    ax.fill_between(theta,alpham-alphamerr,alpham+alphamerr,facecolor='gray',alpha=0.5)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim((1,500))
    plt.ylim(-.1,.49)
    ax.set_xticklabels([])
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'$\alpha$')

    ax=plt.subplot(3,2,5+n)
    plt.errorbar(theta,leakagem,yerr=leakagemerr,marker='',linestyle='-',color='k')
    ax.fill_between(theta,leakagem-leakagemerr,leakagem+leakagemerr,facecolor='gray',alpha=0.5)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim((1,500))
    plt.ylim(-.1,.49)
    plt.xlabel(r'$\theta$ (arcmin)')
    if n==1:
      ax.set_yticklabels([])
    else:
      plt.ylabel(r'Leakage into $\xi_-$')

    return


  @staticmethod
  def Fig1_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 1'
    print '...........'

    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,1,i3.snr,'Signal-to-Noise')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,2,ng.snr,'Signal-to-Noise')
    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,3,i3.radius,'Size')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,4,ng.radius,'Size')
    SVA1.Fig_2ptPaper_subplot(i3,1,3,2,5,i3.colour,'Colour')
    SVA1.Fig_2ptPaper_subplot(ng,1,3,2,6,ng.colour,'Colour')
    plt.figure(1)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig1p.png', bbox_inches='tight')
    plt.close(1)   
    plt.figure(2)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig1m.png', bbox_inches='tight')
    plt.close(2)  

    print 'edge final',i3.edgest,ng.edgest
    print 'mean final',i3.meanst,ng.meanst
    print 'Astp final',i3.astp,ng.astp
    print 'Astm final',i3.astm,ng.astm
    print 'chi2p final',i3.chi2stp,ng.chi2stp
    print 'chi2m final',i3.chi2stm,ng.chi2stm

    tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
    tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
    tmpi3=np.around(tmpi3, decimals=2)
    tmpng=np.around(tmpng, decimals=2)
    np.set_printoptions(suppress=True)
    print tmpi3
    print tmpng    


    return

  @staticmethod
  def Fig2_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 2'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,1,i3.ra,'RA')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,2,ng.ra,'RA')
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,3,i3.dec,'Dec')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,4,ng.dec,'Dec')
    SVA1.Fig_2ptPaper_subplot(i3,2,3,2,5,i3.ebv,'E(B-V)')
    SVA1.Fig_2ptPaper_subplot(ng,2,3,2,6,ng.ebv,'E(B-V)')
    plt.figure(2)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig2p.png', bbox_inches='tight')
    plt.close(2)   
    plt.figure(3)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig2m.png', bbox_inches='tight')
    plt.close(3)   

    print 'edge final',i3.edgest,ng.edgest
    print 'mean final',i3.meanst,ng.meanst
    print 'Astp final',i3.astp,ng.astp
    print 'Astm final',i3.astm,ng.astm
    print 'chi2p final',i3.chi2stp,ng.chi2stp
    print 'chi2m final',i3.chi2stm,ng.chi2stm

    tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
    tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
    tmpi3=np.around(tmpi3, decimals=2)
    tmpng=np.around(tmpng, decimals=2)
    np.set_printoptions(suppress=True)
    print tmpi3
    print tmpng    


    return

  @staticmethod
  def Fig3_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 3'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,1,i3.exptime,'Exposure Time')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,2,ng.exptime,'Exposure Time')
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,3,i3.maglimit,'Mag Limit')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,4,ng.maglimit,'Mag Limit')
    SVA1.Fig_2ptPaper_subplot(i3,3,3,2,5,i3.airmass,'Air Mass')
    SVA1.Fig_2ptPaper_subplot(ng,3,3,2,6,ng.airmass,'Air Mass')
    plt.figure(3)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig3p.png', bbox_inches='tight')
    plt.close(3)   
    plt.figure(4)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig3m.png', bbox_inches='tight')
    plt.close(4)  

    print 'edge final',i3.edgest,ng.edgest
    print 'mean final',i3.meanst,ng.meanst
    print 'Astp final',i3.astp,ng.astp
    print 'Astm final',i3.astm,ng.astm
    print 'chi2p final',i3.chi2stp,ng.chi2stp
    print 'chi2m final',i3.chi2stm,ng.chi2stm

    tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
    tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
    tmpi3=np.around(tmpi3, decimals=2)
    tmpng=np.around(tmpng, decimals=2)
    np.set_printoptions(suppress=True)
    print tmpi3
    print tmpng    

    return


  @staticmethod
  def Fig4_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 4'
    print '...........'

  
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,1,i3.skysigma,'Sky Sigma')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,2,ng.skysigma,'Sky Sigma')
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,3,i3.skybrite,'Sky Brightness')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,4,ng.skybrite,'Sky Brightness')
    SVA1.Fig_2ptPaper_subplot(i3,4,3,2,5,i3.fwhm,'FWHM')
    SVA1.Fig_2ptPaper_subplot(ng,4,3,2,6,ng.fwhm,'FWHM')
    plt.figure(4)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig4p.png', bbox_inches='tight')
    plt.close(4)   
    plt.figure(5)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig4m.png', bbox_inches='tight')
    plt.close(5)   
        
    print 'edge final',i3.edgest,ng.edgest
    print 'mean final',i3.meanst,ng.meanst
    print 'Astp final',i3.astp,ng.astp
    print 'Astm final',i3.astm,ng.astm
    print 'chi2p final',i3.chi2stp,ng.chi2stp
    print 'chi2m final',i3.chi2stm,ng.chi2stm

    tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
    tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
    tmpi3=np.around(tmpi3, decimals=2)
    tmpng=np.around(tmpng, decimals=2)
    np.set_printoptions(suppress=True)
    print tmpi3
    print tmpng    


    return

  @staticmethod
  def Fig5_2ptPaper(i3,ng):


    print '...........'
    print 'Figure 5'
    print '...........'

    
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,1,i3.psf1,r'PSF $e_1$')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,2,ng.psf1,r'PSF $e_1$')
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,3,i3.psf2,r'PSF $e_2$')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,4,ng.psf2,r'PSF $e_2$')
    SVA1.Fig_2ptPaper_subplot(i3,5,3,2,5,i3.psffwhm,'PSF FWHM')
    SVA1.Fig_2ptPaper_subplot(ng,5,3,2,6,ng.psffwhm,'PSF FWHM')
    plt.figure(5)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig5p.png', bbox_inches='tight')
    plt.close(5)   
    plt.figure(6)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig5m.png', bbox_inches='tight')
    plt.close(6)   
        
    print 'edge final',i3.edgest,ng.edgest
    print 'mean final',i3.meanst,ng.meanst
    print 'Astp final',i3.astp,ng.astp
    print 'Astm final',i3.astm,ng.astm
    print 'chi2p final',i3.chi2stp,ng.chi2stp
    print 'chi2m final',i3.chi2stm,ng.chi2stm

    tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
    tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
    tmpi3=np.around(tmpi3, decimals=2)
    tmpng=np.around(tmpng, decimals=2)
    np.set_printoptions(suppress=True)
    print tmpi3
    print tmpng    


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

  @staticmethod
  def Fig2_ShearPaper(i3,ng):


    print '...........'
    print 'Figure 2'
    print '...........'

    SVA1.Fig_ShearPaper_subplot2(i3,ng,i3,1,0,'im3shape')
    SVA1.Fig_ShearPaper_subplot2(ng,ng,i3,1,1,'ngmix')
    plt.figure(1)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_wt-'+str(i3.wt)+'_Fig2_p.png', bbox_inches='tight')
    plt.close(1)
    plt.figure(2)
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_Fig2_m.png', bbox_inches='tight')
    plt.close(2)
    
    return

