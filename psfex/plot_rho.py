#! /usr/bin/env python
# Program to plot rho statistics on PSFEx outputs.

import astropy.io.fits as pyfits
import numpy
import json
import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os
import sys

plt.style.use('/astro/u/mjarvis/.config/matplotlib/stylelib/supermongo.mplstyle')
 
def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits*',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')

    args = parser.parse_args()
    return args


def parse_file_name(file_name):
    """Parse the PSFEx file name to get the root name and the chip number
    """
    import os

    base_file = os.path.split(file_name)[1]
    if os.path.splitext(base_file)[1] == '.fz':
        base_file=os.path.splitext(base_file)[0]
    if os.path.splitext(base_file)[1] != '.fits':
        raise ValueError("Invalid file name "+file)
    root = os.path.splitext(base_file)[0]

    ccdnum = int(root.split('_')[-1])
    return root, ccdnum

def plot_rho(meanr, rhop, sigp, sqrtn, rhom=None, sigm=None):
    print 'meanr = ',meanr
    print 'rhop = ',rhop
    print 'sigp = ',sigp
    plt.plot(meanr, rhop, color='blue')
    plt.plot(meanr, -rhop, color='blue', ls=':')
    plt.errorbar(meanr[rhop>0], rhop[rhop>0], yerr=sigp[rhop>0]/sqrtn, color='blue', ls='')
    plt.errorbar(meanr[rhop>0], rhop[rhop>0], yerr=sigp[rhop>0], color='blue', lw=0.1, ls='')
    plt.errorbar(meanr[rhop<0], -rhop[rhop<0], yerr=sigp[rhop<0]/sqrtn, color='blue', ls='')
    plt.errorbar(meanr[rhop<0], -rhop[rhop<0], yerr=sigp[rhop<0], color='blue', lw=0.1, ls='')
    lp = plt.errorbar(-meanr, rhop, yerr=sigp, color='blue')

    if rhom is not None:
        plt.plot(meanr, rhom, color='green')
        plt.plot(meanr, -rhom, color='green', ls=':')
        plt.errorbar(meanr[rhom>0], rhom[rhom>0], yerr=sigm[rhom>0]/sqrtn, color='green', ls='')
        plt.errorbar(meanr[rhom>0], rhom[rhom>0], yerr=sigm[rhom>0], color='green', lw=0.1, ls='')
        plt.errorbar(meanr[rhom<0], -rhom[rhom<0], yerr=sigm[rhom<0]/sqrtn, color='green', ls='')
        plt.errorbar(meanr[rhom<0], -rhom[rhom<0], yerr=sigm[rhom<0], color='green', lw=0.1, ls='')
        lm = plt.errorbar(-meanr, rhom, yerr=sigm, color='green')

    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.xlabel(r'$\theta$ (arcmin)')

    if rhom is not None:
        return [ lp, lm ]
    else:
        return [ lp ]

def pretty_rho1(meanr, rho, sig, sqrtn, rho3=None, rho4=None):
    import matplotlib.patches as mp
    if False:
        # This is all handwavy arguments about what the requirements are.  
        # I'm taking Cathering's CFHTLS xi+ values of 1.e-4 at 1 arcmin, 2e-6 at 40 arcmin.
        # Then I'm saying our requirements on rho need to be about 0.16 times this for SV (S/N=6),
        # but more like 0.03 times this for Y5.
        t1 = 0.5
        t2 = 300.
        xa = numpy.log(1.0)
        xb = numpy.log(40.0)
        ya = numpy.log(1.e-4)
        yb = numpy.log(2.e-6)
        x1 = numpy.log(t1)
        x2 = numpy.log(t2)
        y1 = (yb * (x1-xa) + ya * (xb-x1)) / (xb-xa)
        y2 = (yb * (x2-xa) + ya * (xb-x2)) / (xb-xa)
        xi1 = numpy.exp(y1)
        xi2 = numpy.exp(y2)

        sv_req = plt.fill( [t1, t1, t2, t2], [0., xi1 * 0.16, xi2 * 0.16, 0.], 
                           color = '#FFFF82')
        y5_req = plt.fill( [t1, t1, t2, t2], [0., xi1 * 0.03, xi2 * 0.03, 0.], 
                           color = '#BAFFA4')
    else:
        # Use Alexandre's xi file as a requirement.
        theta = [0.5, 1.49530470404, 2.91321328166, 5.67019971545, 11.0321639144,
                 21.4548924111, 41.6931936543, 80.8508152859, 156.285886576,
                 297.92139021, 300.]

        dxi = [1.4e-5, 3.82239654447e-06, 2.36185315415e-06, 1.26849547074e-06, 6.3282672138e-07,
               3.25623661098e-07, 1.747852053e-07, 8.75326181278e-08, 3.60247306537e-08,
               1.13521735321e-08, 1.125e-8]
        req = [ x / 5.9 for x in dxi ]
        sv_req = plt.fill( [theta[0]] + theta + [theta[-1]], [0.] + req + [0.],
                           color = '#FFFF82')
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig[rho3>0]/sqrtn, color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig[rho3<0]/sqrtn, color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig, color='green', marker='s')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig[rho4>0]/sqrtn, color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig[rho4<0]/sqrtn, color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig, color='red', marker='^')
    sv_req = mp.Patch(color='#FFFF82')
    if rho3 is not None and rho4 is not None:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$', r'$\rho_4(\theta)$'],
                   loc='upper right', fontsize=24)
        plt.ylim( [1.e-9, 5.e-6] )
    elif True:
        plt.legend([rho1_line, sv_req],
                   [r'$\rho_1(\theta)$', r'Requirement'],
                   loc='upper right')
        plt.ylim( [1.e-9, 5.e-6] )
    else: # For talk
        plt.legend([rho1_line, sv_req],
                   [r'$\rho_1(\theta)$',
                    r'Requirements for $d\sigma_8/\sigma_8 < 0.03$'],
                   loc='upper right')
        plt.ylim( [1.e-9, 3.e-6] )
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\rho(\theta)$')
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, sqrtn, rho5=None):
    import matplotlib.patches as mp
    # The requirements on rho2 are less stringent.  They are larger by a factor 1/alpha.
    # Let's use alpha = 0.03.
    alpha = 0.03
    if False:
        t1 = 0.5
        t2 = 300.
        xa = numpy.log(1.0)
        xb = numpy.log(40.0)
        ya = numpy.log(1.e-4)
        yb = numpy.log(2.e-6)
        x1 = numpy.log(t1)
        x2 = numpy.log(t2)
        y1 = (yb * (x1-xa) + ya * (xb-x1)) / (xb-xa)
        y2 = (yb * (x2-xa) + ya * (xb-x2)) / (xb-xa)
        xi1 = numpy.exp(y1)
        xi2 = numpy.exp(y2)

        sv_req = plt.fill( [t1, t1, t2, t2], [0., xi1 * 0.16 / alpha, xi2 * 0.16 / alpha, 0.], 
                           color = '#FFFF82')
        y5_req = plt.fill( [t1, t1, t2, t2], [0., xi1 * 0.03 / alpha, xi2 * 0.03 / alpha, 0.], 
                           color = '#BAFFA4')

    else:
        # Use Alexandre's xi file as a requirement.
        theta = [0.5, 1.49530470404, 2.91321328166, 5.67019971545, 11.0321639144,
                 21.4548924111, 41.6931936543, 80.8508152859, 156.285886576,
                 297.92139021, 300.]

        dxi = [1.4e-5, 3.82239654447e-06, 2.36185315415e-06, 1.26849547074e-06, 6.3282672138e-07,
               3.25623661098e-07, 1.747852053e-07, 8.75326181278e-08, 3.60247306537e-08,
               1.13521735321e-08, 1.125e-8]
        req = [ x / (2.4 * alpha) for x in dxi ]
        sv_req = plt.fill( [theta[0]] + theta + [theta[-1]], [0.] + req + [0.],
                           color = '#FFFF82')
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho2_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='green')
        plt.plot(meanr*1.03, -rho5, color='green', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig[rho5>0]/sqrtn, color='green', ls='', marker='s')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig[rho5<0]/sqrtn, color='green', ls='', marker='s')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig, color='green', marker='s')
    sv_req = mp.Patch(color='#FFFF82')
    if rho5 is not None:
        plt.legend([rho2_line, rho5_line],
                   [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                   loc='upper right', fontsize=24)
        plt.ylim( [1.e-7, 5.e-4] )
    elif True: # For paper
        plt.legend([rho2_line, sv_req],
                   [r'$\rho_2(\theta)$', r'Requirement'],
                   loc='upper right')
        plt.ylim( [1.e-7, 5.e-4] )
    else:
        plt.legend([rho2_line, sv_req],
                   [r'$\rho_2(\theta)$',
                    r'Requirements for $d\sigma_8/\sigma_8 < 0.03$'],
                   loc='upper right')
        plt.ylim( [1.e-7, 3.e-4] )
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\rho(\theta)$')
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho3(meanr, rho, sig, sqrtn):
    import matplotlib.patches as mp
    t1 = 0.5
    t2 = 300.
    y_sv = 1.1e-2**2
    #y_y5 = 2.2e-3**2

    sv_req = plt.fill( [t1, t1, t2, t2], [0., y_sv, y_sv, 0.],
                        color = '#FFFF82')
    #y5_req = plt.fill( [t1, t1, t2, t2], [0., y_y5, y_y5, 0.],
                        #color = '#BAFFA4')
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='')
    rho3_line, = plt.plot(-meanr, rho, color='blue')
    sv_req = mp.Patch(color='#FFFF82')
    #y5_req = mp.Patch(color='#BAFFA4')
    plt.legend([rho3_line, sv_req],
               [r'$\rho_3(\theta)$', 'SV Requirements'])
    plt.xlim( [0.5,300.] )
    plt.ylim( [3.e-6, 3.e-4] )
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\rho_3$')
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()



def main():
    import os
    import glob
    import galsim

    args = parse_args()

    datadir = '/astro/u/astrodat/data/DES'

    work = os.path.expanduser(args.work)
    print 'work dir = ',work

    if args.file != '':
        print 'Read file ',args.file
        with open(args.file) as fin:
            data = [ line.split() for line in fin ]
        runs, exps = zip(*data)
    else:
        runs = args.runs
        exps = args.exps

    nexp = len(exps)

    if False:
        ccd_meanlogr = numpy.empty( (nexp*62,37) )
        ccd_rho1p = numpy.empty( (nexp*62,37) )
        ccd_rho1m = numpy.empty( (nexp*62,37) )
        ccd_rho2p = numpy.empty( (nexp*62,37) )
        ccd_rho2m = numpy.empty( (nexp*62,37) )
        ccd_rho3 = numpy.empty( (nexp*62,37) )
        exp_meanlogr = numpy.empty( (nexp,53) )
        exp_rho1p = numpy.empty( (nexp,53) )
        exp_rho1m = numpy.empty( (nexp,53) )
        exp_rho2p = numpy.empty( (nexp,53) )
        exp_rho2m = numpy.empty( (nexp,53) )
        exp_rho3 = numpy.empty( (nexp,53) )
        exp_var1 = numpy.empty( (nexp,53) )
        exp_var2 = numpy.empty( (nexp,53) )
        exp_var3 = numpy.empty( (nexp,53) )
        exp_var4 = numpy.empty( (nexp,53) )
        desdm_meanlogr = numpy.empty( (nexp,53) )
        desdm_rho1p = numpy.empty( (nexp,53) )
        desdm_rho1m = numpy.empty( (nexp,53) )
        desdm_rho2p = numpy.empty( (nexp,53) )
        desdm_rho2m = numpy.empty( (nexp,53) )
        desdm_rho3 = numpy.empty( (nexp,53) )
        desdm_var1 = numpy.empty( (nexp,53) )
        desdm_var2 = numpy.empty( (nexp,53) )
        desdm_var3 = numpy.empty( (nexp,53) )
        desdm_var4 = numpy.empty( (nexp,53) )

        meande1 = 0
        meande2 = 0
        varde1 = 0
        varde2 = 0
        nde = 0
        histde1 = numpy.zeros(200)  # bin size = 1.e-3
        histde2 = numpy.zeros(200)  # bin size = 1.e-3
        histnstars = numpy.zeros(200)  # bin size = 10
        listnstars = []
        meannstars = 0
        ngoodccd = 0

        iexp = 0
        iccd = 0
        for run,exp in zip(runs,exps):

            print 'Start work on run, exp = ',run,exp
            expnum = int(exp[6:])
            print 'expnum = ',expnum,'  ',iexp,'/',nexp,'  ',iccd

            exp_dir = os.path.join(work,exp)

            cat_file = os.path.join(exp_dir, exp + "_psf.fits")
            with pyfits.open(cat_file) as pyf:
                data = pyf[1].data
            ccdnums = numpy.unique(data['ccdnum'])
            for ccdnum in ccdnums:
                nstars = ((data['ccdnum'] == ccdnum) & (data['flag'] == 0)).sum()
                if nstars > 0 and nstars < 2000:
                    histnstars[ int(numpy.floor(nstars/10)) ] += 1
                if nstars > 0:
                    meannstars += nstars
                    ngoodccd += 1
                    listnstars.append(nstars)
            mask = data['flag'] == 0
            de1 = data['e1'][mask] - data['psfex_e1'][mask]
            de2 = data['e2'][mask] - data['psfex_e2'][mask]
            meande1 += numpy.sum(de1)
            meande2 += numpy.sum(de2)
            varde1 += numpy.sum(de1*de1)
            varde2 += numpy.sum(de2*de2)
            nde += len(de1)
            histde1 += numpy.histogram(de1, bins=200, range=(-1.e-1,1.e-1))[0]
            histde2 += numpy.histogram(de2, bins=200, range=(-1.e-1,1.e-1))[0]

            stat_file = os.path.join(exp_dir, exp + ".json")

            # Read the json file 
            if not os.path.exists(stat_file):
                print stat_file,' not found'
                print 'No JSON file for this exposure.  Skipping.'
                continue
            with open(stat_file,'r') as f:
                stats = json.load(f)

            print "len stats = ",len(stats)
            ( expnum, 
              rho1_meanlogr,
              rho1_xip,
              rho1_xip_im,
              rho1_xim,
              rho1_xim_im,
              rho1_varxi,
              rho2_xip,
              rho2_xip_im,
              rho2_xim,
              rho2_xim_im,
              rho2_varxi,
              rho3_xi,
              rho3_varxi,
              drho1_meanlogr,
              drho1_xip,
              drho1_xip_im,
              drho1_xim,
              drho1_xim_im,
              drho1_varxi,
              drho2_xip,
              drho2_xip_im,
              drho2_xim,
              drho2_xim_im,
              drho2_varxi,
              drho3_xi,
              drho3_varxi,
            ) = stats[-1]
            exp_meanlogr[iexp,:] = rho1_meanlogr
            exp_rho1p[iexp,:] = rho1_xip
            exp_rho1m[iexp,:] = rho1_xim
            exp_rho2p[iexp,:] = rho2_xip
            exp_rho2m[iexp,:] = rho2_xim
            exp_rho3[iexp,:] = rho3_xi
            exp_var1[iexp,:] = rho1_varxi
            exp_var2[iexp,:] = rho2_varxi
            exp_var3[iexp,:] = rho3_varxi
            desdm_meanlogr[iexp,:] = drho1_meanlogr
            desdm_rho1p[iexp,:] = drho1_xip
            desdm_rho1m[iexp,:] = drho1_xim
            desdm_rho2p[iexp,:] = drho2_xip
            desdm_rho2m[iexp,:] = drho2_xim
            desdm_rho3[iexp,:] = drho3_xi
            desdm_var1[iexp,:] = drho1_varxi
            desdm_var2[iexp,:] = drho2_varxi
            desdm_var3[iexp,:] = drho3_varxi
            iexp += 1
 
            for s in stats[:-1]:

                ( ccdnum, 
                  rho1_meanlogr,
                  rho1_xip,
                  rho1_xim,
                  rho2_xip,
                  rho2_xim,
                  rho3_xi,
                ) = s

                ccd_meanlogr[iccd,:] = rho1_meanlogr
                ccd_rho1p[iccd,:] = rho1_xip
                ccd_rho1m[iccd,:] = rho1_xim
                ccd_rho2p[iccd,:] = rho2_xip
                ccd_rho2m[iccd,:] = rho2_xim
                ccd_rho3[iccd,:] = rho3_xi
                iccd += 1

        print '\nFinished processing all exposures'
        nexp = iexp

        nccd = iccd
        # Compute some stats and plot histograms
        meande1 /= nde
        meande2 /= nde
        varde1 -= nde * meande1**2
        varde2 -= nde * meande2**2
        varde1 /= nde
        varde2 /= nde
        print 'nde = ',nde
        print 'mean de1 = ',meande1
        print 'sigma = ',numpy.sqrt(varde1)
        print 'mean de2 = ',meande2
        print 'sigma = ',numpy.sqrt(varde2)

        plt.clf()
        left = [ (i-100) * 1.e-3 for i in range(200) ]
        plt.bar( left, histde1, width=1.e-3 )
        plt.xlim( [-1.e-1,1.e-1] )
        plt.xlabel(r'$e1_{psf} - e1_{model}$')
        plt.ylabel(r'$N_{stars}$')
        plt.title('Distribution of PSF e1 residuals')
        import matplotlib.mlab as mlab
        plt.plot(left,numpy.sum(histde1)*1.e-3*mlab.normpdf(left,meande1,numpy.sqrt(varde1)),
                 color='red')
        #plt.savefig('de1hist.png')
        plt.savefig('de1hist.pdf')

        plt.clf()
        plt.bar( left, histde2, width=1.e-3 )
        plt.xlim( [-1.e-1,1.e-1] )
        plt.xlabel(r'$e2_{psf} - e2_{model}$')
        plt.ylabel(r'$N_{stars}$')
        plt.title('Distribution of PSF e2 residuals')
        plt.plot(left,numpy.sum(histde2)*1.e-3*mlab.normpdf(left,meande2,numpy.sqrt(varde2)),
                 color='red')
        #plt.savefig('de2hist.png')
        plt.savefig('de2hist.pdf')

        plt.clf()
        left = [ 10*i for i in range(200) ]
        plt.bar( left, histnstars, width=10 )
        plt.xlim( [0,700] )
        plt.xlabel(r'stars per CCD')
        plt.ylabel(r'$N_\mathrm{CCDs}$')
        #plt.title('Distribution of number of stars per chip')
        #plt.savefig('nstarshist.png')
        plt.tight_layout()
        plt.savefig('nstarshist.pdf')

        imax = numpy.argmax(histnstars)
        meannstars /= ngoodccd
        print 'mode nstars = ',(imax+0.5)*10.
        print 'mean nstars = ',meannstars
        listnstars.sort()
        print 'mean nstars = ',sum(listnstars)/len(listnstars)
        print 'median nstars = ',listnstars[len(listnstars)//2]


    if False:
        # Plots for CCDs
        print 'nccd = ',nccd
        sqrtn = numpy.sqrt(nccd)
        meanr = numpy.exp(numpy.mean(ccd_meanlogr[:nccd,:], axis=0))
        rho1p = numpy.mean(ccd_rho1p[:nccd,:], axis=0)
        rho1m = numpy.mean(ccd_rho1m[:nccd,:], axis=0)
        rho2p = numpy.mean(ccd_rho2p[:nccd,:], axis=0)
        rho2m = numpy.mean(ccd_rho2m[:nccd,:], axis=0)
        rho3 = numpy.mean(ccd_rho3[:nccd,:], axis=0)
        sig_rho1p = numpy.std(ccd_rho1p[:nccd,:], axis=0)
        sig_rho1m = numpy.std(ccd_rho1m[:nccd,:], axis=0)
        sig_rho2p = numpy.std(ccd_rho2p[:nccd,:], axis=0)
        sig_rho2m = numpy.std(ccd_rho2m[:nccd,:], axis=0)
        sig_rho3 = numpy.std(ccd_rho3[:nccd,:], axis=0)
        print 'meanr = ',meanr
        print 'rho1p = ',rho1p
        print 'sig_rho1p = ',sig_rho1p
        plt.rc('font', family='serif')

        plt.clf()
        plt.title(r'SPTE $\rho_1$ (i.e. $\langle de de \rangle$) for individual CCDs')
        lines = plot_rho(meanr, rho1p, sig_rho1p, sqrtn, rho1m, sig_rho1m)
        plt.legend(lines, [r'$\rho_1(\theta)+$', r'$\rho_1(\theta)-$'] )
        plt.xlim( [0.5,20] )
        plt.ylabel(r'$\rho_1$')
        #plt.savefig('ccd_rho1.png')
        plt.savefig('ccd_rho1.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_2$ (i.e. $\langle e de \rangle$) for individual CCDs')
        lines = plot_rho(meanr, rho2p, sig_rho2p, sqrtn, rho2m, sig_rho2m)
        plt.legend(lines, [r'$\rho_2(\theta)+$', r'$\rho_2(\theta)-$'] )
        plt.xlim( [0.5,20] )
        plt.ylabel(r'$\rho_2$')
        #plt.savefig('ccd_rho2.png')
        plt.savefig('ccd_rho2.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_3$ (i.e. $\langle ds ds \rangle$) for individual CCDs')
        lines = plot_rho(meanr, rho3, sig_rho3, sqrtn)
        plt.legend(lines, [r'$\rho_3(\theta)$'] )
        plt.xlim( [0.5,20] )
        plt.ylabel(r'$\rho_3$')
        #plt.savefig('ccd_rho3.png')
        plt.savefig('ccd_rho3.pdf')

    if False:
        # Plots for exposures:
        print 'nexp = ',nexp
        sqrtn = numpy.sqrt(nexp)
        meanr = numpy.exp(numpy.mean(exp_meanlogr[:nexp,:], axis=0))
        rho1p = numpy.mean(exp_rho1p[:nexp,:], axis=0)
        rho1m = numpy.mean(exp_rho1m[:nexp,:], axis=0)
        rho2p = numpy.mean(exp_rho2p[:nexp,:], axis=0)
        rho2m = numpy.mean(exp_rho2m[:nexp,:], axis=0)
        rho3 = numpy.mean(exp_rho3[:nexp,:], axis=0)
        sig_rho1p = numpy.std(exp_rho1p[:nexp,:], axis=0)
        sig_rho1m = numpy.std(exp_rho1m[:nexp,:], axis=0)
        sig_rho2p = numpy.std(exp_rho2p[:nexp,:], axis=0)
        sig_rho2m = numpy.std(exp_rho2m[:nexp,:], axis=0)
        sig_rho3 = numpy.std(exp_rho3[:nexp,:], axis=0)
        print 'meanr = ',meanr
        print 'rho1p = ',rho1p
        print 'sig_rho1p = ',sig_rho1p
        plt.rc('font', family='serif')

        plt.clf()
        plt.title(r'SPTE $\rho_1$ (i.e. $\langle de de \rangle$) for full exposures')
        lines = plot_rho(meanr, rho1p, sig_rho1p, sqrtn, rho1m, sig_rho1m)
        plt.legend(lines, [r'$\rho_1(\theta)+$', r'$\rho_1(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_1$')
        #plt.savefig('exp_rho1.png')
        plt.savefig('exp_rho1.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_2$ (i.e. $\langle e de \rangle$) for full exposures')
        lines = plot_rho(meanr, rho2p, sig_rho2p, sqrtn, rho2m, sig_rho2m)
        plt.legend(lines, [r'$\rho_2(\theta)+$', r'$\rho_2(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_2$')
        #plt.savefig('exp_rho2.png')
        plt.savefig('exp_rho2.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_3$ (i.e. $\langle ds ds \rangle$) for full exposures')
        lines = plot_rho(meanr, rho3, sig_rho3, sqrtn)
        plt.legend(lines, [r'$\rho_3(\theta)$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_3$')
        #plt.savefig('exp_rho3.png')
        plt.savefig('exp_rho3.pdf')

        # Prettier plots for Erin's talk
        plt.clf()
        pretty_rho1(meanr, rho1p, sig_rho1p, sqrtn)
        plt.savefig('rho1.pdf')
        #plt.savefig('rho1.png')

        plt.clf()
        pretty_rho2(meanr, rho2p, sig_rho2p, sqrtn)
        plt.savefig('rho2.pdf')
        #plt.savefig('rho2.png')

        plt.clf()
        pretty_rho3(meanr, rho3, sig_rho3, sqrtn)
        plt.savefig('rho3.pdf')
        #plt.savefig('rho3.png')

    k10arcmin = int(round(numpy.log(10 / 0.5)/0.1))
    if False:
        # Plots for worst rho1 exposure:
        # Find worst exposure based on rho2 at theta = 10 arcmin
        i = numpy.argmax(numpy.abs(exp_rho1p[:nexp,k10arcmin]), axis=0)
        print 'k10arcmin = ',k10arcmin
        print 'rho1[k] = ',exp_rho1p[:nexp,k10arcmin]
        print 'rho2[k] = ',exp_rho2p[:nexp,k10arcmin]
        print 'i = ',i
        print 'rho1[i] = ',exp_rho1p[i,:]
        meanr = numpy.exp(exp_meanlogr[i,:])
        rho1p = exp_rho1p[i,:]
        rho1m = exp_rho1m[i,:]
        rho2p = exp_rho2p[i,:]
        print 'rho2p = ',rho2p
        rho2m = exp_rho2m[i,:]
        rho3 = exp_rho3[i,:]
        sig_rho1p = numpy.sqrt(exp_var1[i,:])
        sig_rho1m = numpy.sqrt(exp_var1[i,:])
        sig_rho2p = numpy.sqrt(exp_var2[i,:])
        sig_rho2m = numpy.sqrt(exp_var2[i,:])
        sig_rho3 = numpy.sqrt(exp_var3[i,:])

        plt.clf()
        plt.title(r'$\rho_1$ for exposure with worst $\rho_1$ at 10 arcmin')
        lines = plot_rho(meanr, rho1p, sig_rho1p, 1, rho1m, sig_rho1m)
        plt.legend(lines, [r'$\rho_1(\theta)+$', r'$\rho_1(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_1$')
        #plt.savefig('w1_rho1.png')
        plt.savefig('w1_rho1.pdf')

        plt.clf()
        plt.title(r'$\rho_2$ for exposure with worst $\rho_1$ at 10 arcmin')
        lines = plot_rho(meanr, rho2p, sig_rho2p, 1, rho2m, sig_rho2m)
        plt.legend(lines, [r'$\rho_2(\theta)+$', r'$\rho_2(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_2$')
        #plt.savefig('w1_rho2.png')
        plt.savefig('w1_rho2.pdf')

        plt.clf()
        plt.title(r'$\rho_3$ for exposure with worst $\rho_1$ at 10 arcmin')
        lines = plot_rho(meanr, rho3, sig_rho3, 1)
        plt.legend(lines, [r'$\rho_3(\theta)+$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_3$')
        #plt.savefig('w1_rho3.png')
        plt.savefig('w1_rho3.pdf')

        # Plots for worst rho2 exposure:
        # Find worst exposure based on rho2 at theta = 10 arcmin
        i = numpy.argmax(numpy.abs(exp_rho2p[:nexp,k10arcmin]), axis=0)
        print 'k10arcmin = ',k10arcmin
        print 'rho1[k] = ',exp_rho1p[:nexp,k10arcmin]
        print 'rho2[k] = ',exp_rho2p[:nexp,k10arcmin]
        print 'i = ',i
        print 'rho2[i] = ',exp_rho2p[i,:]
        meanr = numpy.exp(exp_meanlogr[i,:])
        rho1p = exp_rho1p[i,:]
        rho1m = exp_rho1m[i,:]
        rho2p = exp_rho2p[i,:]
        print 'rho2p = ',rho2p
        rho2m = exp_rho2m[i,:]
        rho3 = exp_rho3[i,:]
        sig_rho1p = numpy.sqrt(exp_var1[i,:])
        sig_rho1m = numpy.sqrt(exp_var1[i,:])
        sig_rho2p = numpy.sqrt(exp_var2[i,:])
        sig_rho2m = numpy.sqrt(exp_var2[i,:])
        sig_rho3 = numpy.sqrt(exp_var3[i,:])

        plt.clf()
        plt.title(r'$\rho_1$ for exposure with worst $\rho_2$ at 10 arcmin')
        lines = plot_rho(meanr, rho1p, sig_rho1p, 1, rho1m, sig_rho1m)
        plt.legend(lines, [r'$\rho_1(\theta)+$', r'$\rho_1(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_1$')
        #plt.savefig('w2_rho1.png')
        plt.savefig('w2_rho1.pdf')

        plt.clf()
        plt.title(r'$\rho_2$ for exposure with worst $\rho_2$ at 10 arcmin')
        lines = plot_rho(meanr, rho2p, sig_rho2p, 1, rho2m, sig_rho2m)
        plt.legend(lines, [r'$\rho_2(\theta)+$', r'$\rho_2(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_2$')
        #plt.savefig('w2_rho2.png')
        plt.savefig('w2_rho2.pdf')

        plt.clf()
        plt.title(r'$\rho_3$ for exposure with worst $\rho_2$ at 10 arcmin')
        lines = plot_rho(meanr, rho3, sig_rho3, 1)
        plt.legend(lines, [r'$\rho_3(\theta)+$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_3$')
        #plt.savefig('w2_rho3.png')
        plt.savefig('w2_rho3.pdf')

    if False:
        # Plots for desdm:
        print 'nexp = ',nexp
        sqrtn = numpy.sqrt(nexp)
        meanr = numpy.exp(numpy.mean(desdm_meanlogr[:nexp,:], axis=0))
        rho1p = numpy.mean(desdm_rho1p[:nexp,:], axis=0)
        rho1m = numpy.mean(desdm_rho1m[:nexp,:], axis=0)
        rho2p = numpy.mean(desdm_rho2p[:nexp,:], axis=0)
        rho2m = numpy.mean(desdm_rho2m[:nexp,:], axis=0)
        rho3 = numpy.mean(desdm_rho3[:nexp,:], axis=0)
        sig_rho1p = numpy.std(desdm_rho1p[:nexp,:], axis=0)
        sig_rho1m = numpy.std(desdm_rho1m[:nexp,:], axis=0)
        sig_rho2p = numpy.std(desdm_rho2p[:nexp,:], axis=0)
        sig_rho2m = numpy.std(desdm_rho2m[:nexp,:], axis=0)
        sig_rho3 = numpy.std(desdm_rho3[:nexp,:], axis=0)
        print 'meanr = ',meanr
        print 'rho1p = ',rho1p
        print 'sig_rho1p = ',sig_rho1p
        plt.rc('font', family='serif')

        plt.clf()
        plt.title(r'SPTE $\rho_1$ (i.e. $\langle de de \rangle$) for DESDM PSFEx solution')
        lines = plot_rho(meanr, rho1p, sig_rho1p, sqrtn, rho1m, sig_rho1m)
        plt.legend(lines, [r'$\rho_1(\theta)+$', r'$\rho_1(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_1$')
        #plt.savefig('desdm_rho1.png')
        plt.savefig('desdm_rho1.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_2$ (i.e. $\langle e de \rangle$) for DESDM PSFEx solution')
        lines = plot_rho(meanr, rho2p, sig_rho2p, sqrtn, rho2m, sig_rho2m)
        plt.legend(lines, [r'$\rho_2(\theta)+$', r'$\rho_2(\theta)-$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_2$')
        #plt.savefig('desdm_rho2.png')
        plt.savefig('desdm_rho2.pdf')

        plt.clf()
        plt.title(r'SPTE $\rho_3$ (i.e. $\langle ds ds \rangle$) for DESDM PSFEx solution')
        lines = plot_rho(meanr, rho3, sig_rho3, sqrtn)
        plt.legend(lines, [r'$\rho_3(\theta)$'] )
        plt.xlim( [0.5,100] )
        plt.ylabel(r'$\rho_3$')
        #plt.savefig('desdm_rho3.png')
        plt.savefig('desdm_rho3.pdf')

        plt.clf()
        pretty_rho1(meanr, rho1p, sig_rho1p, sqrtn)
        #plt.savefig('desdm_pretty_rho1.png')
        plt.savefig('desdm_pretty_rho1.pdf')

        plt.clf()
        pretty_rho2(meanr, rho2p, sig_rho2p, sqrtn)
        #plt.savefig('desdm_pretty_rho2.png')
        plt.savefig('desdm_pretty_rho2.pdf')

        plt.clf()
        pretty_rho3(meanr, rho3, sig_rho3, sqrtn)
        #plt.savefig('desdm_pretty_rho3.png')
        plt.savefig('desdm_pretty_rho3.pdf')

    if False:
        # Do some counts of how many exposures have high rho2
        count5 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 5.e-4).sum()
        count4 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 4.e-4).sum()
        count3 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 3.e-4).sum()
        count2 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 2.e-4).sum()
        count1 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 1.e-4).sum()
        count05 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 5.e-5).sum()
        count03 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 3.e-5).sum()
        count02 = (numpy.abs(exp_rho2p[:nexp,k10arcmin]) > 2.e-5).sum()

        print 'Exposure outliers:'
        print 'N with |rho2| > 5e-4 = ',count5
        print 'N with |rho2| > 4e-4 = ',count4
        print 'N with |rho2| > 3e-4 = ',count3
        print 'N with |rho2| > 2e-4 = ',count2
        print 'N with |rho2| > 1e-4 = ',count1
        print 'N with |rho2| > 5e-5 = ',count05
        print 'N with |rho2| > 3e-5 = ',count03
        print 'N with |rho2| > 2e-5 = ',count02

        count100 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 1.e-2).sum()
        count50 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 5.e-3).sum()
        count30 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 3.e-3).sum()
        count20 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 2.e-3).sum()
        count10 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 1.e-3).sum()
        count5 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 5.e-4).sum()
        count4 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 4.e-4).sum()
        count3 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 3.e-4).sum()
        count2 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 2.e-4).sum()
        count1 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 1.e-4).sum()
        count05 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 5.e-5).sum()
        count03 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 3.e-5).sum()
        count02 = (numpy.abs(ccd_rho2p[:nexp,k10arcmin]) > 2.e-5).sum()

        print 'CCD outliers:'
        print 'N with |rho2| > 1e-2 = ',count100
        print 'N with |rho2| > 5e-3 = ',count50
        print 'N with |rho2| > 3e-3 = ',count30
        print 'N with |rho2| > 2e-3 = ',count20
        print 'N with |rho2| > 1e-3 = ',count10
        print 'N with |rho2| > 5e-4 = ',count5
        print 'N with |rho2| > 4e-4 = ',count4
        print 'N with |rho2| > 3e-4 = ',count3
        print 'N with |rho2| > 2e-4 = ',count2
        print 'N with |rho2| > 1e-4 = ',count1
        print 'N with |rho2| > 5e-5 = ',count05
        print 'N with |rho2| > 3e-5 = ',count03
        print 'N with |rho2| > 2e-5 = ',count02
  
    if True:

        print 'Plot overall rho stats'

        #for key in ['griz', 'riz', 'g', 'r', 'i', 'z']:
        for key in ['riz']:
            stat_file = os.path.join(work, "rho_" + key + ".json")

            # Read the json file 
            with open(stat_file,'r') as f:
                stats = json.load(f)

            ( meanlogr,
              rho1p,
              rho1p_im,
              rho1m,
              rho1m_im,
              var1,
              rho2p,
              rho2p_im,
              rho2m,
              rho2m_im,
              var2,
              rho3p,
              rho3p_im,
              rho3m,
              rho3m_im,
              var3,
              rho4p,
              rho4p_im,
              rho4m,
              rho4m_im,
              var4,
              rho5p,
              rho5p_im,
              rho5m,
              rho5m_im,
              var5,
            ) = stats[-1]

            meanr = numpy.exp(meanlogr)
            rho1p = numpy.array(rho1p)
            rho1m = numpy.array(rho1m)
            rho2p = numpy.array(rho2p)
            rho2m = numpy.array(rho2m)
            rho3p = numpy.array(rho3p)
            rho3m = numpy.array(rho3m)
            rho4p = numpy.array(rho4p)
            rho4m = numpy.array(rho4m)
            rho5p = numpy.array(rho5p)
            rho5m = numpy.array(rho5m)
            sig_rho1 = numpy.sqrt(var1)
            sig_rho2 = numpy.sqrt(var2)
            sig_rho3 = numpy.sqrt(var3)
            sig_rho4 = numpy.sqrt(var4)
            sig_rho5 = numpy.sqrt(var5)
            sqrtn = 1

            cols = numpy.array((meanr, rho1p, sig_rho1, rho2p, sig_rho2,
                                rho1m, sig_rho1, rho2m, sig_rho2)).T
            numpy.savetxt('rho.dat', cols, fmt='%.6e',
                          header='meanr  rho1  sig_rho1  rho2  sig_rho2  rho1_xim  sig_rho1_xim  rho2_xim  sig_rho2_xim')
 
            plt.clf()
            pretty_rho1(meanr, rho1p, sig_rho1, sqrtn, rho3p, rho4p)
            #plt.savefig('rho1_' + key + '.png')
            plt.savefig('rho1_' + key + '.pdf')

            plt.clf()
            pretty_rho2(meanr, rho2p, sig_rho2, sqrtn, rho5p)
            #plt.savefig('rho2_' + key + '.png')
            plt.savefig('rho2_' + key + '.pdf')

            #plt.clf()
            #pretty_rho3(meanr, rho3, sig_rho3, sqrtn)
            #plt.savefig('rho3_' + key + '.png')
            #plt.savefig('rho3_' + key + '.pdf')
 
            if False:
                plt.clf()
                plt.title(r'$\rho_3$ (i.e. $\langle (e dT/T)* (e dT/T) \rangle$)')
                lines = plot_rho(meanr, rho3p, sig_rho3, sqrtn, rho3m, sig_rho3)
                plt.legend(lines, [r'$\rho_3(\theta)+$', r'$\rho_3(\theta)-$'] )
                plt.xlim( [0.5,300] )
                plt.ylabel(r'$\rho_3$')
                plt.tight_layout()
                plt.savefig('rho3.pdf')

                plt.clf()
                plt.title(r'$\rho_4$ (i.e. $\langle de* (e dT/T) \rangle$)')
                lines = plot_rho(meanr, rho4p, sig_rho4, sqrtn, rho4m, sig_rho4)
                plt.legend(lines, [r'$\rho_4(\theta)+$', r'$\rho_4(\theta)-$'] )
                plt.xlim( [0.5,300] )
                plt.ylabel(r'$\rho_4$')
                plt.tight_layout()
                plt.savefig('rho4.pdf')

                plt.clf()
                plt.title(r'$\rho_5$ (i.e. $\langle e* (e dT/T) \rangle$)')
                lines = plot_rho(meanr, rho5p, sig_rho5, sqrtn, rho5m, sig_rho5)
                plt.legend(lines, [r'$\rho_5(\theta)+$', r'$\rho_5(\theta)-$'] )
                plt.xlim( [0.5,300] )
                plt.ylabel(r'$\rho_5$')
                plt.tight_layout()
                plt.savefig('rho5.pdf')
 
if __name__ == "__main__":
    main()
