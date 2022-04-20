import piff
import numpy as np
import yaml
import galsim

info_name = 'meds/des-pizza-slices-y6-v12/pizza_cutter_info/DES0147-6122_g_pizza_cutter_info.yaml'
with open(info_name) as f:
    info = yaml.safe_load(f)

print(len(info['src_info']))
print('g: ',np.sum([src['band']=='g' for src in info['src_info']]))
print('r: ',np.sum([src['band']=='r' for src in info['src_info']]))
print('i: ',np.sum([src['band']=='i' for src in info['src_info']]))
print('z: ',np.sum([src['band']=='z' for src in info['src_info']]))

nsrc = None
stepx = stepy = 128
stepcolor = 0.2
use_bands = ['g']

data = []

colors = np.arange(0, 3.01, stepcolor)
xs = np.arange(stepx/2, 2048, stepx)
ys = np.arange(stepy/2, 4096, stepy)

ncolor = len(colors)
nx = len(xs)
ny = len(ys)

def output(n,x,y,color,im):
    fname = 'HSM_fail_{}_{}_{}_{:0.1f}.fits'.format(n,int(x),int(y),color)
    im.write(fname)

for n, src in enumerate(info['src_info'][:nsrc]):
    band = src['band']
    if band not in use_bands:
        continue
    image_path = src['image_path']
    piff_path = src['piff_path']
    piff_info = src['piff_info']  # ccdnum, expnum, nstar, desdm_flags, fwhm_cen,
                                  # exp_star_t_mean, exp_star_t_std, star_t_mean, star_t_std
    ccdnum = piff_info['ccdnum']
    
    print(n,band,piff_path)
    psf = piff.read(piff_path)

    allT = np.zeros((ncolor, nx, ny), dtype=float)
    allg1 = np.zeros((ncolor, nx, ny), dtype=float)
    allg2 = np.zeros((ncolor, nx, ny), dtype=float)

    for icolor, color in enumerate(colors):
        for ix, x in enumerate(xs):
            for iy, y in enumerate(ys):
                im = psf.draw(x=x, y=y, GI_COLOR=color, chipnum=ccdnum)
                # Do this now so the image is trivially writable.
                im.wcs = im.wcs.jacobian(image_pos=galsim.PositionD(x,y))
                try:
                    shape_data = im.FindAdaptiveMom()
                except Exception as e:
                    print('Caught ',e)
                    output(n,x,y,color,im)
                    allT[icolor,ix,iy] = -1
                    allg1[icolor,ix,iy] = -1
                    allg2[icolor,ix,iy] = -1
                    continue
                else:
                    if shape_data.moments_status != 0:
                        print('moments_status = ',shape_data.moments_status)
                        output(n,x,y,color,im)
                        allT[icolor,ix,iy] = -1
                        allg1[icolor,ix,iy] = -1
                        allg2[icolor,ix,iy] = -1
                        continue
                    else:
                        e1 = shape_data.observed_shape.e1
                        e2 = shape_data.observed_shape.e2
                        s = shape_data.moments_sigma

                        # Apply WCS jacobian
                        # TODO: This should really be a GalSim function...
                        M = np.matrix( [[ 1 + e1, e2 ], [ e2, 1 - e1 ]] ) * s*s
                        jac = im.wcs
                        J = jac.getMatrix()
                        M = J * M * J.T

                        e1 = (M[0,0] - M[1,1]) / (M[0,0] + M[1,1])
                        e2 = (2.*M[0,1]) / (M[0,0] + M[1,1])
                        T = M[0,0] + M[1,1]

                        shear = galsim.Shear(e1=e1, e2=e2)
                        g1 = shear.g1
                        g2 = shear.g2

                allT[icolor,ix,iy] = T
                allg1[icolor,ix,iy] = g1
                allg2[icolor,ix,iy] = g2

    medianT = np.median(allT)
    maxT = np.max(allT)
    minT = np.min(allT[allT>=0])
    nfail = np.sum(allT < 0)
    bad = np.where((allT>=0) & ((allT < medianT - 0.15) | (allT > medianT + 0.15)))
    print(n,medianT,maxT-minT,nfail,len(bad[0]))
    if len(bad[0]) > 0:
        print()
        print('Found some PSFs that exceed 0.15 in |T - T_median|')
        print('  median = ',medianT)
        print('  bad = ',allT[bad])
        print('  at ',bad)
        print()
    data.append( (n,medianT,minT,maxT,nfail,bad,allT,allg1,allg2) )

np.savez('piff_stats.npz', data=data)

for icolor, color in enumerate(colors):
    T_list = np.array([ d[6][icolor,:,:] for d in data ]).ravel()
    print(color, len(T_list))
    mask = T_list >= 0
    print('T: ',np.mean(T_list[mask]), np.min(T_list[mask]), np.max(T_list), np.sum(~mask) )
    print()

