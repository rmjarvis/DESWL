import piff
import numpy as np
import yaml
import galsim
import matplotlib.pyplot as plt

info_name = 'meds/des-pizza-slices-y6-v12/pizza_cutter_info/DES0147-6122_g_pizza_cutter_info.yaml'
with open(info_name) as f:
    info = yaml.safe_load(f)

# This is the only one that tripped Matt's 0.15 test.
n = 88

stepx = stepy = 128
stepcolor = 0.2

colors = np.arange(0, 3.01, stepcolor)
xs = np.arange(stepx/2, 2048, stepx)
ys = np.arange(stepy/2, 4096, stepy)

ncolor = len(colors)
nx = len(xs)
ny = len(ys)

src = info['src_info'][n]

band = src['band']
image_path = src['image_path']
piff_path = src['piff_path']
piff_info = src['piff_info']  # ccdnum, expnum, nstar, desdm_flags, fwhm_cen,
                              # exp_star_t_mean, exp_star_t_std, star_t_mean, star_t_std
ccdnum = piff_info['ccdnum']

print(n,band,piff_path)
psf = piff.read(piff_path)

print(piff_info)
print('chisq = ',psf.chisq)
print('dof = ',psf.dof)
print('nremoved = ',psf.nremoved)
print('nused = ',len(psf.stars))
print('colors = ',[s['GI_COLOR'] for s in psf.stars])

allT = np.zeros((ncolor, nx, ny), dtype=float)
allg1 = np.zeros((ncolor, nx, ny), dtype=float)
allg2 = np.zeros((ncolor, nx, ny), dtype=float)

for icolor, color in enumerate(colors):
    print('color = ',color)
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            im = psf.draw(x=x, y=y, GI_COLOR=color, chipnum=ccdnum)
            # Do this now so the image is trivially writable.
            im.wcs = im.wcs.jacobian(image_pos=galsim.PositionD(x,y))
            try:
                shape_data = im.FindAdaptiveMom()
            except Exception as e:
                print('Caught ',e)
                allT[icolor,ix,iy] = -1
                allg1[icolor,ix,iy] = -1
                allg2[icolor,ix,iy] = -1
                continue
            else:
                if shape_data.moments_status != 0:
                    print('moments_status = ',shape_data.moments_status)
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
print(n,medianT,minT,maxT,nfail)

nx = int(np.ceil(np.sqrt(ncolor)))
ny = int(np.ceil(ncolor / nx))
fig, ax = plt.subplots(nx,ny, figsize=(15,8))

x,y = np.meshgrid(xs,ys)

icolor = 0
for iy in range(ny):
    for ix in range(nx):
        if icolor >= ncolor:
            break
        c = ax[iy,ix].pcolormesh(x, y, allT[icolor].T, vmin=minT, vmax=maxT)
        ax[iy,ix].set_title('g-i = %.1f'%(colors[icolor]))
        if iy == ny-1:
            ax[iy,ix].set_xlabel('x')
        else:
            ax[iy,ix].axes.xaxis.set_ticklabels([])
        if ix == 0:
            ax[iy,ix].set_ylabel('y')
        else:
            ax[iy,ix].axes.yaxis.set_ticklabels([])
        icolor += 1

fig.suptitle('D00378642_g_c32_r5738p01')
fig.colorbar(c, ax=ax.ravel().tolist())

fig.savefig('D00378642_g_c32_r5738p01.png', bbox_inches='tight')
