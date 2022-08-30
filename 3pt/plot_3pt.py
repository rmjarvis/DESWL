import treecorr
import numpy as np
import matplotlib.pyplot as plt

rmin = 1
rmax = 25
nr = 10

narrow = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
              min_u=0.0, max_u=1, nubins=20,
              min_v=0.0, max_v=0.1, nvbins=1, verbose=2)

wide = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
            min_u=0.9, max_u=1, nubins=1,
            min_v=0.0, max_v=0.8, nvbins=20, verbose=2)

wider = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
             min_u=0.9, max_u=1, nubins=1,
             min_v=0.8, max_v=0.95, nvbins=20, verbose=2)

widest = dict(min_sep=rmin, max_sep=rmax, sep_units='arcmin', nbins=nr,
              min_u=0.9, max_u=1, nubins=1,
              min_v=0.95, max_v=1.0, nvbins=20, verbose=2)

ggg1 = treecorr.GGGCorrelation(narrow)
ggg1.read('narrow.hdf')
ggg2 = treecorr.GGGCorrelation(wide)
ggg2.read('wide.hdf')
ggg3 = treecorr.GGGCorrelation(wider)
ggg3.read('wider.hdf')
ggg4 = treecorr.GGGCorrelation(widest)
ggg4.read('widest.hdf')

all_g_ttt = []
all_sig_ttt = []
all_meanr = []
all_phi = []

for ggg in [ggg1, ggg2, ggg3, ggg4]:

    g_ttt = -0.25 * (ggg.gam0 + ggg.gam1 + ggg.gam2 + ggg.gam3).real
    var_ttt = 0.25**2 * (ggg.vargam0 + ggg.vargam1 + ggg.vargam2 + ggg.vargam3)

    _nr, nu, nv = g_ttt.shape
    print(nr,nu,nv)
    assert _nr == nr
    assert nv % 2 == 0
    nv //= 2
    assert nu == 1 or nv == 1

    d1 = ggg.meand1
    d2 = ggg.meand2
    d3 = ggg.meand3
    if nu == 1:
        # if nu==1, then u=1, so d2 = d3, and phi is between d2 and d3
        phi = np.arccos( (d2**2 + d3**2 - d1**2) / (2*d2*d3) )
        meanr = np.array([np.mean([d2[ir], d3[ir]]) for ir in range(nr)])
    else:
        # if nv==1, then v=0, so d1 = d2, and phi is between d1 and d2
        phi = np.arccos( (d1**2 + d2**2 - d3**2) / (2*d1*d2) )
        meanr = np.array([np.mean([d1[ir], d2[ir]]) for ir in range(nr)])
    phi *= 180/np.pi

    # We don't care about v>0 vs v<0, so combine them.
    phi = (phi[:,:,nv-1::-1] + phi[:,:,nv:]) / 2
    g_ttt = (g_ttt[:,:,nv-1::-1] + g_ttt[:,:,nv:]) / 2
    var_ttt = (var_ttt[:,:,nv-1::-1] + var_ttt[:,:,nv:]) / 4
    sig_ttt = var_ttt**0.5

    print('shapes:')
    print('phi: ',phi.shape)
    print('g_ttt: ',g_ttt.shape)
    print('sig_ttt: ',sig_ttt.shape)
    print('meanr: ',meanr.shape)

    print('meanr =  ',meanr)

    if nu == 1:
        phi = phi[:,0,:]
        g_ttt = g_ttt[:,0,:]
        sig_ttt = sig_ttt[:,0,:]
    else:
        phi = phi[:,:,0]
        g_ttt = g_ttt[:,:,0]
        sig_ttt = sig_ttt[:,:,0]

    print('shapes ->')
    print('phi: ',phi.shape)
    print('g_ttt: ',g_ttt.shape)
    print('sig_ttt: ',sig_ttt.shape)

    all_phi.append(phi)
    all_g_ttt.append(g_ttt)
    all_sig_ttt.append(sig_ttt)
    all_meanr.append(meanr)

phi = np.concatenate(all_phi, axis=1)
g_ttt = np.concatenate(all_g_ttt, axis=1)
sig_ttt = np.concatenate(all_sig_ttt, axis=1)
meanr = np.concatenate(all_meanr, axis=0)

fig, ax = plt.subplots()

lines = []
for ir in range(nr):
    print('ir = ',ir)
    print('meanr = ',meanr[ir])
    print('phi = ',phi[ir])
    print('g = ',g_ttt[ir])
    print('sig = ',sig_ttt[ir])
    print()

    line = ax.errorbar(phi[ir], g_ttt[ir], sig_ttt[ir])
    lines.append((line, 'd2 ~= %.1f'%meanr[ir]))

ax.legend(*(list(zip(*lines))))
ax.set_xlabel(r'Opening angle $\phi$ [deg]')
ax.set_ylabel(r'$\gamma_{\rm ttt}$ isoceles')

fig.set_tight_layout(True)
plt.savefig('isoc.pdf')
plt.savefig('isoc.png')
