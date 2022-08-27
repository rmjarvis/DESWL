import treecorr
import numpy as np
import matplotlib.pyplot as plt

params_3pt = dict(min_sep=1, max_sep=5, sep_units='arcmin', nbins=5,
                min_u=0.9, max_u=1, nubins=1,
                min_v=0.7, max_v=1.0, nvbins=10, verbose=2)
ggg = treecorr.GGGCorrelation(params_3pt, var_method='bootstrap')
ggg.read('ggg.hdf')

u = ggg.meanu
v = ggg.meanv
d1 = ggg.meand1
d2 = ggg.meand2
d3 = ggg.meand3

g_ttt = -0.25 * (ggg.gam0 + ggg.gam1 + ggg.gam2 + ggg.gam3).real
var_ttt = 0.25**2 * (ggg.vargam0 + ggg.vargam1 + ggg.vargam2 + ggg.vargam3)

phi = np.arccos( (d2**2 + d3**2 - d1**2) / (2*d2*d3) )
phi *= 180/np.pi

nr, nu, nv = phi.shape

print(nr,nu,nv)
assert nu == 1
assert nv % 2 == 0
nv //= 2

fig, ax = plt.subplots()

lines = []
for ir in range(nr):
    meanr = np.mean([d2[ir], d3[ir]])
    _phi = phi[ir][0]

    # We don't care about v>0 vs v<0, so combine them.
    _phi = (phi[ir][0][nv-1::-1] + phi[ir][0][nv:]) / 2
    _g = (g_ttt[ir][0][nv-1::-1] + g_ttt[ir][0][nv:]) / 2
    _var = (var_ttt[ir][0][nv-1::-1] + var_ttt[ir][0][nv:]) / 4
    _sig = _var**0.5

    print('ir = ',ir)
    print('meanr = ',meanr)
    print('phi = ',_phi)
    print('g = ',_g)
    print('sig = ',_sig)
    print()

    line = ax.errorbar(_phi, _g, _sig)
    lines.append((line, 'd2 ~= %.1f'%meanr))

ax.legend(*(list(zip(*lines))))
ax.set_xlabel(r'Opening angle $\phi$ [deg]')
ax.set_ylabel(r'$\gamma_{\rm ttt}$ isoceles')

fig.set_tight_layout(True)
plt.savefig('isoc.pdf')
plt.savefig('isoc.png')
