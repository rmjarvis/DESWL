import numpy
import astropy.io.fits as pyfits
from scipy.spatial import KDTree

dir = 'e2e_v4'
truth = pyfits.open(dir+'/end2end-truth.fits')[1].data
new = pyfits.open(dir+'/DES0436-5748_r_cat.fits')[1].data

print 'truth has %d columns, %d entries'%(len(truth.columns), len(truth))
print 'new has %d columns, %d entries'%(len(new.columns), len(new))

mean_dec = truth['dec'].mean()
cosdec = numpy.cos(mean_dec * numpy.pi/180.)

true_pos = numpy.empty( (len(truth), 2) )
true_pos[:,0] = truth['ra'] * cosdec
true_pos[:,1] = truth['dec']

tree = KDTree(true_pos)

new_pos = numpy.empty( (len(new), 2) )
new_pos[:,0] = new['ALPHAWIN_J2000'] * cosdec
new_pos[:,1] = new['DELTAWIN_J2000']

dist, index = tree.query(new_pos)

ok_dist = (dist < 0.1 / 3600.)
ok_sex = (new['FLAGS'] == 0)
ok_truth = (truth['flags'][index] == 0)

ok = ok_dist & ok_sex & ok_truth

match_file = dir+'/match.fits'
columns = [ pyfits.Column(name='index', format='J', array=index),
            pyfits.Column(name='ok', format='J', array=ok)
            ]
coldefs = pyfits.ColDefs(columns)
table = pyfits.new_table(coldefs)
table.writeto(match_file, clobber=True)

