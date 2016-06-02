import os, glob
import fitsio

meds_dir = '/astro/u/mjarvis/DES/meds'
meds_tag = 'y1a1-spt-002'
dirs = [ os.path.join(meds_dir, meds_tag) ]

outfile = 'y1all-002'
bands = ['g','r','i','z']
#dir = '/astro/u/mjarvis/DES/meds/tb-y1a1-v01'
#outfile = 'tb-y1a1'

for b in bands:
    out_b = outfile + '_' + b
    print out_b
    exps = set()
    for dir in dirs:
        pat = os.path.join(dir,'*','*-'+b+'-meds-stubby-%s.fits'%meds_tag)
        print pat
        for srclist in glob.glob(pat):
            #if 'DES2305-0124' not in srclist:
                #continue
            print srclist
            with fitsio.FITS(srclist) as f:
                files = f['image_info']['image_path'][:]
                for file in files:
                    if 'coadd' in file: continue
                    # file name looks like ${DESDATA}/OPS/red/$run/red/$exp/filename
                    tokens = file.split('/')
                    run = tokens[-4]
                    exp = tokens[-2]
                    #print run,exp
                    exps.add((run,exp))

    with open(out_b,'w') as out:
        for run,exp in sorted(list(exps)):
            out.write(run + ' ' + exp + '\n')
