import os, glob
import fitsio

meds_dir = '/astro/u/mjarvis/DES/meds'
meds_tag = 'y3v02'
dirs = [ os.path.join(meds_dir, meds_tag) ]

outfile = meds_tag
#bands = ['g','r','i','z']
bands = ['r']

for b in bands:
    out_b = outfile + '_' + b
    print out_b
    exps = set()
    for dir in dirs:
        pat = os.path.join(dir,'*','*_'+b+'_meds-*.fits*')
        print pat
        for srclist in sorted(glob.glob(pat)):
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
