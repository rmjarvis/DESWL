import os, glob
import fitsio

meds_dir = '/astro/u/mjarvis/DES/meds'
meds_tag = 'y3v02a'
dirs = [ os.path.join(meds_dir, meds_tag) ]

outfile = meds_tag
bands = ['g','r','i','z','Y']
#bands = ['r']

all_exps = set()

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
                    # file name looks like /somedir/D$exp_$band_c$ccd_r$run_immasked.fits
                    file_name = file.split('/')[-1]
                    tokens = file_name.split('_')
                    exp = tokens[0][1:]
                    band = tokens[1]
                    ccd = tokens[2][1:]
                    run = tokens[3][1:]
                    print(exp,band,ccd,run)
                    exps.add(str(int(exp)))
                    all_exps.add(str(int(exp)))

    with open(out_b,'w') as out:
        for exp in sorted(list(exps)):
            out.write(exp + '\n')

out_all = outfile + '_' + ''.join(bands)
with open(out_all,'w') as out:
    for exp in sorted(list(all_exps)):
        out.write(exp + '\n')
