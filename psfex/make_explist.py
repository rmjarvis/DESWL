import os, glob

dirs = [ 
        '/astro/u/mjarvis/DES/meds/spt-e-gold',
        '/astro/u/mjarvis/DES/meds/y1a1-alpha',
        '/astro/u/mjarvis/DES/meds/y1a1-beta',
        '/astro/u/mjarvis/DES/meds/y1a1-gamma',
        '/astro/u/mjarvis/DES/meds/y1a1-delta',
       ]

outfile = 'y1all'
bands = ['g','r','i','z']
#dir = '/astro/u/mjarvis/DES/meds/tb-y1a1-v01'
#outfile = 'tb-y1a1'

for b in bands:
    out_b = outfile + '_' + b
    print out_b
    exps = set()
    for dir in dirs:
        print os.path.join(dir,'*','*-'+b+'-meds-srclist*.dat')
        for srclist in glob.glob(os.path.join(dir,'*','*-'+b+'-meds-srclist*.dat')):
            print srclist
            with open(srclist,'r') as f:
                for line in f:
                    vals = line.split()
                    #print vals[2]
                    file = vals[2]
                    # file name looks like /astro/u/astrodat/data/DES/OPS/red/$run/red/$exp/filename
                    tokens = file.split('/')
                    run = tokens[-4]
                    exp = tokens[-2]
                    #print run,exp
                    exps.add((run,exp))

    with open(out_b,'w') as out:
        for run,exp in sorted(list(exps)):
            out.write(run + ' ' + exp + '\n')
