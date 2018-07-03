#!/usr/bin/env python
# Program to loop over exposures and run a given command
#
# The njobs argument will bundle the exposures into batches to make a total
# of njobs jobs to submit to wq.
# 
# The file argument expects a file containing a list of run,expnum such as:
#     
#    20130711064512_20121202 DECam_00157539
# 
# The command argument is the name of exectuation that is expected 
# to take a --exps and --runs to specify input directories
# 
# I only have two scripts that do this currently run_psfex.py and run_findstars.py
# To run them you can do something like
#    ./run_wq_exp.py --njobs 10 --file test_exp --cmd="./run_psfex.py --exp_match \"*_[0-9][0-9].fits.fz\" --use_findstars 1 --mag_cut 3.0"

# 
# See https://github.com/esheldon/wq for information about using wq.

from __future__ import print_function

import argparse,os,re
import time
import numpy
import datetime
import numpy as np

# Note: I originally had the wrong shell set.  I ran bash in my .login
# rather than setting it correctly by mailing RT-RACF-UserAccounts@bnl.gov.
# If you have the same problem, you can have the command look like:
#    cd /astro/u/mjarvis/rmjarvis/DESWL/psf
#    bash -l -c '{cmd}'
# instead of 
#    cd /astro/u/mjarvis/rmjarvis/DESWL/psf
#    source /astro/u/mjarvis/.bashrc
#    {cmd}

top_txt="""
command: |
    cd /astro/u/mjarvis/rmjarvis/DESWL/psf
    source /astro/u/mjarvis/.bashrc
    conda activate py2.7
    export OMP_NUM_THREADS={cores_per_job}
    export MKL_NUM_THREADS={cores_per_job}
    {cmd}

job_name: {name}

# this is the type of node/host selection. bynode means select entire
# nodes.  bycore1 means all N cores from the same node.
mode: bycore1
N: {cores_per_job}

# Select from this group(s)
# I've had trouble with OS operations from a31, a32, so avoid them.
group: [gen7, gen6, gen4, a22, a24, a25, a26, a27, a28, a29, a30, gen5]
"""


parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--njobs', default=100, type=int,
                    help='How many jobs to run')
parser.add_argument('--cores_per_job', default=1, type=int,
                    help='How many cores to use per job')
parser.add_argument('--file', default='',
                    help='list of run/exposures')
parser.add_argument('--submit_dir',default=None,
                    help='where to put submit files')
parser.add_argument('--cmd', default='./run_findstars.py --condor 0',
                    help='command to run on the exposures')
parser.add_argument('--debug', default=False, action='store_const', const=True,
                    help='Set priority to high for debugging run')
parser.add_argument('--tag', default=None,
                    help='Some kind of tag for the auto-generated submit dir')


args = parser.parse_args()

if args.submit_dir is None:
    date = datetime.date.today()
    if args.tag is None:
        tag = args.file
    else:
        tag = args.tag
    submit_dir = '~/work/submit_%d%02d%02d_%s'%(date.year,date.month,date.day,tag)
else:
    submit_dir = args.submit_dir

submit_dir = os.path.expanduser(submit_dir)
print('submit_dir = ',submit_dir)
if not os.path.isdir(submit_dir): os.makedirs(submit_dir)

# Read in the runs, exps from the input file
print('Read file ',args.file)
with open(args.file) as fin:
    lines = [ line for line in fin if line[0] != '#' ]
nexps = len(lines)
print('Exposures file has %d lines'%nexps)
lines = np.unique(lines)
if nexps != len(lines):
    nexps = len(lines)
    print('Warning: the exposure list has duplicate entries!')
    print('Only %d unique exposures'%nexps)

if args.njobs != 1:
    # Shuffle the order so we don't have all the LMC exposures in the same job.
    print('first 3 lines of input file are ',lines[0:3])
    numpy.random.shuffle(lines)
    print('After shuffling, first 3 lines of input file are ',lines[0:3])

data = [ line.split() for line in lines ]

if len(data[0]) == 2:
    runs, exps = zip(*data)
else:
    runs = None
    exps = list(zip(*data))[0]
#print('runs = ',runs)
#print('exps = ',exps)

if args.njobs > nexps:
    args.njobs = nexps

import math
n_per_job = int(math.ceil(float(nexps) / float(args.njobs)))
print('njobs = ',args.njobs)
print('total n = ',nexps)
print('n_per_job = ',n_per_job)

submit_list = []

end = 0
for job in range(args.njobs):
    start = end
    end = int(math.ceil(float(nexps) * (job+1) / float(args.njobs)))
    if end > nexps: 
        end = nexps
    if end == start:
        continue

    if args.njobs == 1:
        cmd=args.cmd+' --file %s'%(args.file)
    else:
        # Make single string with a list of the runs and exps for this job:
        if runs is not None:
            s_runs = " ".join(runs[start:end])
            cmd=args.cmd+' --runs %s'%(s_runs)
        s_exps = " ".join(exps[start:end])
        cmd=args.cmd+' --exps %s'%(s_exps)

    job_name = args.file + '_' + str(job)
    job_submit = top_txt.format(name=job_name, cores_per_job=args.cores_per_job,
                                cmd=cmd)
    if args.debug:
        job_submit += "priority: high\n"

    submit_file='%s/submit_%s'%(submit_dir,job_name)

    with open(submit_file,'w') as fout:
        fout.write(job_submit)
    submit_list.append(submit_file)

time.sleep(0.1)
s_sub = " ".join(submit_list)
cmd = 'nohup wq sub -b %s >& %s/wq_sub_%s.out'%(s_sub,submit_dir,args.file)
#print(cmd)
print('Note: This will take %d seconds to run, since wq waits 1 second'%len(submit_list))
print('      between each job submission.')
os.system(cmd)

