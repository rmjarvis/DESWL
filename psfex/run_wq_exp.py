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
# and --output for output
# 
# I only have two scripts that do this currently run_psfex.py and run_findstars.py
# To run them you can do something like
#    ./run_wq_exp.py --njobs 10 --file test_exp --output test2 --cmd="./run_psfex.py --exp_match \"*_[0-9][0-9].fits*\" --use_findstars 1 --mag_cut 3.0"

# 
# See https://github.com/esheldon/wq for information about using wq.

import argparse,os,re
import time
import numpy

# Note: I originally had the wrong shell set.  I ran bash in my .login
# rather than setting it correctly by mailing RT-RACF-UserAccounts@bnl.gov.
# If you have the same problem, you can have the command look like:
#    cd /astro/u/mjarvis/rmjarvis/DESWL/psfex
#    bash -l -c '{cmd}'
# instead of 
#    cd /astro/u/mjarvis/rmjarvis/DESWL/psfex
#    source /astro/u/mjarvis/.bashrc
#    {cmd}

top_txt="""
command: |
    cd /astro/u/mjarvis/rmjarvis/DESWL/psfex
    source /astro/u/mjarvis/.bashrc
    {cmd}

job_name: {name}

# this is the type of node/host selection. bynode means select entire
# nodes.
mode: bycore

# Select from this group(s)
group: [new, new2, new3]
"""


parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--njobs', default=100, type=int,
                    help='How many jobs to run')
parser.add_argument('--file', default='',
                    help='list of run/exposures')
parser.add_argument('--submit_dir',default='submit',
                    help='where to put submit files')
parser.add_argument('--output', default='/astro/u/mjarvis/work',
                    help='where to put output files')
parser.add_argument('--cmd', default='./run_findstars.py --condor 0',
                    help='command to run on the exposures')
parser.add_argument('--debug', default=False, action='store_const', const=True,
                    help='Set priority to high for debugging run')


args = parser.parse_args()


if not os.path.isdir(args.submit_dir): os.makedirs(args.submit_dir)

# open it once to count the number of lines
# I'm sure there is a better way to do this
print 'Read file ',args.file
with open(args.file) as fin:
    data = [ line.split() for line in fin ]

# Shuffle the order so we don't have all the LMC exposures in the same job.
print 'first 3 lines of input file are ',data[0:3]
numpy.random.shuffle(data)
print 'After shuffling, first 3 lines of input file are ',data[0:3]

runs, exps = zip(*data)

import math
n_per_job = int(math.ceil(float(len(runs)) / float(args.njobs)))
print 'njobs = ',args.njobs
print 'total n = ',len(runs)
print 'n_per_job = ',n_per_job

submit_list = []

for job in range(args.njobs):
    start = job * n_per_job
    end = start + n_per_job
    if end > len(runs): 
        end = len(runs)
    if end == start:
        continue

    # Make single string with a list of the runs and exps for this job:
    s_runs = " ".join(runs[start:end])
    s_exps = " ".join(exps[start:end])

    cmd=args.cmd+' --runs %s --exps %s --output %s'%(s_runs,s_exps,args.output)

    job_submit = top_txt.format(runs=runs, exps=exps, output=args.output,name=str(job),cmd=cmd)
    if args.debug:
        job_submit += "priority: high\n"

    submit_file='%s/submit_%d'%(args.submit_dir,job)

    with open(submit_file,'w') as fout:
        fout.write(job_submit)
    submit_list.append(submit_file)

time.sleep(0.1)
s_sub = " ".join(submit_list)
cmd = 'nohup wq sub -b %s >& %s/wq_sub.out'%(s_sub,args.submit_dir)
print cmd
print 'Note: This will take %d seconds to run, since wq waits'%len(submit_list)
print '      1 second between each job submission.'
os.system(cmd)

