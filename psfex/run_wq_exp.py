#!/usr/bin/env python
# Program to loop over exposures and run a given command
#
# the njob argument will batch the exposures into njobs/per wq job
# 
# the file argument  expects a file containing a list of run,expnum such as:
#     
#    20130711064512_20121202 DECam_00157539
# 
# the command argument is the name of exectuation that is expected 
# to take a --exps and --runs to specify input directories
# and --output for output
# 
# I only have two scripts that do this currently run_psfex.py and run_findstars.py
# To run them you can do something like
#    ./run_wq_exp.py --njob 1 --file test_exp --output test2 --cmd="./run_psfex.py --exp_match \"*_[0-9][0-9].fits*\" --use_findstars 1 --mag_cut 3.0"

# 

import argparse,os,re
import time
def file_len(file):
    lines = 0
    for line in file:
        lines += 1
    return lines
                    
top_txt="""
# these are the commands to be run.  if you only have a 
# single command, you can use a single line such as 
# command: ./script

# The config files are in the production directory
# Old version uses:
#     source /astro/u/mjarvis/.bashrc
#     {cmd}
# I switched to bash -l -c '{cmd}' because my shell seems to be set up wrong.
command: |
    cd /astro/u/mjarvis/rmjarvis/DESWL/psfex
    bash -l -c '{cmd}'

# show this name in job listings instead of the command
job_name: {name}

# this is the type of node/host selection. bynode means select entire
# nodes.
mode: bycore

# Since the mode is bynode, this means 5 full nodes
# N: 5

# Select from this group(s)
group: [new, new2, new3]

# Do not select from this set of groups
notgroup: [slow,crappy]

# require at least this many cores
# min_cores: 8

# used by MPI jobs
# hostfile: auto

# If we have 5 full nodes of 12 cores each,
# there is 60 cores in total. Threads:4 ensures each
# host is listed 3 times. So the command above will
# run 15 MPI nodes of 4 threads each
# threads: 4

"""


parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--njob', default=10, type=int,
                    help='exposures per job')
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
files=open(args.file)
nline=file_len(files)
files.close()

files=open(args.file)

# list of runs
runs=''
exps=''

ijob=0
i=0
iline=0

for line in files:

    iline+=1
    vals=line.split()
    run=vals[0]
    exp=vals[1]

    runs+=run+' '
    exps+=exp+' '
    
    if (i+1)%args.njob==0 or iline==nline:
        
        cmd=args.cmd+' --runs %s --exps %s --output %s'%(runs,exps,args.output)

        job_submit = top_txt.format(runs=runs, exps=exps, output=args.output,name=str(ijob),cmd=cmd)
        if args.debug:
            job_submit += "priority: high\n"

        submit_file='%s/submit_%d'%(args.submit_dir,ijob)
        submit_out='%s/submit_%d.out'%(args.submit_dir,ijob)
        file=open(submit_file,'w')
        file.write(job_submit)
        file.close()
        time.sleep(0.1)
        runs=''
        exps=''
        os.system('nohup wq sub -b %s >& %s'%(submit_file,submit_out))
        time.sleep(0.1)
        ijob+=1

    i+=1
