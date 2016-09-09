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

import argparse,os,re
import time
import numpy
import datetime

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
# nodes.  bycore1 means all N cores from the same node.
mode: bycore1
N: {cores_per_job}

# Select from this group(s)
# I've had trouble with OS operations from neww2, so avoid that one.
# astro0001 in particular, so possibly could try the other ones in new2.
# I didn't try them individually.
group: [new, new2, new3]
"""


parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--njobs', default=100, type=int,
                    help='How many jobs to run')
parser.add_argument('--cores_per_job', default=1, type=int,
                    help='How many cores to use per job')
parser.add_argument('--submit_dir',default=None,
                    help='where to put submit files')
parser.add_argument('--exp_match',
                    help='the exp_match argument with a %d that will be replaced by the job number')
parser.add_argument('--cmd', help='command to run on the exposures (required)')
parser.add_argument('--debug', default=False, action='store_const', const=True,
                    help='Set priority to high for debugging run')
parser.add_argument('--tag', default=None,
                    help='Some kind of tag for the auto-generated submit dir')


args = parser.parse_args()

if args.submit_dir is None:
    date = datetime.date.today()
    if args.tag is None:
        tag = 'alt'
    else:
        tag = args.tag
    submit_dir = '~/work/submit_%d%02d%02d_%s'%(date.year,date.month,date.day,tag)
else:
    submit_dir = args.submit_dir

submit_dir = os.path.expanduser(submit_dir)
print 'submit_dir = ',submit_dir
if not os.path.isdir(submit_dir): os.makedirs(submit_dir)

submit_list = []

end = 0
for job in range(args.njobs):

    cmd=args.cmd+' --exp_match=' + args.exp_match%job
    #print 'cmd = ',cmd

    job_name = args.tag + '_' + str(job)
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
cmd = 'nohup wq sub -b %s >& %s/wq_sub_%s.out'%(s_sub,submit_dir,args.tag)
print cmd
print 'Note: This will take %d seconds to run, since wq waits 1 second'%len(submit_list)
print '      between each job submission.'
os.system(cmd)

