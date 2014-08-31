#!/usr/bin/env python
import argparse,os,re

parser = argparse.ArgumentParser(description='Run findstars on multiple exposures')
parser.add_argument('--runs',default='',nargs='+',
                    help='list of runs')

parser.add_argument('--exps',default='',nargs='+',
                    help='list of exposures')
parser.add_argument('--output',default='',
                    help='output directory of the products')

# This option has never worked well for me, don't use
parser.add_argument('--condor',default=0,type=int,
                    help='use condor')

parser.add_argument('--find',default=1,type=int,
                    help='output directory')
parser.add_argument('--wl_soft',default='/astro/u/rarmst/findstars/wl_trunk/src/',
                    help='software directory')
parser.add_argument('--input_dir',default='',
                    help='input directory of the data')

args = parser.parse_args()
hdfs='/astro/astronfs03/mapr/hadoop/hadoop-0.20.2/bin/hadoop fs -fs maprfs://astronfs03.rcf.bnl.gov:7222'

# For condor it may better to work locally and then cp back at the end of processing
local='./'
if args.condor:
      local=os.environ['_CONDOR_SCRATCH_DIR']+'/'


for run,exp in zip(args.runs,args.exps):


      if not args.condor:
            odir=args.output+'/'+exp+'/'
            logfile=odir+'/log.'+exp
      else:
            odir=local+'/'+exp+'/'
            logfile=odir+'/log.'+exp

      # create output path
      if not os.path.exists(odir):
            os.makedirs(odir)

      # create input dir if not given from des data dir
      if args.input_dir=='':
            if args.condor:
                  datadir='/data/esheldon/desdata/OPS'
            else:
                  datadir='/astro/u/astrodat/data/DES/OPS'
      
            input_dir='%s/red/%s/red/%s/'%(datadir,run,exp)
      else:
            input_dir=args.input_dir

      for ccd in range(1,63):
            # skip chip 61
            if ccd==61:continue
            
            # copy the files using hadoop command
            if args.condor:
                  cmd='%s -get %s/red/%s/red/%s/*%02d.fits* %s'%(hdfs,datadir,run,exp,ccd,local)
                  # print cmd
                  result=os.system(cmd)
                  cmd='%s -get %s/red/%s/red/%s/*%02d_cat.fits* %s'%(hdfs,datadir,run,exp,ccd,local)
                  # print cmd
                  result=os.system(cmd)
                  input_dir=local

          
            cmd="{wl_soft}/findstars wl.config +wl_desdm.config +wl_finalcut.config input_prefix={input_dir} root={exp}_{ccd:02d} output_prefix={odir} >>{logfile} 2>&1".format(wl_soft=args.wl_soft,run=run,exp=exp,ccd=ccd,odir=odir,logfile=logfile,input_dir=input_dir)
            # print cmd
            os.system(cmd)
      
            # run measure psf
            # cmd="/astro/u/rarmst/findstars/wl_trunk/src/measurepsf wl.config +wl_desdm.config +wl_finalcut.config input_prefix={input_dir} root={exp}_{ccd:02d} output_prefix={odir}  >>{logfile} 2>&1".format(run=run,exp=exp,ccd=ccd,odir=odir,logfile=logfile,input_dir=input_dir)
            # print cmd
            # os.system(cmd)


      if args.condor:
            #print odir,args.output,os.path.isdir(args.output)
            oodir=args.output+'/'+exp+'/'
            cmd='%s -rm %s/*fz'%(hdfs,local)
            os.system(cmd)
            #print cmd
            cmd='%s -rm %s/*cat.fits'%(hdfs,local)
            os.system(cmd)
           #print cmd
            
            if os.path.isdir(oodir):
                  #print 'mv %s/* %s'%(odir,oodir)
                  os.system('%s -mv %s/* %s'%(hdfs,odir,oodir))
            else:
                  #print 'mv %s %s'%(odir,args.output)
                  os.system('%s -mv %s %s'%(hdfs,args.output))

