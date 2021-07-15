#!/bin/bash

#script to submit beamsplitter scans to the queue

numcalcs=$1;
ppn=$2;
mem=$3;

if [ -z $1 ]; then 
  echo "USAGE: push_job_tog <numclacs> [<#cores>] [<memory>]";
  echo "defaults: #cores=16; memory=6000mb"
  exit;
fi

if [ -z $2 ]; then
  ppn="16"
fi

if [ -z $3 ]; then
  mem="6000mb"
fi

WHOAMI=$(whoami)
WORK=$PWD #/home/$WHOAMI/runs/beamsplitter_scans

for (( c=0; c<$numcalcs; c++ ))
do

script=metall_glass_setup$c.pbs

echo "
  #!/bin/bash
  #PBS -l nodes=1:ppn=$ppn -l mem=$mem

  # nasty trick to make the .bashrc recognize this script as an interactive shell
  PS1='\u@\h$ '
  source ~/.bashrc

  date
  hostname

  mkdir -p $WORK

#####################################################
### run python
#####################################################
  
  cd $PWD
  python metallsubstrat.py -j $ppn -c $c -o '$PWD/' -l '$WORK/' -t 60.0 -s '/home/ledwon/fmm_source/'

  date 
  hostname

" > $script

qsub $script

#more $script

rm $script

done
