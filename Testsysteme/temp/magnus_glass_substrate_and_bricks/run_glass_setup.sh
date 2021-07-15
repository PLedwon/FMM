#!/bin/bash

#script to submit beamsplitter scans to the queue

ppn=$1;
mem=$2;


if [ -z $1 ]; then
  echo "USAGE: run_glass_setup <#cores> [<memory>]";
  echo "defaults: memory=6000mb"
  exit;
fi

if [ -z $2 ]; then
  mem="6000mb"
fi

WHOAMI=$(whoami)
WORK=$PWD #/home/$WHOAMI/pbs_runs/glass_setup

script=glass_setup.pbs

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
  python calc_glass_setup.py -j $ppn -o '$PWD/' -l '$WORK/' -t 30.0 -s '/home/ledwon/fmm_source/'

  date 
  hostname

" > $script

qsub $script

#more $script

rm $script

