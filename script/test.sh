#!/bin/bash

# script to put testing code that has the typical setup

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found"
msg "continuing anyway"

msg "testing with these $# args:" $@

msg "pgmdir: $hifimitie_pgmdir"

msg
msg $PYTHONPATH
msg

mito_analyze.py

msg "show numrecs"
numrecs /home/jhenderson/CAJQFG01.fasta

echo ${exec_path}/${consensus_dir}
