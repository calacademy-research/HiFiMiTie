#!/bin/bash

# main job is to run trf and analyze any tandem repeats for annotation and heteroplasmy discovery
# then choose one or more representative CR sequences (next step uses this to create final assembly -- which may be the megahit or msa)

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

trf_path=$(realpath $trf_dir)
[ ! -d $trf_path ] && mkdir $trf_path

declare -i stddev=0
stddev_tolerance=20

function log_use_msa {
   update_setting assemble_CR_result "msa"
   msglog_module "No tandem repeats found in the $numseqs Control Region sequences and length stddev ${stddev} is less than $stddev_tolerance, consensus will be used."
}

# all the trf scripts want the input and output to be in the current dir
# this is trf_dir or trf_path which is in the cr_analysis directory
# cr_analysis is where the other results are collected

# make soft link to CR sequences if not already there -- eg split_sequences/Thr_CR_Pro.fasta

cd $cr_dir

CR_fasta_file=$(ls ../split_sequences/*_CR_*.fasta | head -1)
CR_fasta_name=$(basename $CR_fasta_file)
[ ! -s $CR_fasta_name ] && ln -s $CR_fasta_file

numseqs=$(numrecs $CR_fasta_file)

cd $trf_path

[ -s trf_stats.txt ] && msg "trf_output/trf_stats.txt already created" && exit

run_trf.sh ../$CR_fasta_name
echo "" # no final eol by trf program

prefix=$(get_setting_or_default Primary_CR "CR1")
stddev=$(get_setting_or_default ${prefix}_stddev 0)

CR_seqs_w_repeats=$(grep "^Sequences" trf_stats.txt -c)
update_setting_if_changed ${prefix}_seqs_w_repeats $CR_seqs_w_repeats

[ $CR_seqs_w_repeats -eq 0 ] && [ $stddev -lt $stddev_tolerance ] && log_use_msa && exit
