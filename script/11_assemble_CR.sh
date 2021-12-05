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
   update_setting_if_changed assemble_CR_result "msa"

   local m="No tandem repeats found in the $numseqs Control Region sequences and length stddev ${stddev} is less than $stddev_tolerance, consensus will be used."
   if [ -z $trf_already_run ]; then # log it if we ran trf here
      msglog_module $m
   else # otherwise output to screen but not the log file
      msg $m
   fi
}

# all the trf scripts want the input and output to be in the current dir
# this is trf_dir or trf_path which is in the cr_analysis directory
# cr_analysis is where the other results are collected

function search_for_tandem_repeats {
   local trf_report=trf_stats_report.tsv

   trf_report_path=trf_output/$trf_report  # trf_report_path will be used later, so do not make it local
   [ -s $trf_report_path ] && msg "$trf_report_path already created" && trf_already_run=1 && return

   pushd $trf_path >/dev/null

   run_trf.sh ../$CR_fasta_name
   echo "" # no final eol by trf program

   popd >/dev/null
}

function handle_settings {
   prefix=$(get_setting_or_default Primary_CR "CR1")
   stddev=$(get_setting_or_default ${prefix}_stddev 0)

   CR_seqs_w_repeats=$(grep -v "^#" $trf_report_path -c)
   update_setting_if_changed ${prefix}_seqs_w_repeats $CR_seqs_w_repeats
}

function make_softlinks {  # make soft link to CR sequences if not already there -- eg split_sequences/Thr_CR_Pro.fasta
   CR_fasta_file=$(ls ../split_sequences/*_CR_*.fasta | head -1)
   CR_fasta_name=$(basename $CR_fasta_file)
   generic_name=CR_recs_w_flanks.fasta

   [ ! -s $CR_fasta_name ] && ln -s $CR_fasta_file
   [ ! -s $generic_name ]  && ln -s $CR_fasta_name $generic_name

   # make softlinks to the CR consensus frm the msa work so we can compare it with the tandem repeat versions
   CR_consensus_path=${exec_path}/${consensus_dir}/*_CR_*.consensus.fa
   CR_consensus_path=$(echo $CR_consensus_path | sed "s|.*hfmt_[0-9]*/|../|") # make it relative instead of absolute
   CR_consensus=$(basename $CR_consensus_path)
   generic_consensus=CR_recs_w_flanks.consensus.fa

   [ ! -s $CR_consensus ] && ln -s $CR_consensus_path
   [ ! -s $generic_consensus ] && ln -s $CR_consensus $generic_consensus
}

function do_prep {
   # prep the records associated with each copynum of the repeats by creating dirs with relevaant reads under tr_copynum_crs dir
   # the prep shows its work in the cr_dir file tr_copynum_info.tsv, so if it exists we have already done at least the prep work prseumably

   if [ ! -s tr_copynum_info.tsv ] && [ $CR_seqs_w_repeats -gt 0 ]; then
      ${script_dir}/tr_rec_analysis_prep.sh
      show_prep_results
   elif [ -s tr_copynum_info.tsv ]; then
      msg "tr_copynum_info.tsv already created"
   fi
}
function show_prep_results {
   if [ -s tr_copynum_info.tsv ]; then
      msglog_module "$(grep -v ^# tr_copynum_info.tsv -c) sets of tandem repeat records to analyze -- this includes one set of records with no tandem repeats\n"
      msglog_file tr_copynum_info.tsv
      msglog ""
   else
      msglog_module "warning: tr_copynum_info.tsv was not created, indicating no tandem repeat records but we thought there were $CR_seqs_w_repeats records from the trf run"
   fi
}

function make_cr_tr_consensus_anno {
   [ ! -d tr_copynum_crs ] && return

   pushd $exec_path >/dev/null  # exec_path is where we originally executed the hifimitie (hfmt) command

   ${script_dir}/cr_tr_copynum_consensus_anno.sh

   popd >/dev/null
}

cd $cr_dir

make_softlinks  # also sets CR_fasta_file CR_fasta_name generic_name vars
numseqs=$(numrecs $CR_fasta_file)

search_for_tandem_repeats
handle_settings # also sets CR_seqs_w_repeats from the trf report file

[ $CR_seqs_w_repeats -eq 0 ] && [ $stddev -lt $stddev_tolerance ] && log_use_msa && exit

do_prep  # if we get here we should have some seqs with repeats to deal with

make_cr_tr_consensus_anno  # from each set of records with the same tr copy num, make a consensus fasta, annotate it, and compare it with the full consensus CR
