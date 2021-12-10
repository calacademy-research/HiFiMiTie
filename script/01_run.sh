#!/bin/bash

# 01 run the pipeline. if no -t <num_threads> option specified ask about this

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msg "\n    The hfmt_<num> working directory not found.\n    Prepare by running ${Blue}hifimitie init${NC}\n" && exit 2

[ ! -z $1 ] && [ $1 == "-t" ] && is_number $2 && tnum=$2

function thread_setting {
   cur_threads=$threads
   if [ -z $tnum ]; then

      tnum=$threads; [ -z $tnum ] && tnum=$default_threads

      mach_total=$(getconf _NPROCESSORS_ONLN)
      msg "\n    No -t <num_threads> argument specified for how many maximum threads can be used.\n    There are $mach_total machine threads and the current maximum is $tnum.\n"
      read -p "    Is $tnum maximum thread usage OK [Y|N]? " resp; firstCharUC_resp
      [ "$resp" != "Y" ] && msg "\n    Specify maximum number of threads in the ${Blue}hifimitie run -t <num_threads>${NC} parameter\n" && exit 1

   else
      threads=$tnum
   fi

   if [ "$threads" != "$cur_threads" ]; then
      update_setting "threads" $threads
      msglog_module "maximum thread usage set to $threads"
   fi
}

function show_pipeline_time {
   local run_time="$(get_setting run_time)"
   local long_form_date_start="$(get_setting start)"

   if [ ! -z "$long_form_date_start" ]; then
      start_seconds=$(date -d "$long_form_date_start" +%s)
      time_to_now=$(pprt_till_now $start_seconds)
      [ ! -z "$run_time" ] && time_to_now="$run_time"

      finished="$(get_setting finished)"
      if [ -z "$finished" ]; then
         last_step_completed="$(get_setting step_completed)"
         finished="$last_step_completed"
         update_setting finished "$finished"
      fi

      complete_msg="$time_to_now to complete $finished"

      if [ -z "$run_time" ]; then
         update_setting run_time "$time_to_now"
         msglog ""
         msglog_module "$complete_msg"
      else
         msg ""
         msg "$complete_msg"
      fi
   fi
   msg ""
}

function report_step {
   local step_descr="Step $step_num -- $1"

   # only log it if this is the first time for this step, always output to screen

   max_step=$(get_setting_or_default "step" 1)

   if [ $max_step -ge $step_num ]; then  # only show on screen do not log it
      msg ""
      msg $step_descr
   else
      msglog ""
      msglog_module $step_descr
      update_setting "step" $step_num
   fi
}

# stop_run_file is defined in shared.sh
function check_to_stop_run { # check for the file stop_run in the hfmt wdir, if it is there the caller will top the pipeline from running
   [ -s $wdir/$stop_run_file ] && return 0
   return 1
}

step_num=1

function run_step {
   step=$1
   let "step_num++"

   # run the step
   report_step $step
   local start_step=$(date +%s)
   $(get_script $step)

   if check_to_stop_run; then
      msglog_module "pipeline stopped after step $step due to:\n"
      msglog_file $wdir/$stop_run_file
      msglog ""
      exit $step_num
   fi

   if [ $max_step -lt $step_num ]; then # first time to do this step put time it took into log
      msglog_module "Step $step_num completed in $(pprt_till_now $start_step)"
      update_setting "step_completed" "$(date)"
   fi
}

function begin_run {

   clear_stop_run
   thread_setting

   local max_step=$(get_setting_or_default "step" 1)
   if [ $max_step -gt 1 ]; then # show on screen do not log it
      local num_files=$(wc -l ${wdir_path}/files_to_search.fofn | awk '{print $1}')
      local fstr="file"; [ "$num_files" -gt 1 ] && fstr="files"

      msg "\nStep 1 -- Setup: taxid and HiFi file(s) to use"
      msg "$(get_setting taxid) $(get_setting taxname) and $num_files $fstr chosen with maximum of $threads threads"
   fi
}

##############################################################
#                     run the pipeline                       #
##############################################################

# 01 set thread by cmd line or ask user, show msg if re-running
begin_run

# 02 run blast on the files in fofn list and create tsv files in hifi_mito_matches directory
run_step blast_to_mito

# 03 create the fasta files with the candidate mito records
run_step pull_fofn_cand_recs

# 04 create top_match_feature_sequences.fasta with a record for each mito feature: trna, rrna, gene
run_step select_mito_features

# 05 blast the individual features from the closely related organism to our mito record database
run_step blast_features

# 06 use cmsearch for the trna via mitfi and for the 12S and 16S rna, also add goosehairpins to the output cm_anno file
run_step rna_search

# 07 CR_analysis
run_step CR_analysis

# 08 split sequences based on 05 and 06 annotations and the 07 analysis
run_step split_recs_into_sets

# 09 will become assemble script that calls the other 2 assemble scripts
run_step assemble_mito

# 10 compare assemblies
run_step compare_assemblies

# 11 assemble_CR, stick with msa version or choose one as representative. do some heteroplasmy anaylsis
run_step assemble_CR

# 12 complete the assembly
run_step complete

# 13 show the time it took
show_pipeline_time
