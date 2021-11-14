#!/bin/bash

# 01 run the pipeline. if no -t <num_threads> option specified ask about this

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msg "\n    The hfmt_<num> working directory not found.\n    Prepare by running ${Blue}hifimitie init${NC}\n" && exit 2

[ ! -z $1 ] && [ $1 == "-t" ] && is_number $2 && tnum=$2

function thread_setting {
   if [ -z $tnum ]; then

      tnum=$threads; [ -z $tnum ] && tnum=$default_threads

      mach_total=$(getconf _NPROCESSORS_ONLN)
      msg "\n    No -t <num_threads> argument specified for how many maximum threads can be used.\n    There are $mach_total machine threads and the current maximum is $tnum.\n"
      read -p "    Is $tnum maximum thread usage OK [Y|N]? " resp; firstCharUC_resp
      [ "$resp" != "Y" ] && msg "\n    Specify maximum number of threads in the ${Blue}hifimitie run -t <num_threads>${NC} parameter\n" && exit 1
   else
      threads=$tnum
   fi

   update_setting "threads" $threads
   msglog_module "maximum thread usage set to $threads"
   msg ""
}

step_num=1
export stop_run="continue"

function run_step {
   step=$1

   # run the step
   $(get_script $step)

   let "step_num++"
   if [ "$stop_run" = "stop"  ]; then
      msglog_module "pipeline stopped after step ${step_num}: $step"
      exit $step_num
   fi
}

thread_setting

##############################################################
#                     run the pipeline                       #
##############################################################

# 02 run blast on the files in fofn list and create tsv files in hifi_mito_matches directory
run_step blast_to_mito

# 03 create the fasta files with the candidate mito records
run_step pull_fofn_cand_recs

# 04 create mito_feature_sequences.fasta with a record for each mito feature: trna, rrna, gene
run_step select_mito_features

# 05 blast the individual features from the closely related organism to our mito record database
run_step blast_features

# 06 use cmsearch for the trna via mitfi and for the 12S and 16S rna, also add goosehairpins to the output cm_anno file
run_step rna_search
exit
# 07

# 08
$(get_script assemble_w_megahit)

# 09

# 10
