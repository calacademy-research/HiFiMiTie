#!/bin/bash

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, wdir_path, script_dir, msglog, msglog_module, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

# add create_dist_file and get_stats_from_dist_file and their support functions
source ${script_dir}/cr_shared_funcs.sh
source ${script_dir}/first_trna_from_features_or_taxname.sh

CR_str="Control Region"

######################################################################
#                       function definitions                         #
######################################################################

function update_CR_settings {
   beg_trna=$1; end_trna=$2

   (( num_CRs++ ))

   update_setting_if_changed "num_CRs" $num_CRs
   crtag="CR${num_CRs}"

   type="$CR_str"
   if (( $mean < 1000 )); then
      type=remnant
      (( num_remnant++ ))
   elif (( (num_CRs - num_remnant) > 1 )); then
      type="Duplicate $CR_str"
      update_setting_if_changed "CR1_type" "$type"
   fi

   if (( mean > largest_mean )); then
      largest_mean=$mean
   fi

   update_setting_if_changed "Primary_CR" "$crtag"
   update_setting_if_changed ${crtag}_flanks "$beg_trna $end_trna"
   update_setting_if_changed ${crtag}_mean $mean
   update_setting_if_changed ${crtag}_stddev $stddev
   update_setting_if_changed ${crtag}_recs_w_CR $recs
   update_setting_if_changed ${crtag}_type "$type"

   msglog_module "$crtag has type $type, with mean length ${mean} bp found from $recs reads containing it. flanking trnas $beg_trna and $end_trna"
}

function log_stat_file {
   # log the stats
   [ -s $stat_file ] && msglog "" && msglog $(basename $stat_file) && msglog_file $stat_file && msglog ""
}

function create_wrap_around_cr_dist_file {
   create_dist_file $last_rna $first_rna
}

function get_wrap_around_cr_stats {  # 03Jan2023 change from trna to rna to accomodate 12S CR flank
   get_stats_from_dist_file $last_rna $first_rna
   update_CR_settings $last_rna $first_rna

   log_stat_file
}

function set_trna_vars {
   first_trna=$(get_setting_or_default "first_trna" "F")
   last_trna=$(get_setting_or_default "last_trna" "P")

   # working toward allowing 12S as last_trna
   first_rna=$(get_setting_or_default "first_rna" $first_trna)
   last_rna=$(get_setting_or_default "last_rna" $last_trna)

   gh_prev_trna=$(get_setting_or_default "gh_prev_trna" $last_trna)
   gh_succ_trna=$(get_setting_or_default "gh_succ_trna" $first_trna)
}

# if the goosehairpin flanking trnas are different from the wrap around trnas, make a distance file and stats for it
function gh_trnas_are_first_last {
   set_trna_vars

   [ $gh_prev_trna = $last_trna ] && [ $gh_succ_trna = $first_trna ] && return 0  # return true since they are same as first last

   false
}

function create_gh_cr_dist_file {
   create_dist_file $gh_prev_trna $gh_succ_trna
}

function get_gh_cr_stats {
   get_stats_from_dist_file $gh_prev_trna $gh_succ_trna
   update_CR_settings $gh_prev_trna $gh_succ_trna

   log_stat_file
}

#######################################################################
#                                                                     #
#   call the functions if the files they create are not in $cr_dir    #
#                                                                     #
#######################################################################

declare -i num_remnant=0
declare -i num_CRs=0  # incremented in update_CR_settings
declare -i largest_mean=0

[ ! -d $cr_dir ] && mkdir $cr_dir && msglog_module "created $(basename $cr_dir) to hold control region(s) analysis files"
[ ! -d $cr_dir ] && msglog_module "directory $cr_dir could not be created" && exit 5

# this uses the first trna in the feature list of the best match (usually F but can be others; eg in Arthropoda Lepidoptera is M)
get_first_trna_setting
set_trna_vars  # 03Jan2023 also set first_rna and last_rna

set_cr_filenames $last_rna $first_rna  # 03Jan2023 change to use rna instread of trna -- tho unless last_rna is 12S they are the same
run_if_no_file create_wrap_around_cr_dist_file  $dist_file
run_if_no_file get_wrap_around_cr_stats         $stat_file

if ! gh_trnas_are_first_last; then
   set_cr_filenames $gh_prev_trna $gh_succ_trna
   run_if_no_file create_gh_cr_dist_file  $dist_file
   run_if_no_file get_gh_cr_stats         $stat_file
fi

[ $num_CRs = 1 ] && update_setting_if_changed "CR1_type" "$CR_str"
