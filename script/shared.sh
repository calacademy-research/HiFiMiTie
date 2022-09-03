#!/bin/bash

# source this to get shared functions and variables set

# this is how to use this
: source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, wdir_path, script_dir, msglog, msglog_module, is_number functions and usage among other things
# old way:src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

###### variables ########################################################

src_dir=$(dirname $(realpath $(which hifimitie) ) )

shared_loaded="true"

# replace these the lines in usage.sh when the title or version number changes
hfmt_version=0.03
hfmt_version_date=05-Dec-2021
hfmt_title="HiFiMiTie version $hfmt_version -- Find & Analyze Metazoan Mitochondria from HiFi reads"

script_dir=${src_dir}/script
module=$(basename $(realpath $0) .sh | sed "s/^[0-9][0-9]_//")

data_dir=${src_dir}/data

# when this file is created in wdir via function stop_run, the run command will stop after the current step
stop_run_file="STOP_RUN_DUE_TO_ERROR"

: path from which we run this is in $exec_path

#########################################################################

###### functions ########################################################

# script name can either have .sh extension or not and be in src_dir or script_dir
# also let a plain script name match .sh extesnion and also a 2 digit prefix followed by underscore -- blast_to_mito matches 02_blast_to_mito.sh
function get_script {
  local scr=$1
  local scrsh=${scr}.sh
  local numscrsh="[0-1][0-9]_*${scrsh}"

  [ -f ${script_dir}/$scr ] && echo ${script_dir}/$scr && return
  [ -f ${script_dir}/$scrsh ] && echo ${script_dir}/$scrsh && return
  [ -f ${script_dir}/$numscrsh ] && echo ${script_dir}/$numscrsh && return

  [ -f ${src_dir}/$scr ] && echo ${src_dir}/$scr && return
  [ -f ${src_dir}/$scrsh ] && echo ${src_dir}/$scrsh && return
  [ -f ${src_dir}/$numscrsh ] && echo ${src_dir}/$numscrsh && return

  echo ""
}
function msg { # write a msg to stderr
   >&2 echo -e "$@"
}
function log { # write a msg to a log in the working dir
   [ ! -z $wdir_path ] && echo -e "$@" >> ${wdir_path}/hfmt.log
}
function msglog { # write to stderr and to logfile
   msg "$@"
   log "$@"
}
function msglog_module { # prefix the msg with the time and the name of the module
   set_now
#   change to not include module name just time
#   msglog "[$now]\t$module:" "$@"
   msglog "[$now]\t$@"
   return 0
}
function empty_log {
   [ ! -z $wdir_path ] && printf "" > ${wdir_path}/hfmt.log
}
function msglog_file { # output file lines to stderr and to logfile
   local file=$1
   [ ! -s $file ] && return 1

   while IFS= read -r line; do
     msglog "$line"
   done < "$file"
}

function set_wdir {  # wdir is just the directory name, wdir_path is the fullpath upto and including the directory name without a slash at the end
   unset wdir; unset wdir_path

   # check if we are in the hfmt dir and move up to its parent, so things work, and we are not penalized for starting a run in the hfmt directory
   if [[ $(basename $(pwd)) =~ hfmt_[0-9]* ]]; then  # need double brackets so asterisk is used as wildcard in string and not for file globbing
      cd ..
      msg moving up to $(pwd)
   fi

   exec_path=$(pwd)

   cand=$(ls -td hfmt_[0-9]*/ 2>/dev/null | head -1)  # get most recent of this format, ie hfmt_ followed by a number
   [ ! -z "$cand" ] && [ -d $cand ] && wdir=$cand

   [[ "${wdir}" == */ ]] && wdir="${wdir: : -1}"  # remove slash if it is there

   # wdir is just the directory name. set wdir_path to the complete path of wdir
   [ ! -z $wdir ] && wdir_path=$(realpath $wdir)
}

function nocomment {
   awk '/^#/{next}{print}' $@
}
function noncomment {
   nocomment $@
}

function is_number { # checks if argument is a positive integer, return value is true if a number false otherwise
   re='^[0-9]+$'
   [[ "$1" =~ $re ]]
}

function firstCharUC_resp { # return in the var $resp the first char of argument $1 uppercased
   resp=${resp:0:1}; resp=$(echo $resp |tr '[:lower:]' '[:upper:]')
}

# remove all the usual ends on a fasta fastq gzipped or not file
function set_fastx_basename {
   unset fastx_basename
   [ -z $1 ] && return
   fastx_basename=$(basename $1 .gz)   # get rid of any gzip suffix
   fastx_basename=$(basename $fastx_basename .fa)
   fastx_basename=$(basename $fastx_basename .fq)
   fastx_basename=$(basename $fastx_basename .fna)
   fastx_basename=$(basename $fastx_basename .fasta)
   fastx_basename=$(basename $fastx_basename .fastq)
}
function get_fastx_basename {  # so we can use $() to assign result to var with 1 line
   set_fastx_basename $1
   echo $fastx_basename
}

function is_fastx {  # set fastx var if it is a fastq or fasta looking file
  [ -z "$1" ] && fastx="" && return 0
  fastx=$(bawk '$name!=""{print "fastx"}{exit}' $1)
}

function numrecs {
   bioawk -c fastx '
     NR==1{LastFilename = FILENAME}
     FNR==1{if(NR>1){print LastFilename": " FileRecs} FileRecs=0; LastFilename=FILENAME}

     {FileRecs++}

     END{if(FNR!=NR){printf "%s: ", FILENAME}; print FileRecs}' "$@"
}

function sorttab { # can never remember how to use the tab as the separator, so put in function
   sort -t$'\t' $@
}

# function to set blast major and minor version number variables
# we can use this iwth version 2.12 or greater to use the new -mt_mode 1 option
# this will use the threads among the queries instead of partitioning the database search
# with a -taxidlist in use this is probably faster

function set_blast_major_minor {
   blast_version=$(blastn -version | awk '/[1-9]/{sub("[^\\.0-9]*",""); print; exit}')
   blast_major=$(echo $blast_version | awk '{split($0,ar,"\\."); major=int(ar[1]); print major; exit}')
   blast_minor=$(echo $blast_version | awk '{split($0,ar,"\\."); minor=int(ar[2]); print minor; exit}')
}

function one_letter_AA { # call with 1 letter or 3 letter version of AA, returns the 1 letter version
   AAsyms.sh $1 | awk '{print $1}'
}
function three_letter_AA { # call with 3 letter or 1 letter version of AA, returns the 3 letter version
   AAsyms.sh $1 | awk '{print $2}'
}

# called with mitfi file path in $1 and single letter tRNA code in $2 and column value to retrieve in $3
function get_mitfi_entry_col {
   awk -v letter=$2 -v col=$3 'BEGIN{FS="\t"}
      /^#/{next}
      {++entry}
      $7==letter {
         if(col==0){print entry}
         else{print $col}
      exit
   }' $1
}
function get_mitfi_entry_info { # called with $1 mitfi file path, $2 trna letter code
   unset tbeg tend tentry
   tbeg=$(get_mitfi_entry_col $1 $2 2)
   tend=$(get_mitfi_entry_col $1 $2 3)
   tentry=$(get_mitfi_entry_col $1 $2 0)
}

function search_OL {  # specify fasta file and output dir if not the cm_search result
   local fasta=$1
   local cm_model=${src_dir}/models/OL.cm
   [ -z $fasta ] && fasta=${wdir_path}/mito_hifi_recs.fasta

   local outdir=$2
   [ -z $outdir ] && outdir=$cm_dir

   # 08Aug2022 change from --cpu 1 to --cpu $threads
   cmsearch --noF4b --cpu $threads --notextw -E 0.01 --mxsize 80000 --tblout ${outdir}/OL.tbl $cm_model $fasta >${outdir}/OL.cmout
}

function make_mito_rec_candidates_tsv { # make the combined candidate tsv from the individual tsvs in hifi_mito_matches
   local tsv_recs="$@"
   local out_dir=${wdir}/hifi_mito_matches

   [ -z "$tsv_recs" ] && tsv_recs=${out_dir}/*.tsv
   [[ ! "$tsv_recs" == *"hifi_mito_matches/"* ]] && tsv_recs=${out_dir}/*.tsv

   sorttab -k1,1 -k12,12nr $tsv_recs | tophit | awk '$NF>9' | sorttab -k17,17nr > ${wdir}/mito_rec_candidates.tsv
}

# give duration of seconds in HMS or DHMS format
function pprt_seconds {
   awk -v seconds=$1 '
      function pprt_tm(S    ,s,m,h) {
          s=(S%60);  tm = (s>9) ? s"s" : "0"s"s"
          if(S >= 60) {
              m=int(S/60)
              if(m >= 60) {
                 h=int(m/60)
                 if(int(h/24)>0) { # add days to format
                    d = int(h/24)"d"; h = h - (d*24)
                 } else d = ""
                 m = d h "h" (m % 60)
              }
              else { m="0h"m}
              tm = m"m" tm
          } else tm = "0h0m" tm
          return tm
      }
      BEGIN{ print pprt_tm(seconds) }
   '
}
function pprt_till_now {
  now=$(date +%s)
  start=$1
  let duration=$now-$start
  pprt_seconds $duration
}

function set_now {
   now=$(date "+%Y-%m-%d %T")
}

# set default_threads var to 1/4 of the machine threads up to a max of 32
function set_default_threads {
   default_threads=1

   mach_total=$(getconf _NPROCESSORS_ONLN)
   is_number $mach_total && default_threads=$(awk -v mach_total=$mach_total 'BEGIN{ t = int( (mach_total+3)/4 ); if(t>32)t=32; print t}')
}

function settings_file {
   echo ${wdir_path}/settings.tsv
}

function update_setting {
   local key=$1; local value=$2
   [ -z "$key" ] && return 1

   local file=$(settings_file)
   [ ! -f $file ] && echo -e "$key\t$value" >$file && return

   awk -v key="$key" -v value="$value" '
      BEGIN{FS="\t"; OFS="\t"}
      $1==key{print key, value; reset=1; next}
      { print }
      END{if(!reset){print key, value}}  # was not in there, add it at the end of the file
   ' $file > tmp_settings
   mv tmp_settings $file
}

function get_setting {
   key=$1
   [ -z "$key" ] && return 1

   file=$(settings_file)
   [ ! -f $file ] && return

   awk -v key="$key" 'BEGIN{FS="\t"; OFS="\t"}
        $1==key{print $2; exit}
   ' $file
}

function get_setting_or_default { # $1 is name of setting to get and $2 is default to return if it is not found
   key=$1
   [ -z "$key" ] && return 1

   file=$(settings_file)
   [ ! -f $file ] && return

   awk -v key="$key" -v dflt="$2" '
        BEGIN{value = dflt; FS="\t"; OFS="\t"}
        $1==key{value = $2; exit}
        END{if(value) print value}
   ' $file
}

function update_setting_if_changed { # checks current value and if it is missing or different from value in $2 call update_setting, otherwise do nothing
   key=$1; value=$2
   [ -z "$key" ] && return 1

   cur_val=$(get_setting $key)
   [ "$cur_val" == "$value" ] && return

   update_setting $key "$value"
}

function get_threads {
   set_default_threads
   get_setting_or_default "threads" $default_threads
}

# function to safely call if a function is needed as an arg
function noop {
   return 0
}

# if there is an existing file in arg 2, we output a msg and return, otherwise run arg 1
function run_if_no_file {
   to_run="$1"
   [ ! -z $2 ] && [ -s $2 ] && msg $(basename $2) already created && return

   $to_run
}

function mkdir_if_needed {
  [ ! -d $1 ] && mkdir $1
}

function feature_sequence_file {
   echo ${wdir_path}/top_match_feature_sequences.fasta
}

function anno_start_item_distribution_file {
   echo ${cm_dir}/higher_qual_anno_start_item_distribution.txt
}

function set_mitodb { # 11Nov2021 use new vars for sytem settings to get mito db path
   local mitodir=$(get_setting mitodb_dir)
   local mitoname=$(get_setting mitodb_name)

   mitodb_path=${mitodir}/${mitoname}
}

function set_taxonomy_vars {
   taxdmp_dir=$(get_setting taxonomy_dir)   #/ccg/db_sets/taxdump
   [ -z $taxdmp_dir ] && return  # this should be valid after the init run that creates the working hfmt directory

   fullnamelineage=${taxdmp_dir}/fullnamelineage.dmp
   taxidlineage=${taxdmp_dir}/taxidlineage.dmp
   taxnodes=${taxdmp_dir}/nodes.dmp
}

function set_dir_vars {  # eventually should set them all here, for now just some
   blast_dir=${wdir}/blast_results
   cr_dir=${wdir}/cr_analysis
   cm_dir=${wdir}/cm_results
   splitseq_dir=${wdir}/split_sequences
   megahit_dir=${wdir}/megahit_out
   alignasm_dir=${wdir}/msa_assembly
   msa_dir=${alignasm_dir}/msa
   consensus_dir=${alignasm_dir}/consensus
   aligncm_dir=${alignasm_dir}/cm_mitfi
   final_dir=${wdir}/final_results
   compare_dir=${wdir}/compare_megahit_msa
   trf_dir=${cr_dir}/trf_output
   complete_dir=${wdir}/complete
}

function stop_run {
   reason="$@"
   [ -z "$reason" ] && reason="unspecified reason to stop"

   echo "$reason" > $wdir/$stop_run_file
}
function clear_stop_run {  # delete the file so we can try agin
   [ -f $wdir/$stop_run_file ] && rm $wdir/$stop_run_file
}

############################################################################

set_wdir
set_taxonomy_vars
set_dir_vars

threads=$(get_threads)
update_setting_if_changed "threads" $threads

usage=$(get_script usage.sh)

source $(get_script clr_vars.sh)
