#!/bin/bash

# source this to get shared functions and variables set

# this is how to use this
: source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, wdir_path, script_dir, msglog, msglog_module, is_number functions and usage among other things
# old way:src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

###### variables ########################################################

src_dir=$(dirname $(realpath $(which hifimitie) ) )

share_loaded="true"

script_dir=${src_dir}/script
module=$(basename $(realpath $0) .sh | sed "s/^[0-9][0-9]_//")

taxdmp_dir=/ccg/db_sets/taxdump

fullnamelineage=${taxdmp_dir}/fullnamelineage.dmp
taxidlineage=${taxdmp_dir}/taxidlineage.dmp

#########################################################################

###### functions ########################################################

# script name can either have .sh extension or not and be in src_dir or script_dir
# also let a plain script name match .sh extesnion and also a 2 digit prefix followed by underscore -- blast_to_mito matches 02_blast_to_mito.sh
function get_script {
  scr=$1
  scrsh=${scr}.sh
  numscrsh="[0-1][0-9]_*${scrsh}"

  [ -f ${script_dir}/$scr ] && echo ${script_dir}/$scr && return
  [ -f ${script_dir}/$scrsh ] && echo ${script_dir}/$scrsh && return
  [ -f ${script_dir}/$numscrsh ] && echo ${script_dir}/$numscrsh && return

  [ -f ${src_dir}/$scr ] && echo ${src_dir}/$scr && return
  [ -f ${src_dir}/$scrsh ] && echo ${src_dir}/$scrsh && return
  [ -f ${src_dir}/$numscrsh ] && echo ${src_dir}/$numscrsh && return

  echo ""
}
function msg { # write a msg to stderr
   >&2 echo -e $@
}
function log { # write a msg to a log in the working dir
   [ ! -z $wdir_path ] && echo -e $@ >> ${wdir_path}/hfmt.log
}
function msglog { # write to stderr and to logfile
   msg $@
   log $@
}
function msglog_module { # prefix the msg with the time and the name of the module
   set_now
   msglog "$now $module: " $@
   return 0
}
function empty_log {
   [ ! -z $wdir_path ] && printf "" > ${wdir_path}/hfmt.log
}

function set_wdir {  # wdir is just the directory name, wdir_path is the fullpath uptp and including the directory name wothout a slash at the end
   unset wdir; unset wdir_path

   # check if we are in the hfmt dir and move up to its parent, so things work, and we are not penalized for starting a run in the hfmt directory
   if [[ $(basename $(pwd)) =~ hfmt_[0-9]* ]]; then  # need double brackets so asterisk is used as wildcard in string and not for file globbing
      cd ..
   fi

   cand=$(ls -td hfmt_[0-9]*/ 2>/dev/null | head -1)  # get most recent of this format, ie hfmt_ followed by a number
   [ ! -z "$cand" ] && [ -d $cand ] && wdir=$cand

   [[ "${wdir}" == */ ]] && wdir="${wdir: : -1}"  # remove slash if it is there

   # wdir is just the directory name. set wdir_path to the complete path of wdir
   [ ! -z $wdir ] && wdir_path=$(realpath $wdir)
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
   fastx_basename=$(basename $1 .gz)   # get rid of any gzip suffix
   fastx_basename=$(basename $fastx_basename .fa)
   fastx_basename=$(basename $fastx_basename .fq)
   fastx_basename=$(basename $fastx_basename .fna)
   fastx_basename=$(basename $fastx_basename .fasta)
   fastx_basename=$(basename $fastx_basename .fastq)
}

function is_fastx {  # set fastx var if it is a fastq or fasta looking file
  [ -z "$1" ] && fastx="" && return 0
  fastx=$(bawk '$name!=""{print "fastx"}{exit}' $1)
}

# function to set blast major and mino version number variables
# we can use this iwth version 2.12 or greater to use the new -mt_mode 1 option
# this will use the threads among the queries instead of partitioning the database search
# with a -taxidlist in use this is probably faster

function set_blast_major_minor {
   blast_version=$(blastn -version | awk '/[1-9]/{sub("[^\\.0-9]*",""); print; exit}')
   blast_major=$(echo $blast_version | awk '{split($0,ar,"\\."); major=int(ar[1]); print major; exit}')
   blast_minor=$(echo $blast_version | awk '{split($0,ar,"\\."); minor=int(ar[2]); print minor; exit}')
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
   key=$1; value=$2
   [ -z "$key" ] && return 1

   file=$(settings_file)
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

############################################################################

set_wdir
threads=$(get_threads)

usage=$(get_script usage.sh)

source ${script_dir}/clr_vars.sh
