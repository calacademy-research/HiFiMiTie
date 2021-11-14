#!/bin/bash

# source this to get shared functions and variables set

# this is how to use this
: source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, script_dir, msglog, msglog_module, is_number functions and usage among other things
# old way:src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

###### variables ########################################################

src_dir=$(dirname $(realpath $(which hifimitie) ) )

share_loaded="true"

script_dir=${src_dir}/script
module=$(basename $(realpath $0) .sh)

taxdmp_dir=/ccg/db_sets/taxdump

fullnamelineage=${taxdmp_dir}/fullnamelineage.dmp
taxidlineage=${taxdmp_dir}/taxidlineage.dmp

#########################################################################

###### functions ########################################################

# script name can either have.sh extension or not and be in src_dir or script_dir
function get_script {
  scr=$1
  scrsh=${scr}.sh

  [ -f ${script_dir}/$scr ] && echo ${script_dir}/$scr && return
  [ -f ${script_dir}/$scrsh ] && echo ${script_dir}/$scrsh && return

  [ -f ${src_dir}/$scr ] && echo ${src_dir}/$scr && return
  [ -f ${src_dir}/$scrsh ] && echo ${src_dir}/$scrsh && return

  echo ""
}
function msg { # write a msg to stderr
   >&2 echo -e $@
}
function log { # write a msg to a log in the working dir
   [ ! -z $wdir ] && echo -e $@ >> ${wdir}/hfmt.log
}
function msglog { # write to stderr and to logfile
   msg $@
   log $@
}
function  msglog_module { # prefix the msg with the name of the module
   msglog "[$module]: " $@
   return 0
}
function empty_log {
   [ ! -z $wdir ] && printf "" > ${wdir}/hfmt.log
}

function set_wdir {
   unset wdir
   cand=$(ls -td hfmt_[0-9]*/ 2>/dev/null | head -1)  # get most recent of this format, ie hfmt_ followed by a number
   [ ! -z "$cand" ] && [ -d $cand ] && wdir=$cand
   [[ "${wdir}" == */ ]] && wdir="${wdir: : -1}"  # remove slash if it is there
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

############################################################################

set_wdir
usage=$(get_script usage.sh)

source ${script_dir}/clr_vars.sh
