#!/bin/bash

# no existing hfmt dir needed
# check the various scrips and programs that are needed are accessible

# a work in progress

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things

Red="\033[0;31m";    BRed="\033[1;31m";    URed="\033[4;31m"
Green="\033[0;32m";  BGreen="\033[1;32m";  UGreen="\033[4;32m"
Blue="\033[0;34m";   BBlue="\033[1;34m";   UBlue="\033[4;34m"
NC="\033[0m" # No Color (reset)


declare -i pgms_found=0
declare -i pgms_missing=0

function check_for_program {
   loc=$(which $1)
   [ -z $loc ] && echo -e "${BRed}Could not find${NC}: $1" && (( pgms_missing++  )) && return 13

   [ ! -z "$loc" ] && echo -e $(basename $loc)"     \t"$(dirname $loc) && (( pgms_found++ ))
}

# this one can run into dependencies problems itself, check it for these
function check_python_pgm {
   local pgm=$1
   shift
   args="$@"

   loc=$(which $pgm)
   [ -z $loc ] && echo -e "${BRed}Could not find${NC}: $pgm" && (( pgms_missing++ )) && return 13

   printf $(basename $loc)"     \t"$(dirname $loc)
   if $pgm $args >/dev/null 2>/dev/null; then
      printf "\n"
      (( pgms_found++ ))
   else
      printf " was found but can not run due to:\n"
      $pgm $args
      (( pgms_found++ ))
      (( pgms_missing++ ))
   fi
}

echo -e "\nCheck for dependencies needed to run HiFiMiTie (work in progress)\n"

check_for_program hifimitie
ls ${this_dir}/shared.sh >/dev/null && echo -e shared.sh"   \t"${this_dir}
check_for_program blastn
check_for_program cmsearch
check_for_program mitfi.sh
check_for_program mafft
check_for_program bioawk_cas
check_for_program bawk
# check_for_program nocomment  # moved to shared.sh
check_for_program seqkit
check_for_program seqfold
# check_for_program numrecs  # moved to shared.sh
check_for_program trf
check_for_program multi_mitfi.sh
check_python_pgm  edlib_str.py "actg" "actg"
check_python_pgm  mito_analyze.py
check_for_program mitodb_update.sh
check_for_program convert_hfmt_anno_to_gff.sh
check_for_program consensus_from_fasta_alignment.sh
# check_for_program missing_program_example

echo -e "\nFound $pgms_found with $pgms_missing missing or not able to run\n"
