#!/bin/bash

# create hfmt_randnum dir if need be
# handle putting files into files_to_search.fofon in that directory
# handle getting taxid number and when we get it create the taxid list

# this doubles as hifimitie addfiles when the first argument is "-"

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg function and usage among other things

# see if we should just give the usage and exit
[ -z "$1" ] && $usage init && exit 1
[ "$1" == "-" ] && [ "$#" -eq 1 ]  && $usage addfiles && exit 1



###### functions ################################
function get_workdir { # create if need be but return hfmt work dir name in var wdir
   set_wdir # function in shared.sh
   if [ -z $wdir ]; then  # could not find it, create it
      cand=$(echo hfmt_$(echo $RANDOM))
      mkdir $cand && [ -d $cand ] && wdir=$cand
   fi
}

# set pacbio var if it is a PacBio fastq or fasta looking file -- very silly test of movie format name hint.(bawk handles gzip too. also sets fastx if is so
function is_PacBio {
  unset pacbio; unset fastx
  [ -z "$1" ] && return 0; [ ! -f "$1" ] && return 0
   fastx=$(bawk '$name!=""{print "fastx"}{exit}' $1)
  pacbio=$(bawk '$name ~ "^m[0-9]+"_""{print "pacbio"}{exit}' $1) && return true
  return false
}
function handle_file_args {
   addedfiles=""
   for fname in $@; do
      is_PacBio $fname
      [ -z "$pacbio" ] && msg "$fname is not a PacBio fasta or fastq file" && continue
      fullname=$(realpath $fname)
      echo $fullname >>$fofn
      addedfiles="yes"
   done

   # make sure we get rid of any dups since we may call addfiles more than once
   [ ! -z $addedfiles ] && (sort -V $fofn | uniq > $tmpfofn) && mv $tmpfofn $fofn
}
###############################################

###### do the work ############################

taxid_inf=$1 # this could be a number or could be a word or could be "-" whoch means to ignore it for now
shift

# first take care of creating directory if need be, then add files to fofn if they exist and are fastq fasta file that look like from PacBio

get_workdir
[ -z $wdir ] && msg "can not find or create hfmt_<num> directory" && exit 2

fofn=${wdir}/files_to_search.fofn
tmpfofn=${wdir}/tmp.fofn
empty_log


# any additional args after the first should be files

if [ ! -z $1 ]; then  # something that might look like files in our args
   handle_file_args $@
fi

# after file handling get the taxid, if it is a "-" we skip this
[ $taxid_inf == "-" ] && exit 0

# call other script passing it the original first arg
${src_dir}/gettaxid.sh $taxid_inf
