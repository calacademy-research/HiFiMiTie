#!/bin/bash

# create hfmt_todaysdate dir if need be
# handle putting files into files_to_search.fofon in that directory
# handle getting taxid number and when we get it create the taxid list

# this doubles as hifimitie addfiles when the first argument is "-"

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg function and usage among other things

# see if we should just give the usage and exit
[ -z "$1" ] && $usage init && exit 1

# see if we should show addfiles usage and exit if there are no file arguments after the dash "-" arg
[ "$1" == "-" ] && [ "$#" -eq 1 ] && $usage addfiles && exit 1


###### functions ###########################################################################

function get_workdir { # create if need be but return hfmt work dir name in var wdir
   if [ -z $wdir ]; then  # could not find it, create it
      cand=$(echo hfmt_$(date +"%m%d%y"))
      mkdir $cand && [ -d $cand ] && wdir=$cand

      [[ "${wdir}" == */ ]] && wdir="${wdir: : -1}"  # remove slash if it is there
      [ -d $wdir ] && wdir_path=$(realpath $wdir)  # need wdir_path set for msglog_module to work

      msglog_module HiFiMiTie directory $wdir created
   fi
}

function log_fofn_files {
   if [ -s $fofn ]; then
      msglog_module "HiFi file(s) to search for mitochondrial reads:"
      for file in $(cat $fofn); do
         msglog_module "    $file"
      done
   fi
}

# set pacbio var if it is a PacBio fastq or fasta looking file -- very silly test of movie format name hint.(bawk handles gzip too. also sets fastx if is so
# additional test if that fails to see if there are records greater than 10000 bases in the first 1000 records
function is_PacBio {
  unset pacbio; unset fastx
  [ -z "$1" ] && return 0; [ ! -f "$1" ] && return 0
    fastx=$(bawk '$name!=""{print "fastx"}{exit}' $1)
   pacbio=$(bawk 'length($seq)>10000 || $name ~ "^m[0-9]+"_""{print "pacbio"; exit}NR>1000{exit}' $1)
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

   # put the names in our log file
   log_fofn_files
}

function whats_next_msg {
   todo=""
   [ ! -s ${wdir}/taxidlist ] && todo=" choose the taxid of the mitogenomes to search"
   [ ! -z "$todo" ] && todo="${todo} and"
   [ ! -s $fofn ] && todo=" choose HiFi file(s) to search"

   if [ -z "$todo" ]; then
      msg "\nTo run the pipeline:\n\n  \t${Blue}hifimitie run -t <num_threads>${NC}\n\n(hfmt is shorthand for the hifimitie command)\n"
   else
      msg "\n\n  Still need to${BBlue}$todo${NC}:"
      $usage init
   fi
}
##############################################################################################

###### do the work ###########################################################################

taxid_inf=$1 # this could be a number or could be a word or could be "-" which means to ignore it for now
shift

# first take care of creating directory if need be, then add files to fofn if they exist and are fastq fasta file that look like from PacBio

get_workdir
[ -z $wdir ] && msg "can not find or create hfmt_<num> directory" && exit 2

fofn=${wdir}/files_to_search.fofn
tmpfofn=${wdir}/tmp.fofn

# any additional args after the first should be files

if [ ! -z $1 ]; then  # something that might look like files in our args
   handle_file_args $@
fi

# after file handling, get the taxid, however if it is a "-" we skip this
[ $taxid_inf == "-" ] && whats_next_msg && exit 0

# call other script passing it the original first arg
$(get_script gettaxid) $taxid_inf

whats_next_msg
