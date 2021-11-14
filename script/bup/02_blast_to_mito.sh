#!/bin/bash

# create the blast tsv files blasting the relevant mitogenome to each of the previously selected PacBio HiFi files
# the files are written to the dir hifi_mito_matches/ in he working dir

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, script_dir, msglog, msglog_module, is_number functions and usage among other things

# check for dirs and files we need
[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2
fofn=${wdir}/files_to_search.fofn
[ ! -f $fofn ] && msg "$fofn file could not be found" && exit 3

# create out_dir if necessary
out_dir="${wdir}/hifi_mito_matches"
[ ! -d $out_dir ] && mkdir $out_dir && msglog_module "created $out_dir to hold the blast output tsv files"
[ ! -d $out_dir ] && msglog_module "directory $out_dir could not be created" && exit 4

# see if the blasting has already been done and exit with a success code of 0 if so, code 25 otherwise
function blasting_completed {
   for tsv in ${out_dir}/*.tsv; do

      prefix=$(basename $tsv .tsv)
      dn=${out_dir}/${prefix}.done

      [ ! -s $dn ] && return 25

   done
   true
}

function set_blast_extra_args {
   unset blast_extra_args
   exargs="-mt_mode 1"

   set_blast_major_minor
   [ $threads -lt 4 ] && return
   [ $blast_major -gt 2 ] && blast_extra_args=$exargs
   [ $blast_major -eq 2 ] && $blast_minor -ge 12 ] && blast_extra_args=$exargs

   unset exargs
}

function do_fofn_blast {
   # check tax id info
   toptax=$(head -1 ${wdir}/toptaxid | cut -f1)

   taxidlist=${wdir}/taxidlist
   taxlidlist_arg="-taxidlist $taxidlist"
   [ ! -f $taxidlist ] && msg "$taxidlist could not be found. Searching all the mitogenomes in the mito db" && taxlidlist_arg=""

   # set up for blasting
   blastncmd=$(get_script blastmax5)
   eval_cutoff=1e-10
   db=mito

   start=$(date +%s)

   # loop through the files in the fofn list and do a blast on each, creating a tsv file in the out_dir for each
   f=0  # add a file counter at end of output file to reduce chance of name collision
   while IFS= read -r file; do

      let "f+=1"
      set_fastx_basename $file
      out_base=${out_dir}/mito_${fastx_basename}_${f}
      out=${out_base}.tsv
      blast_start=$(date +%s)

      msglog_module blastn -db $db -query $file -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -max_target_seqs 5 -num_threads $threads -taxidlist taxidlist -evalue $eval_cutoff $blast_extra_args

      $blastncmd $db $file $threads $taxlidlist_arg -evalue $eval_cutoff $blast_extra_args >$out
      completion_code=$?

      finmsg=$(echo "completed in" $(pprt_till_now $blast_start) "with completion code $completion_code")
      suffix="err" &&[ $completion_code -eq 0 ] && suffix="done"
      echo $finmsg >${out_base}.$suffix

      msglog_module $finmsg

   done < $fofn

   echo
   msglog_module "$f blasts completed in" $(pprt_till_now $start)
}

function make_candidate_set {
   ####### call out to make_mito_rec_cand_tsv.sh with the names of the blast tsv file to create the candidate set
   cand_script=$(get_script make_mito_rec_cand_tsv)
   $cand_script ${out_dir}/mito_*.tsv
}

# if we have the tsv files already complete, then we will skip the rest
blasting_completed && msg "blast tsv files of HiFi reads already created" && exit 0

# do the work, blast HiFi reads to mitogenomes
set_blast_major_minor
do_fofn_blast

# based on matches of the reads to the mitogenomes pull out the records that are likely this organism mitogenome reads
make_candidate_set
