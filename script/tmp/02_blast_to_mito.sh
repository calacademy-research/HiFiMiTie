#!/bin/bash

# create the blast tsv files blasting the relevant mitogenome to each of the previously selected PacBio HiFi files
# the files are written to the dir hifi_mito_matches/ in he working dir

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, script_dir, msglog, msglog_module, is_number functions and usage among other things

[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2
fofn=${wdir}/files_to_search.fofn
[ ! -f $fofn ] && msg "$fofn file could not be found" && exit 3

toptax=$(head -1 ${wdir}/toptaxid | cut -f1)

taxidlist=${wdir}/taxidlist
taxlidlist_arg="-taxidlist $taxidlist"
[ ! -f $taxidlist ] && msg "$taxidlist could not be found. Searching all the mitogenomes in the mito db" && taxlidlist_arg=""

out_dir="${wdir}/hifi_mito_matches"
[ ! -d $out_dir ] && mkdir $out_dir && msg "created $out_dir to hold the blast output tsv files"
[ ! -d $out_dir ] && msg "directory $out_dir could not be created" && exit 4

blastncmd=$(get_script blastmax5)
eval_cutoff=1e-10

start=$(date +%s)
date "+%T %d%b%Y"

f=0  # add a file counter at end of output file to reduce chance of name collision
while IFS= read -r file; do

   let "f+=1"
   set_fastx_basename $file
   out_base=${out_dir}/mito_${fastx_basename}_${f}
   out=${out_base}.tsv
   blast_start=$(date +%s)

   $blastncmd mito $file $taxlidlist_arg -evalue $eval_cutoff >$out
   completion_code=$?

   finmsg=$(echo "completed in" $(pprt_till_now $blast_start) "with completion code $completion_code")
   suffix="err" &&[ $completion_code -eq 0 ] && suffix="done"
   echo $finmsg >${out_base}.$suffix
   echo $finmsg

done < $fofn

echo
date "+%T %d%b%Y"
echo "$f blasts completed in" $(pprt_till_now $start)

####### call out to make_mito_rec_cand_tsv.sh with the names of the blast tsv file to create the candidate set
${src_dir}/make_mito_rec_cand_tsv.sh ${out_dir}/mito_*.tsv
