#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=60

# pull the records named in the file mito_rec_candidates.tsv from the fasta or fastq files list in files_to_search.fofn
# blast_to_mito.sh should have been called first so that mito_rec_candidates.tsv has been created

tsv_recs=${wdir}/mito_rec_candidates.tsv
[ ! -f $tsv_recs ] && msg "\nCould not find file $tsv_recs that has the mito candidate record ids\n" && exit 3

fofn=${wdir}/files_to_search.fofn
[ ! -f $fofn ] && msg "\n$fofn containing the names of fasta or fastq files is missing\n" && exit 3

out_fasta=${wdir}/mito_rec_candidates.fasta

function pull_from_file {
   msg pulling records from $1

   bawk -v FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=$FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE '
      FNR==NR{ ar[$name]++; RC = ($11 > $12) ? "RC" : "" ; do_RC[$name] = RC; next }

      ! ($name in ar){ next }

      {QPct = $NF}
      QPct < FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE { filtered++; next }

      do_RC[$name] == "RC" {print ">" $name "_RC QPct:" QPct ; print revcomp($seq); next }
      {print ">" $name " QPct:" QPct; print $seq}

   ' <(prefix_gt $tsv_recs) $1 | fold -w 120
}

function fofn_pull_recs
{
   [ -f $out_fasta ] && printf "" >$out_fasta

   # pull records from each file listed
   while read -r fname; do

      [ -z "$fname" ] && continue
      [ ! -f $fname ] && msg "$fname is not a file, skipping it." && continue

      is_fastx $fname
      [ -z "$fastx" ] && msg "$fname is not a fasta or fastq file" && continue

      pull_from_file $fname >>$out_fasta

   done < $fofn
}

function make_cand_rec_db {
   out_dir="${wdir}/rec_db"
   [ ! -d $out_dir ] && mkdir $out_dir && msg "created $out_dir to hold a blast database containing mito candidate fasta records"
   [ ! -d $out_dir ] && msg "directory $out_dir could not be created" && exit 4

   set_fastx_basename $out_fasta
   dbname=${out_dir}/$fastx_basename

   makeblastdb -in $out_fasta -out $dbname -dbtype nucl -parse_seqids
}

# pull the associated fasta records from each file that has a high enough scoring line in the tsv_rec file, flipping records as needed so all are in same mito direction
fofn_pull_recs

# make a blast database from those records in a subdir named rec_db
make_cand_rec_db
