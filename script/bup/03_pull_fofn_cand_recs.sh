#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

# pull the records named in the file mito_rec_candidates.tsv from the fasta or fastq files list in files_to_search.fofn
# blast_to_mito.sh should have been called first so that mito_rec_candidates.tsv has been created

out_fasta=${wdir}/mito_rec_candidates.fasta

function check_input_files {
   tsv_recs=${wdir}/mito_rec_candidates.tsv
   [ ! -f $tsv_recs ] && msglog_module "\nCould not find file $tsv_recs that has the mito candidate record ids\n" && exit 3

   fofn=${wdir}/files_to_search.fofn
   [ ! -f $fofn ] && msglog_module "\n$fofn containing the names of fasta or fastq files is missing\n" && exit 3
}

FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=60

function pull_from_file {
   msglog_module pulling records from $1

   bawk -v FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=$FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE -v wdir=$wdir '

      FNR==NR{ ar[$name]++; QPct[$name]=$NF; RC = ($11 > $12) ? "RC" : "" ; do_RC[$name] = RC; next }

      ! ($name in ar){ next }
      QPct[$name] < FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE { filtered++; next }

      do_RC[$name] == "RC" {print ">" $name "_RC QPct:" QPct[$name] ; print revcomp($seq); pulled++; next }
      {print ">" $name " QPct:" QPct[$name]; print $seq; pulled++}

   END { m = sprintf("   pulled %d records, excluded %d with low coverage.", pulled, filtered)
        # print m > "/dev/stderr"
         pfc_msg = wdir "/pfc_msg"
         print m > pfc_msg
   } ' <(prefix_gt $tsv_recs) $1 | fold -w 120
}

function fofn_pull_recs
{
   [ -f $out_fasta ] && printf "" >$out_fasta
   pfc_file=$wdir/pfc_msg; rm -f $pfc_file

   # pull records from each file listed
   while read -r fname; do

      [ -z "$fname" ] && continue
      [ ! -f $fname ] && msglog_module "$fname is not a file, skipping it." && continue

      is_fastx $fname
      [ -z "$fastx" ] && msglog_module "$fname is not a fasta or fastq file" && continue

      pull_from_file $fname >>$out_fasta
      [ -f $pfc_file ] && msglog_module "$(head -1 $pfc_file)" && rm -f $pfc_file

   done < $fofn
}

function rec_db_exists {
   dbname=$1
   blastdbcmd -db $dbname -info >/dev/null 2>/dev/null
   cmpl_code=$?
   return $cmpl_code
}

function make_cand_rec_db {
   out_dir="${wdir}/rec_db"
   [ ! -d $out_dir ] && mkdir $out_dir && msglog_module "created $out_dir to hold a blast database containing mito candidate fasta records"
   [ ! -d $out_dir ] && msglog_module "directory $out_dir could not be created" && exit 4

   set_fastx_basename $out_fasta
   dbname=${out_dir}/$fastx_basename

   if ! rec_db_exists $dbname; then
      makeblastdb -in $out_fasta -out $dbname -dbtype nucl -parse_seqids
   else
      msg "$fastx_basename database already exists"
   fi
}

# job of these is to create the ${wdir}/mito_rec_candidates.fasta file
if [ ! -s $out_fasta ]; then
   # make sure we have input necessary
   check_input_files

   # pull the associated fasta records from each file that has a high enough scoring line in the tsv_rec file, flipping records as needed so all are in same mito direction
   fofn_pull_recs
else
   msg "$(basename $out_fasta) already created with $(numrecs $out_fasta) records"
fi

# make a blast database from those records in a subdir named rec_db
make_cand_rec_db
