#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

# pull the records named in the file mito_rec_candidates.tsv from the fasta or fastq files list in files_to_search.fofn
# blast_to_mito.sh should have been called first so that mito_rec_candidates.tsv has been created

out_fasta=${wdir}/mito_hifi_recs.fasta
candidate_fasta=${wdir}/mito_hifi_candidates.fasta
cand_trf_dir=trf_cand_check

function check_input_files {
   tsv_recs=${wdir}/mito_rec_candidates.tsv
   [ ! -f $tsv_recs ] && stop_run "Could not find file $tsv_recs that has the mito candidate record ids" && exit 3

   fofn=${wdir}/files_to_search.fofn
   [ ! -f $fofn ] && stop_run "$fofn containing the names of fasta or fastq files is missing" && exit 3
}

# 06Dec start to change what we do here:
# we now will pull all records blasting to something in the tsv and process them as candidates.
# this file will be in same directory as the tsvs

# Any of these recs with at least qcov_min we consider mito recs (this is all we did before)
# To accommodate the fact that this stat is skewed when we have large CR sequences
# in the recs we are processing, we will look for OH and gh (goose hairpins)
# and if either is present include that record
# TTO DO: as above states

qcov_min=60
qcov_at_least=50
FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=$qcov_min
KEEP_RECS_WITH_THIS_QUERY_COV=$qcov_min

function pull_candidates_from_file {
   file_to_pull_from=$1
   msglog_module pulling records from $file_to_pull_from

   bawk -v qcov_at_least=$qcov_at_least -v wdir=$wdir '

      FNR==NR {
         tsv_lines++

         QryCov = $NF
         if (QryCov < qcov_at_least)
            will_be_skipped++

         ar[$name]++; QPct[$name] = QryCov
         RC = ($11 > $12) ? "RC" : "" ; do_RC[$name] = RC

         next
      }

      FNR==1 { # check if no work is needed and we will skip looking through the file
         if (will_be_skipped == tsv_lines)
            exit
      }

      ! ($name in ar){ next }
      QPct[$name] < qcov_at_least {
         filtered++
         next
      }

      do_RC[$name] == "RC" {print ">" $name "_RC QPct:" QPct[$name] ; print revcomp($seq); pulled++
         next
      }

      {print ">" $name " QPct:" QPct[$name]; print $seq; pulled++}

   END {
         m = sprintf("   pulled %d candidate records (filtered %d)", pulled, filtered)
         pfc_msg = wdir "/pfc_msg"
         print m > pfc_msg
   } ' <(prefix_gt $tsv_recs) $file_to_pull_from | fold -w 120
}

function pull_from_file {  # this one throws away too many candidate by using only QryCov percentage
   msglog_module pulling records from $1

   bawk -v FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE=$FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE -v wdir=$wdir '

      FNR==NR {
         tsv_lines++

         QryCov = $NF
         if (QryCov < FILTER_OUT_RECS_LESS_THAN_THIS_QUERY_COVERAGE)
            will_be_skipped++

         ar[$name]++; QPct[$name] = QryCov
         RC = ($11 > $12) ? "RC" : "" ; do_RC[$name] = RC

         next
      }

      FNR==1 { # check if no work is needed and we will skip looking through the file
         if (will_be_skipped == tsv_lines)
            exit
      }

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
   [ -f $candidate_fasta ] && printf "" >$candidate_fasta  # new way candidates here filter to new set later
#   [ -f $out_fasta ] && printf "" >$out_fasta
   pfc_file=$wdir/pfc_msg; rm -f $pfc_file

   # pull records from each file listed
   while read -r fname; do

      [ -z "$fname" ] && continue
      [ ! -f $fname ] && msglog_module "$fname is not a file, skipping it." && continue

      is_fastx $fname
      [ -z "$fastx" ] && msglog_module "$fname is not a fasta or fastq file" && continue

      pull_candidates_from_file $fname >>$candidate_fasta
      [ -f $pfc_file ] && msglog_module "$(head -1 $pfc_file)" && rm -f $pfc_file

#      pull_from_file $fname >>$out_fasta
#      [ -f $pfc_file ] && msglog_module "$(head -1 $pfc_file)" && rm -f $pfc_file

   done < $fofn
}

# use the Qcov percentage, as before, also check for goosehairpin and tandem repeats
function filter_recs_from_candidates {
   local mito_hifi_recs=${out_fasta}

   [ -s $mito_hifi_recs ] && msg "$mito_hifi_recs already created" && return

   awk -v KEEP_RECS_WITH_THIS_QUERY_COV=$KEEP_RECS_WITH_THIS_QUERY_COV '
      function goosehairpin_pos() {
         return match($2, /CCCCCCC[AGT]{1,3}CCCCCCC/)
      }
      function qcov_pct_from_header() { # will be in comment like this QPct:75
         qcp = index($3, "QPct:")
         if (qcp < 1) return 0  # could not QPct in comment
         pct = substr($3, qcp+5)
         return int(pct)
      }
      # file 1 has the names of those recs with tandem repeats in field 1
      FNR && /^#/ { next }
      FNR==NR {
         tr_recs[$1]++
         next
      }

      { qcov_pct=qcov_pct_from_header(); gh_pos=goosehairpin_pos() }

      { keep = qcov_pct >= KEEP_RECS_WITH_THIS_QUERY_COV }
      { keep = keep || gh_pos || $1 in tr_recs }

      keep {
         addtl = (gh_pos) ? " gh:" gh_pos : ""
         if($1 in tr_recs)
            addtl = addtl " trf"
         print ">"$1 " QPct:" qcov_pct addtl
         print $2
      }
   ' <(get_trf_recs) \
     <(bawk '{print $name, $seq, $comment}' $candidate_fasta) >$mito_hifi_recs
     # flatten into tabbed fields since the regex {} syntax for gh search not available in bioawk
}

function get_trf_recs {
   echo "# at least one line must be output for awk to see this as a file"

   unset errmsg
   local trf_path=${wdir}/$cand_trf_dir
   local rpt_file=${trf_path}/trf_stats_report.tsv

   [ ! -d $trf_path ] && errmsg="No $trf_path directory"
   [ ! -s $rpt_file ] && errmsg="$rpt_file could not be found"  # empty file indicates no tandem repeats found in the candidate recs

   [ ! -z "$errmsg" ] && msglog_module "$errmsg" && return

   cat $rpt_file
}

function run_trf_on_candidates {
   local trf_path=${wdir}/$cand_trf_dir

   [ -f $trf_path/trf_stats_report.tsv ] && msg "$cand_trf_dir/trf_stats_report.tsv already created" && return

   mkdir_if_needed $trf_path
   pushd $trf_path >/dev/null

   cand_filename=$(basename $candidate_fasta)
   msglog_module "checking candidate records in $cand_filename for CR tandem repeats"

   run_trf.sh ../$cand_filename
   msg ""

   rm *.html
   popd >/dev/null

   # promote records in ${wdir}/${cand_trf_dir}/trf_stats_report.tsv to mito_hifi_recs
}

# using seqrequester from marbl github for this
function basic_rec_stats {

   fasta_for_stats=$out_fasta

   [ ! -s $fasta_for_stats ] && return

   hifi_mito_recs=$(numrecs $fasta_for_stats)
   [ ! -z $hifi_mito_recs ] && update_setting_if_changed "HiFi_mito_reads" $hifi_mito_recs

   set_fastx_basename $fasta_for_stats
   stat_file=${wdir}/${fastx_basename}.seqinf

   [ -s $stat_file ] && msg "$(basename $stat_file) already created" && return
   [ -z $(which seqrequester) ] && return  # stat program not installed or inaccessible

   seqrequester summarize $fasta_for_stats >$stat_file
   echo -e "\n\nsequence lengths:" >>$stat_file
   seqrequester summarize -simple $fasta_for_stats >>$stat_file

   # also write to log and settings the reason the rec was included
   filter_breakdown=$(get_candidate_filter_stats)
   update_setting_if_changed "HiFi_read_filter_stats" "$filter_breakdown"
   msglog_module "$hifi_mito_recs HiFi reads promoted from $(numrecs $candidate_fasta) candidates based on these filters:"
   msglog_module "    $filter_breakdown"
}

# stats for why we included the canidate in the final mito_hifi_rec set
function set_candidate_filter_stat_vars {

   fasta_for_stats=$out_fasta

   local rslts=$(awk -v qcov_min=$qcov_min '
      $1 >= qcov_min { QPct++; next }
      /gh:/ && /trf/ { gh_trf++; next }
      /gh:/ { gh++; next }
      /trf/ { trf++; next}
      END {
         print int(QPct), int(gh_trf), int(gh), int(trf)
   }' <(grep "^>" $fasta_for_stats |sed "s/^.*QPct://" | sort -k1,1nr))

   QPct=$(echo $rslts   | awk '{print $1}')
   gh_trf=$(echo $rslts | awk '{print $2}')
   gh=$(echo $rslts     | awk '{print $3}')
   trf=$(echo $rslts    | awk '{print $4}')
}

function get_candidate_filter_stats {
   set_candidate_filter_stat_vars
   echo "QPct $QPct, gh and trf $gh_trf, gh $gh, trf $trf"
}

function rec_db_exists {
   dbname=$1
   blastdbcmd -db $dbname -info >/dev/null 2>/dev/null
   cmpl_code=$?
   return $cmpl_code
}

function make_mito_hifi_rec_db {
   out_dir=${wdir}/mito_hifi_rec_db
   [ ! -d $out_dir ] && mkdir $out_dir && msglog_module "created $(basename $out_dir) to hold a blast database containing mito candidate fasta records"
   [ ! -d $out_dir ] && msglog_module "directory $out_dir could not be created" && exit 4

   set_fastx_basename $out_fasta
   dbname=${out_dir}/$fastx_basename

   if ! rec_db_exists $dbname; then
      makeblastdb -in $out_fasta -out $dbname -dbtype nucl -parse_seqids
   else
      msg "$fastx_basename database already exists"
   fi
}

function selections_from_candidates {
   # pull the associated fasta records from each file in the tsv_rec file, flipping records as needed so all are in same mito direction
   if [ ! -s $candidate_fasta ]; then
      fofn_pull_recs
   else
      msg "$(basename $candidate_fasta) already created with $(numrecs $candidate_fasta) records"
   fi

   run_trf_on_candidates  # will use this output to possibly include more of the candidate records in the chosen set
   filter_recs_from_candidates  # this will make the mito_hifi_recs.fasta file for downstream use
}

# job of these is to create the ${wdir}/mito_hifi_recs.fasta file
if [ ! -s $out_fasta ]; then

   check_input_files  # make sure we have input necessary
   selections_from_candidates  # this will pull all the candidates recs and from it make the mito_hifi_recs.fasta file for downstream use, aka $out_fasta

else
   msg "$(basename $out_fasta) already created with $(numrecs $out_fasta) records"
fi

basic_rec_stats

# make a blast database from those records in a subdir named mito_hifi_rec_db
make_mito_hifi_rec_db
