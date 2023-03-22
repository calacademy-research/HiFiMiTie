#!/bin/bash

# source this file in the CR analysis script in the HiFiMiTie script directory

function function_exists() { declare -F "$1" > /dev/null; }
! function_exists msg && function msg { echo -e "$@" >/dev/stderr; }

function cr_analysis_dirname {
   dsuff=$1
   [ -z "$dsuff" ] && dsuff="seqs"

   CR_dir=CR_range_${dsuff}
}

function get_length_sorted_tsv_sequences { # prepend length as first field of each rec
   bawk '
      {print length($seq), $0}
   ' $CR_input | sort -k1,1nr
}

function make_range_files {
   local ext_bp=2

   cawk -t -v ext_bp=$ext_bp -v Range_dir=$Range_dir '
      NR==1 { range_begin = $1; range_end = $1 }

      within_cur_range() {
         append_to_range(); next
      }

      { output_range(NR==2) }  # always output longest sequence, even if only 1 item in range
      END { output_range(1) }  # always output shortest sequence, even if only 1 item in range

      function output_range(force_it) {
         if (force_it || range_items >= 3) {
            outfile = range_filename()
            for (r=1; r <= range_items; r++ ) {
               print ">" recnames[r] " " comments[r] "\n"seqs[r] > outfile
            }
         }

         delete recnames; delete seqs; delete comments
         range_items = 0

         range_begin = $1
         append_to_range()
      }
      function append_to_range() {
         ++range_items
         recnames[range_items] = $2
         seqs[range_items] = $3
         comments[range_items] = $5
         range_end = $1
      }
      function within_cur_range() {
         return (range_end - $1) <= ext_bp
      }
      function range_filename() {
         infix = "_recs_"
         if (range_begin!=range_end) prefix = "cr_lens_"; else prefix = "cr_len_"
         prefix = Range_dir "/" prefix
         if (range_begin==range_end)
            return prefix range_begin infix range_items  ".fas"
         else
            return prefix range_end "-" range_begin infix range_items ".fas"
      }
   ' <(get_length_sorted_tsv_sequences)
}

function make_msa_from_range_files {
   cr_analysis_dirname msa  # sets CR_dir var
   msa_dir=$CR_dir
   cr_altnum=0  # global used in other functions

   for f in $Range_dir/cr*.fas; do
      local nrecs=$(numrecs $f)
      local bn=$(basename $f)
      local nm=$(echo $bn | sed "s/.fas//")

      if [ $nrecs == 1 ]; then  # no need for msa just copy with a new name for the CR record
         mkdir -p $CR_altdir
         (( cr_altnum++ ))  # needs to be global
         local fa=$CR_altdir/CR_alt${cr_altnum}.fa

         bawk -v cr_altnum=$cr_altnum '{print ">cr_alt" cr_altnum "_len_" length($seq) " 1_rec"; print $seq; exit}' $f | seqfold 120 >$fa
      else # msa of all the recs followed by a consensus call resulting new single rec file in CR_altdir
         mkdir -p $msa_dir
         local msa=$msa_dir/$(replace_ext.sh $bn msa)

         [ ! -s $msa ] && mafft --auto $f > $msa
      fi
   done
}

function make_cr_alts_from_msa_files {
   cr_analysis_dirname msa  # sets CR_dir var
   msa_dir=$CR_dir
   [ ! -d $msa_dir ] && return  # no msa dir, nothing to do

   for msa in $msa_dir/*.msa; do
      mkdir -p $CR_altdir
      (( cr_altnum++ ))  # needs to be global -- was incremented for the 1 rec ones that were copied
      local fa=$CR_altdir/CR_alt${cr_altnum}.fa

      consensus_from_fasta_alignment.sh $msa |
      bawk -v cr_altnum=$cr_altnum '{
         print ">cr_alt" cr_altnum "_len_" length($seq) " " $name
         print $seq
      }' | seqfold 120 > $fa
   done

   msglog_module "$cr_altnum Control Region(s) created from ranges of close CR lengths"
   update_setting_if_changed "CR_range_alternates_created" $cr_altnum
}

# all the trf scripts want the input and output to be in the current dir
# this is trf_dir or trf_path which is in the cr_analysis directory
# cr_analysis is where the other results are collected

function cr_alt_tandem_repeats_analysis {
   cr_analysis_dirname trf  # sets CR_dir var
   trf_dir=$CR_dir
   [ ! -d $CR_altdir ] && return  # no alt CRs found

   msglog_module "Annotate alternate Control Regions for Repeats, Goose Hairpin, and OH"

   mkdir -p $trf_dir
   cat $CR_altdir/CR_alt*.fa > $trf_dir/CR_alts.fasta  # trf needs a file with all the sequences in it to evaluate for repeats

   local trf_report=trf_stats_report.tsv
   trf_report_path=$trf_dir/$trf_report  # trf_report_path will be used later, so do not make it local

#   [ -s trf_dir/CR_trf_overview.txt ] && msg "$trf_report_path already created" && trf_already_run=1 && return

   pushd $trf_dir >/dev/null

   run_trf.sh CR_alts.fasta
   echo "" # no final eol by trf program

   rm *.html
   cawk -t '{print $1;for(f=3;f<NF;f+=5)print fldcat(f,f+4)}' $trf_report |
   awk 'NF>1{sub("--"," ",$1)}{print}' > CR_trf_overview.txt

   # write repeats file for each CR_alt with repeats in the CR_trfdir in anno style format -- add any OH or gh latter
   cawk -v CR_altdir=$CR_altdir '
      NF==1{
         recname = $1; nf = split(recname,ar,"_"); fname = "CR_" ar[2] ".repeats"
         len = ar[nf]
         printf("%s\t%d \t%d\t.\t.\t.\tcr\tControl Region\t+\n", recname, 1, len) > fname
         next
      }
      {
         printf("%s\t%s\t%s\t%s\t.\t.\trepeat\tRepeat Region\t+\t%s copies %s bp consensus: %s\n",
                  recname, $1, $2, $6, $8,$12,$NF) > fname
      }
   ' CR_trf_overview.txt

   popd >/dev/null
}

function finish_anno {  # basic anno file in trfdir .repeats files. add gh and OH if found
   [ ! -d $trf_dir ] && return  # nothing to do

   [ -z ${script_dir} ] && script_dir="$hfdir"
   [ -z ${script_dir} ] && return # was not set and we could not find hfdir set either

   for rpt in $trf_dir/*.repeats; do
      bn=$(basename $rpt .repeats)
      fa=$CR_altdir/$bn.fa
      anno=$CR_altdir/$bn.anno
      if [ -s $fa ]; then
         $script_dir/add_goosehairpin_to_mitfi.sh $rpt $fa |
         add_OH_to_anno $fa > $anno  # todo: look for OH in $fa and add it if found
      fi
   done
}

function add_OH_to_anno {
   if [ -s ${data_dir}/OrgRepl/OH.fas ]; then  # this lets us run the script directly without this happening
      local fasta=$1; local bn=$(basename $fasta .fa)
      local OH_blast_path=$trf_dir/${bn}_OH.tsv

      blast_OH_to_fasta $fasta $OH_blast_path
      ${script_dir}/add_OH_to_mitfi.sh $OH_blast_path -
   else
      cat -
   fi
}
function blast_OH_to_fasta {
   qry=$1
   evalue=.00001
   output_file=$2
   OH_db=${data_dir}/OrgRepl/OH.fas

   msglog_module blastn -db OH.fas -query $qry -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads
   blastn -db $OH_db -query $qry -outfmt  "6 std staxid stitle qlen qcovhsp qcovus"  -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads >$output_file
}

function run_cr_range_analysis {
   CR_input=CR_recs_no_flanks.fasta
   CR_altdir=CR_alternates

   cr_analysis_dirname
   Range_dir=$CR_dir

   [ ! -d $Range_dir ] && mkdir $Range_dir

   make_range_files

   make_msa_from_range_files
   make_cr_alts_from_msa_files
   cr_alt_tandem_repeats_analysis
   finish_anno
}
