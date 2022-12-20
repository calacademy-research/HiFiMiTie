#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

# 1_trnF  m64044_201205_132125/179309382/ccs      100.000 68      0       0       1       68      8790    8857    1.82e-30        126     1536898 N/A     68      100     100

: ' this is one where we would pull to the end since the Thr is before the Pro
Thr_NC_007975   m64044_210611_022728/1639385/ccs_RC     91.429  70      4       2       1       70      4176    4243    9.11e-21        93.3    N/A     QPct:97 70      100     100
Pro_NC_007975   m64044_210611_022728/1639385/ccs_RC     100.000 69      0       0       1       69      5673    5741    1.51e-30        125     N/A     QPct:97 69      100     100
ND6_NC_007975   m64044_210611_022728/1639385/ccs_RC     81.118  519     98      0       4       522     5756    6274    7.92e-141       495     N/A     QPct:97 522     99      99
'

# extract sequences from each record starting at the start_trna (which is trna flanking the end of the CR. F or Phe often, P in birds with remnant CR)
# stop at the end of the sequence or before the end_trna (which is the trna flanking the CR beginning) whichever comes first

fasta=mito_hifi_recs.fasta
feature_matches=blast_results/mito_feature_match_to_cand_recs.tsv

MINLEN=250

function create_dirs {
   [ ! -d $splitseq_dir ] && mkdir $splitseq_dir && msglog_module "created $(basename $splitseq_dir) to hold sequences for assembly by alignment"
   [ ! -d $splitseq_dir ] && msglog_module "directory $splitseq_dir could not be created" && return 5
   [ ! -d $splitseq_feat_subdir ] && mkdir $splitseq_feat_subdir && msglog_module "  $(basename $splitseq_feat_subdir) subdir created"
   [ ! -d $splitseq_cm_results_subdir ] && mkdir $splitseq_cm_results_subdir && msglog_module "  $(basename $splitseq_cm_results_subdir) subdir created"
   true
}

function set_trnas_to_pull_from { # set the trna vars and sequence split msg vars

   CR=$(get_setting "Primary_CR")
   [ -z $CR ] && msglog_module "Could not find the Primary Control Region setting" && return 1

   flanks=$(get_setting ${CR}_flanks)
   [ -z "$flanks" ] && msglog_module "Could not find the flanking tRNAs for $CR" && return 1

   # for example we have P F flanking the CR: then F is the starting trna and P is the end
   start_trna=$(echo $flanks | awk '{print $2}')
   end_trna=$(echo $flanks | awk '{print $1}')
   start_trna_3let=$(three_letter_AA $start_trna)
   end_trna_3let=$(three_letter_AA $end_trna)

   msglog_module "tRNAs $flanks flank the Control Region."

   seq_split_msg1="Splitting sequences from $start_trna_3let to sequence end."
   seq_split_msg2="Splitting sequences from sequence beginning to $end_trna_3let."
   seq_split_msg3="Splitting sequences from $end_trna_3let to $start_trna_3let to capture the Control Region and its flanks"

   # these will be the names of the files in the blast feature subdir (we will change dir portion for other subdir)
   outdir=$splitseq_feat_subdir
   outfile1=${outdir}/${start_trna_3let}_to_end.fasta
   outfile2=${outdir}/beg_to_${end_trna_3let}.fasta
   outfile3=${outdir}/${end_trna_3let}_CR_${start_trna_3let}.fasta

   outfile4=${splitseq_dir}/full_no_cr_seq.fasta

   true
}

# 04Sep2022 at some point we should replace using the blst features results for sequence splitting info
# and move to using the results from cm_results but here will we use the cm_results to  throwout records
# that have problems if there are still enough to work with (that test to be done later)
function filter_reads {
   fasta=$1
   excludes="$splitseq_dir/excluded_from_check.txt"

   local one_liner=$cm_dir/one_line_per_rec.cm_anno.srt

   if [ ! -s "$one_liner" ]; then
      cat $fasta  # no file to use for the sieve just move things along
   else # remember those that have a comment line with a {TWN] aftern the space and we will not pass those along
      bawk -v excludes=$excludes '
         FNR==NR && /# [TWN]/ {
            dont_use[$name] = $0
         }
         FNR==NR { next }

         $name in dont_use {
            print "skipping " dont_use[$name] > excludes
            next
         }
         {
            print ">" $name " " $comment
            print $seq
         }
      ' <(prefix_gt $one_liner) $fasta
   fi
}

function extract_seq_after_CR {

   msglog_module $seq_split_msg1

   bawk -v st=$start_trna_3let -v en=$end_trna_3let -v MINLEN=$MINLEN '
      BEGIN{start_trna_3let="^"st; end_trna_3let="^"en}
      FNR==NR { recname=$4 }
      FNR==NR && $1 ~ start_trna_3let {
         pos_trna_start = $11
         seq_start_pos[recname] = pos_trna_start
      }
      FNR==NR && $1 ~ end_trna_3let {
         pos_trna_end = $11
         if ( (recname in seq_start_pos) && pos_trna_end > seq_start_pos[recname])
         {
            seq_end_pos[recname] = pos_trna_end
         }
      }
      FNR==NR{next}

      $name in seq_start_pos {
         slen=length($seq)
         pos = seq_start_pos[$name]
         if ($name in seq_end_pos) {
            endpos = seq_end_pos[$name]
            outlen = endpos-pos+1
            output = substr($seq, pos, outlen)
            endstr=en
         }
         else {
            output = substr($seq, pos)
            endstr="end"
         }
         outlen = length(output)

         if (outlen > MINLEN) {
            printf(">%s_%s_to_%s %d..%d %dnt\n%s\n", $name, st, endstr, pos, slen, outlen, output)
         }
   }' <(prefix_gt ${wdir_path}/$feature_matches) <(filter_reads ${wdir_path}/$fasta) > $outfile1
}

function extract_seq_before_CR {
   msglog_module $seq_split_msg2

   bawk -v st=$start_trna_3let -v en=$end_trna_3let -v MINLEN=$MINLEN '
      BEGIN{start_mtch="^"st; end_mtch="^"en}
      FNR==NR { recname=$4 }
      FNR==NR && (recname in seq_end_pos || recname in skip_rec) { # eiher already found an earlier version of the trna
         next
      }
      FNR==NR && $1 ~ end_mtch { # this is the trna before the CR and the first time we have seen it for the rec
         pos_trna_end = $12
         seq_end_pos[recname] = pos_trna_end
      }
      FNR==NR && $1 ~ start_mtch { # we have run into the trna that comes after the CR, if we run into it before the other trna, skip this rec
         if (recname in seq_end_pos)
            skip_rec[recname]++
      }

      $name in seq_end_pos {
         endpos = seq_end_pos[$name]
         output = substr($seq, 1, endpos)
         outlen = length(output)

         if (outlen > MINLEN) {
            printf(">%s_beg_to_%s %d..%d %dnt\n%s\n", $name, en, 1, outlen, outlen, output)
         }
   }' <(prefix_gt ${wdir_path}/$feature_matches) <(filter_reads ${wdir_path}/$fasta) > $outfile2
}

function extract_CR_and_flanks {
   tstart=$end_trna_3let; tend=$start_trna_3let # we flip the trnas for the CR extraction
   msglog_module $seq_split_msg3

   bawk -v st=$tstart -v en=$tend -v MINLEN=$MINLEN '
      BEGIN{ start_mtch="^"st; end_mtch="^"en }
      FNR==NR { recname=$4 }
      FNR==NR && $1 ~ start_mtch {  # ^Pro before the typical CR placement (often ^Thr for Aves)
         pos_trna_start = $11
         if (recname in seq_start_pos) # ignore it if we have seen it before for this hifi rec
            next
         seq_start_pos[recname] = pos_trna_start
         seq_start_tlen[recname] = $12 - pos_trna_start + 1
      }
      FNR==NR && $1 ~ end_mtch {    # ^Phe after the typical CR placement (often ^Pro here for Aves)
         pos_trna_end = $12
         if (recname in seq_start_pos) { # only insert the ending position if we have already seen the start for this hifi rec
            if (! (recname in seq_end_pos)) { # only add the first time we see this one
               seq_end_pos[recname] = pos_trna_end
               seq_end_len[recname] = $12 - $11 + 1
            }
         }
      }
      FNR==NR{ next }

      $name in seq_start_pos && $name in seq_end_pos {
         begpos = seq_start_pos[$name]
         endpos = seq_end_pos[$name]
         len = endpos-begpos+1
         output = substr($seq, begpos, len)
         outlen = length(output) # should be same as len
         if (outlen > MINLEN) {
            #modstr(output, 1, seq_start_tlen[$name]) # lowercase the trna bases
            #end_trna_start = seq_end_pos[$name] - seq_end_len[$name] + 1
            #modstr(output, end_trna_start, seq_end_len[$name]) # lowercase the ending trna
            printf(">%s_%s_CR_%s %d..%d %dnt\n%s\n", $name, st, en, begpos, endpos, outlen, output)
         }
      }
   '  <(prefix_gt ${wdir_path}/$feature_matches) <(filter_reads ${wdir_path}/$fasta) > $outfile3
}

# 04Sep2022 this will be a backup if the typical extractions are not working. depends on number of recs we get
function extract_any_full_no_nr_seqs {
   [ ! -s $cm_dir/full_no_cr_seq.info ] && return  # nothing to do, no read had the full_no_cr portion in it

   local one_liner=$cm_dir/one_line_per_rec.cm_anno.srt
   local cm_anno=$cm_dir/mito_hifi_recs.cm_anno
   local recs=${wdir_path}/$fasta

   [ ! -s $one_liner ] && msglog_module "Could not extract full_no_cr_seqs.fasta: $one_liner not found" && return 1
   [ ! -s $cm_anno ] && msglog_module "Could not extract full_no_cr_seqs.fasta: $cm_anno not found" && return 2
   [ ! -s $recs ] && msglog_module "Could not extract full_no_cr_seqs.fasta: $recs not found" && return 3

   $script_dir/extract_full_no_cr_seqs.sh $one_liner $cm_anno $start_trna $end_trna  $recs >$outfile4

   if [ -s $outfile4 ]; then
      nr=$(numrecs $outfile4)
      msglog_module "$(basename $outfile4) written with $nr records"
   fi
}

# replace_dir moved to shared.sh

# 29Nov2022 use the cm_anno results to determine how to split the sequences
function create_cm_extract_file {
   local file_to_create=$1
   local type=$2

   cm_file=$(replace_dir $file_to_create $splitseq_cm_results_subdir)
   $script_dir/split_seqs_using_cm_results.sh $start_trna $end_trna $type $cm_dir >$cm_file
}

# 29Nov2022 make split files in a subdir named cm_results using cm_anno file instead of feature blasts
function make_cm_results_extracts {
   cm_donefile=$(replace_dir $donefile $splitseq_cm_results_subdir)

   if [ ! -f $cm_donefile ]; then
      [ -z "$start_trna" ] && set_trnas_to_pull_from

      create_cm_extract_file $outfile1 $start_trna
      create_cm_extract_file $outfile2 $end_trna
      create_cm_extract_file $outfile3 CR

      numrecs ${splitseq_cm_results_subdir}/*.fasta > $cm_donefile
   fi
}

function make_blast_feature_extracts {
   run_if_no_file extract_seq_after_CR  $outfile1
   run_if_no_file extract_seq_before_CR $outfile2
   run_if_no_file extract_CR_and_flanks $outfile3

   feat_donefile=$(replace_dir $donefile $splitseq_feat_subdir)
   numrecs ${splitseq_feat_subdir}/*.fasta > $feat_donefile
}

# softlink to either the blast features files or to the cm_results fasta depending on some criteria (currently hardwired to cm_results)
function make_softlinks {
   pushd $splitseq_dir >/dev/null
   local source_dir=$(basename $splitseq_cm_results_subdir)

   for fa in $source_dir/*.fasta; do
      ln -s $fa
   done

   popd >/dev/null
}

function set_done {
   non_empty_fastas=$((head -1 ${splitseq_dir}/*.fasta | wc -l) 2>/dev/null)
   [ ! -z $non_empty_fastas ] && [ "$non_empty_fastas" -ge 1 ] && numrecs ${splitseq_dir}/*.fasta > $donefile

   [ -s $donefile ] && msglog_module "Split sequences created for alignment assembly and for Control Region assembly & repeat analysis"
}


##########################################################################
#                      main routine to call:                             #
#              creates 2 sets of the split sequences in subdirs          #
#                   softlink to one of the sets                          #
##########################################################################

function make_split_sequence_files {
   if create_dirs && set_trnas_to_pull_from; then
      # 29Nov2022 going to do 2 sets of splits one using the blast of the features and new one using cm_results and split_seqs_using_cm_results.sh script

      make_blast_feature_extracts
      make_cm_results_extracts
      make_softlinks

      run_if_no_file extract_any_full_no_nr_seqs $outfile4

      # 10Dec2022 create CR flank consensus seqs and use them to pull CR from the original input recs
      $script_dir/split_CR_recs_from_input_reads.sh

      set_done # set a done semaphore file when finished creating the extracts
   fi
}

donefile="${splitseq_dir}/sequence_split_files.done"
run_if_no_file make_split_sequence_files $donefile
