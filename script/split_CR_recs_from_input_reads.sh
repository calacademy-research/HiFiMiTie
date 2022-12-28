#!/bin/bash

source $(dirname $(realpath $0))/shared.sh  # sets msg, script_dir is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

maxeditdist=4
cr_seqs_from_input=CR_with_flanks_pulled_from_input.fa

# this creates a subdir in split_seqs dir and given the CR flank names and the
#  split file already made for the CR  take the flanks from this file and for
# each of the 2 make a mafft alignment and then a consensus fasta in the subdir

# use the flank sequence fasta files and read through all of the original input
# reads again and extract the flanks and the sequence (CR) in n=between them

################################################################################
#                              utility functions                               #
################################################################################

function msg_for_CR_search {
   msglog ""
   msglog_module "Control region repeat diversity can cause the mitodb blast to skip some relevant HiFi reads"
   msglog_module "Use CR flanks $cr_begin_flank and $cr_end_flank to search the HiFi reads for these flanks (edit distance $maxeditdist)"
   msglog_module "and pull the flanks and the CR sequence between them into records in a new a fasta file."
   msglog ""
}

# probably move these 2 funcs to shared.sh since 08_split_seqs does type of stuff too

# sets: start_element end_element start_element_3let end_element_3let cr_begin_flank cr_end_flank
function set_CR_flanks_vars {
   local CR=$(get_setting "Primary_CR")
   [ -z $CR ] && msglog_module "Could not find the Primary Control Region setting" && return 1

   local flanks=$(get_setting ${CR}_flanks)
   [ -z "$flanks" ] && msglog_module "Could not find the flanking elements for $CR" && return 1

   # for example we have P F flanking the CR: then F is the starting element and P is the end but P is begin flank and F is end flank
   start_element=$(echo $flanks | awk '{print $2}'); start_element_3let=$(three_letter_AA $start_element)
   end_element=$(echo $flanks | awk '{print $1}');   end_element_3let=$(three_letter_AA $end_element)

   # for backward compatibility when we thought only trna elements could be the flanks
   start_trna_3let=$start_element_3let;  end_trna_3let=$end_element_3let

   cr_begin_flank=$end_element_3let
   cr_end_flank=$start_element_3let

   [ -z $cr_seqs_from_input ] && cr_seqs_from_input=CR_with_flanks_pulled_from_input.fa
   true
}

# this is the split file created from the matching mito reads
function get_file_for_CR_flanks {
   ! set_CR_flanks_vars && return 1

   CR_seqs=$splitseq_dir/${cr_begin_flank}_CR_${cr_end_flank}.fasta
   [ ! -s $CR_seqs ] && msglog_module "Could not find file \"$CR_seqs\" to extract flanks" && return 1
   true
}

function create_cr_seq_search_dir {
   local dir=$splitseq_cr_search_subdir
   [ -d $dir ] && return 0
   mkdir $dir

   [ -d $dir ] && > $dir/__use_cr_flanks_to_search_hifi_inputs_and_pull_seqs__ && return 0

   false
}

function cat_input_files {
   fofn=$wdir/files_to_search.fofn
   msglog_module "Reading files listed in $fofn"

   while IFS= read -r hifi_input_file; do
      [ -z "$hifi_input_file" ] && continue; [ ! -s "$hifi_input_file" ] && continue

      msglog_module $hifi_input_file
      cat $hifi_input_file
   done < $fofn
}

# make a softlink for the generic name of the CR seq file and if it is better than split seq version
# rename split_seq version with first_cut_ prefix and softlink to the XR seq file in the subdir
function make_CR_softlink {  # make a softlink for the generic name of the CR seq file
   [ -z "$CR_seqs" ] && get_file_for_CR_flanks
   [ -z "$CR_seqs" ] && return 1; [ ! -s $CR_seqs ] && return 1
   [ -z "$new_cr_seq_file" ] && return 2; [ ! -s $new_cr_seq_file ] && return 2

   recs_first=$(numrecs $CR_seqs)
   recs_new_search=$(numrecs $new_cr_seq_file)

   if (( recs_new_search > recs_first )); then
      nrecs=$recs_new_search
      orecs=$recs_first
      link_to=$rel_path_new_cr_seq_file

      # add prefix to $CR_seqs
      CR_seqs_base=$(basename $CR_seqs)
      mv $CR_seqs $(dirname $CR_seqs)/first_cut_$CR_seqs_base
      ln -s $link_to $CR_seqs
   else
      nrecs=$recs_first
      orecs=$recs_new_search
      link_to=$(basename $CR_seqs)
   fi

   msglog_module "Creating softlink CR_recs_w_flanks.fasta to $link_to having $nrecs records (from $orecs records)"
   ln -s $link_to $splitseq_dir/CR_recs_w_flanks.fasta
}

################################################################################
#                                main functions                                #
################################################################################

# get a consensus sequence from the flanks in the cm_results CR split sequence file
# it must have the format used in the cm_results split CR file for this to work
# in particular CR seq is uppercase flanks lowercase

function get_begin_flank_seqs {
   bawk '{
      sub("[ACGT].*", "", $seq)
      printf(">begin_flank_%d\n%s\n", ++rec,$seq)
   }' $CR_seqs
}

function get_end_flank_seqs {
   bawk '{
      sub("[acgt]*[ACGT]*", "", $seq)
      printf(">end_flank_%d\n%s\n", ++rec,$seq)
   }' $CR_seqs
}

function get_sequence_consensus {
   local msa=$(replace_ext $1 mafft)
   mafft --auto $1 > $msa

   local consensus=$(replace_ext $1 consensus.fa)
   consensus_from_fasta_alignment.sh $msa >$consensus
}

function make_flank_consensus_files {
   local dir=$splitseq_cr_search_subdir

   get_begin_flank_seqs > $dir/begin_flanks.fa
   get_end_flank_seqs   > $dir/end_flanks.fa

   get_sequence_consensus $dir/begin_flanks.fa > $dir/begin_flanks.consensus.fa
   get_sequence_consensus $dir/end_flanks.fa   > $dir/end_flanks.consensus.fa
}

function run_extract_seq {
   local dir=$splitseq_cr_search_subdir

   # extract_seq_by_flanks.sh <seq_fasta> <seq or file beg seq> <seq or file end seq> [<max edit distance> <begin_flank name> <end_flank name> <seq name>]
   # depends on set_CR_flanks_vars having been called
   cat_input_files | $script_dir/extract_seq_by_flanks.sh -                           \
                          $dir/begin_flanks.consensus.fa $dir/end_flanks.consensus.fa \
                          $maxeditdist                                                \
                          $cr_begin_flank $cr_end_flank CR                            \
       > $new_cr_seq_file
}

function search_for_cr_in_input_reads_using_flanks {
   ! get_file_for_CR_flanks && return 1  # could not get CR split file name

   ! create_cr_seq_search_dir && msglog_module "Problem creating $splitseq_cr_search_subdir for new CR searcu" && return 2

   msg_for_CR_search

   # make a consensus file for the cr_begin_flank and cr_end_flank element
   make_flank_consensus_files

   # using the 2 consensus files for the flank sequences to search in the original hifi reads
   run_extract_seq

   # softlink the generic name for the CR seqs to either first split file or the new one
   make_CR_softlink
}

################################################################################
#                                 starts here                                  #
################################################################################

new_cr_seq_file=$splitseq_cr_search_subdir/$cr_seqs_from_input
rel_path_new_cr_seq_file=$(basename $splitseq_cr_search_subdir)/$cr_seqs_from_input

run_if_no_file   search_for_cr_in_input_reads_using_flanks   $new_cr_seq_file
