#!/bin/bash

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

function get_file_for_flanks {
   CR_seqs=$(ls ../*_CR_*fasta | head -n 1)

   [ ! -s $CR_seqs ] && exit 1
}

function replace_ext {
   local name=$1
   local ext=$2
   local dir=$(dirname $name)

   # show the dir component if there is a slash in the input name
   echo $name | grep "/" >/dev/null && echo -n ${dir}/

   local bname=$(basename $name)
   if echo $bname | grep "\." >/dev/null; then # there is a dot in the name defining the extension
      echo $bname | sed -E "s/(.*)\..*/\1.${ext}/"
   else # no extension append the new one with a dot then the ext
      echo ${bname}.${ext}
   fi
}

get_file_for_flanks
get_begin_flank_seqs >begin_flanks.fa
get_end_flank_seqs   >end_flanks.fa

get_sequence_consensus begin_flanks.fa
get_sequence_consensus end_flanks.fa
