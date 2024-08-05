#!/bin/bash

function prep_seq {
   local cr_seq=$1; local cr_anno=$2

   cawk -t '
      FILENUM==1 { name = $1; sequence = $2; comment = $4 }
      FILENUM==2 && tolower($7) ~ "repeat" {  # lowercase repeat in sequence
         modstr(sequence, $2, $3-$2+1)
      }
      END { print ">" name " " comment; print sequence }
   ' <(bawk '{print}' $cr_seq) $cr_anno
}

function get_msa_CR_anno {
   local cm_mitfi_dir=$(echo $(pwd) | sed -E "s|hfmt_([^/]*/).*|hfmt_\1msa_assembly/cm_mitfi/|")
   local msa_anno=$cm_mitfi_dir/mito_msa.anno

   # pull out repeat lines and subtract beginning CR pos from their start end locs
   [ -s $msa_anno ] && cawk -t ' $8 == "Control Region" { ofs = $2-1 }
                                 ofs { $2-=ofs; $3-=ofs; print }
                               ' $msa_anno
}

function make_all_fasta_file {
   # prep msa CR with anno
   prep_seq $orig_dir/CR_consensus_no_flanks.fasta <(get_msa_CR_anno)

   # get all the alternate CR files with repeat lowercased
   for f in $(ls -S CR*fa); do
      anno=$(replace_ext.sh $f anno)
      prep_seq $f $anno
   done
}

function output_comparison {
   width=160

   if [ "$1" == "stacked" ]; then
      bawk '{print ""}{print ">" $name " " $comment; print $seq}'  | seqfold $width
   else # interleaved
      bawk -v width=$width '
      {
         line = 0; bases_shown = 0
         for (i = 1; i <= length($seq); i += width) {
            slice = substr($seq, i, width); tmp = slice  # count with gsub modifies string
            bases_in_slice = gsub(/[a-zA-Z]/, "A", tmp)
            pos_str = (bases_in_slice) ? sprintf("%d..%d", bases_shown+1, bases_shown + bases_in_slice) : "--..--"
            bases_shown += bases_in_slice

            text = sprintf("%-"width"s  %s %s %s\n", slice, pos_str, $name, $comment)

            line++
            ar[line] = ar[line] text
         }
      }
      END {
         printf("Msa and %d alternate Control Regions compared. Repeats are in lowercase.\n\n", NR-1)  # -1 since msa is in here too

         for (l = 1; l <= line; l++) {
            print ar[l]
         }
      }'
   fi
}

function run_comparison_in_dir {
   unset pushed
   orig_dir=$(pwd)
   [ ! -z $1 ] && [ -d $1 ] && pushd $1 >/dev/null && pushed=pushed

   make_all_fasta_file           |
   mafft --auto --preservecase - |
   output_comparison interleaved > msa_CR_and_alternate_CR.comparisons

  [ ! -z $pushed ] && popd >/dev/null
}


###############################################################################
# if no argument then safe to source and use the functions

[ ! -z $1 ] && [ -d $1 ] && run_comparison_in_dir $1
