#!/bin/bash

# sometimes using mafft to merge the two pieces can be a problem. especially if the consensus for the ending portion encomapasses the beginning trna

# so this module will use mitfi to identify the trnas in the both sequences.
# the last entry in the beginning consensus sequence's mitfi will be searched for in the ending sequence mitfi to get the ending position of it
# then we will append everthing after that to the beginning sequence to create: non_cr_consensus.fasta

################################################################################
#                              utility functions                               #
################################################################################

function usage {
   msg "\n    usage: merge_two_consensus_fa.sh <begin seq file> <end seq file>\n"
   exit 1
}

# if for HiFiMiTe can delete this one
function msg {
   echo -e "$@" >/dev/stderr
}

function get_fld {
   line="$1"
   fld=$2  # field number
   awk -v fld=$fld '{print $fld; exit}' <(echo $line)
}

################################################################################
#                                main functions                                #
################################################################################

function create_mitfis {
   beg_mitfi=$(echo $begseq | sed "s/\.[a-zA-Z]*$/.mitfi/")
   end_mitfi=$(echo $endseq | sed "s/\.[a-zA-Z]*$/.mitfi/")

   [ ! -s $beg_mitfi ] && msg "$begseq" && mitfi.sh -code $code -onlycutoff $begseq > $beg_mitfi
   [ ! -s $end_mitfi ] && msg "$endseq" && mitfi.sh -code $code -onlycutoff $endseq > $end_mitfi

}

function set_last_beg_entry_info {
   last_entry=$(tail -n 1 $beg_mitfi)

   last_entry_cm=$(get_fld "$last_entry" 8)
   last_entry_type=$(get_fld "$last_entry" 7)
   last_entry_beg=$(get_fld "$last_entry" 2)
   last_entry_end=$(get_fld "$last_entry" 3)
   last_entry_len=$(( $last_entry_end - $last_entry_beg + 1 ))

}

function set_endseq_info {
   begseq_len=$(bawk '{print length($seq); exit}' $begseq)
   endseq_len=$(bawk '{print length($seq); exit}' $endseq)

   last_end_mitfi_match=$(grep -e $last_entry_cm -e $last_entry_type $end_mitfi | tail -n 1)
   last_match_endpos=$(get_fld "$last_end_mitfi_match" 3)

   endseq_append_from=$(( $last_match_endpos + 1 ))

   merge_seq_len=$(( last_entry_end + $endseq_len - $endseq_append_from + 1 ))
}

function write_merged_fa {
   # write header for fasta rec from our various variables
   begseq_basename=$(basename $begseq)
   endseq_basename=$(basename $endseq)
   echo ">non_cr_consensus $begseq_basename 1.."$last_entry_end $endseq_basename $endseq_append_from".."$endseq_len ${merge_seq_len}nt

   # write begseq sequence to end of last entry (used -onlycutoff so we can trust it to be complete) and from endseq_append_from to end of the endseq
   bawk -v begpos=$endseq_append_from -v begseq_maxlen=$last_entry_end '
      NR==1 {
         prefix = substr($seq, 1, begseq_maxlen)
         print prefix
      }

      FILENUM > 1 {
         suffix=substr($seq, begpos)
         print suffix
         exit
      }
   ' $begseq $endseq
}

################################################################################
#                                 starts here                                  #
################################################################################

## set input vars ##

begseq=$1
endseq=$2

[ ! -s $begseq ] && usage; [ ! -s $endseq ] && usage

code=$3; [ -z $code ] && code=5

## do the work ##

create_mitfis
set_last_beg_entry_info
set_endseq_info
write_merged_fa
