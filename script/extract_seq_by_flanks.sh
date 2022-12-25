#!/bin/bash

# created to look for CR sequences by the rna flanks of the CR tho could be used elsewhere

DEFAULTMAX=3
debug_info=0

################################################################################
#                              utility functions                               #
################################################################################

function usage {
  msg "
    usage: extract_seq_by_flanks.sh <seq_fasta> <seq or file beg seq> <seq or file end seq> [<max edit distance> <begin_flank name> <end_flank name> <seq name>]

    first arg is the fasta seqs to search. can also be fastq and can be gzipped

    args 2 and 3 are either fasta files or sequence strings
        if a files then we use the first record of the file for the sequence

    optional 4th arg is max edit distance, defaults to $DEFAULTMAX

    optional 5th through 7th can be used to replace generic begin_flank, end_flank, Seq
"
   exit 1
}

function msg {
   echo -e "$@" >&2
}

# validate and get arguments

function set_args { # call with $@
   [ -z $3 ] && usage

   reads_to_check=$1
   [ ! $1 = "-" ] && [ ! -s $reads_to_check ] && usage

   begin_flank=$2
   [ -s $begin_flank ] && begin_flank=$(bawk '{print $seq; exit}' $begin_flank)

   end_flank=$3
   [ -s $end_flank ] && end_flank=$(bawk '{print $seq; exit}' $end_flank)

   maxdist=$4
   [ -z $maxdist ] && maxdist=$DEFAULTMAX
   ! is_int $maxdist && msg "\n    max edit distance must be integer. It is \"$maxdist\"" && usage

   # instead of generic names can use these args to use the specific ones being used
   bf_name=$5; ef_name=$6; seq_name=$7
}

# sets int_value if argument passed to function is an int
# input can have commas like 1,234 or underscores like 100_404 and end in . or .0
function is_int {
   putative_int=$1
   re='^[+-]?[0-9]+$'

   local val_to_check=$(echo $putative_int | sed -e "s/,//g" -e "s/_//g" -e "s/\.0*$//")
   unset int_value

   if [[ $val_to_check =~ $re ]]; then
      int_value=$(awk -v val_to_check=$val_to_check 'BEGIN{print int(val_to_check)}')
      true
   else
      false
   fi
}

################################################################################
#                                main functions                                #
################################################################################

# search for begin and end flanks excerpting as reads them and seq between when approproate
# example matchInfo: 0 70= 11519 11588
# example split seq header: >m64044_201205_132125/23724902/ccs_RC_Pro_CR_Phe 89..4740 4652nt Pro 70 Phe 68 CR 4514

function do_search_and_extract {
   bawk -v begin_flank=$begin_flank -v end_flank=$end_flank -v maxdist=$maxdist -v debug_info=$debug_info -v bf_name="$bf_name" -v ef_name="$ef_name" -v seq_name="$seq_name" '
      BEGIN {
             beg_len = length(begin_flank); end_len = length(end_flank)
             begin_flank = toupper(begin_flank); end_flank = toupper(end_flank)
             lc_begin_flank = tolower(begin_flank); lc_end_flank = tolower(end_flank)

             if(bf_name=="") bf_name = "beginflank"; if(ef_name=="") ef_name="endflank"; if(seq_name=="") seq_name="Seq"
             flank_cmnt = sprintf("_%s_%s_%s", bf_name, seq_name, ef_name) # "_begflank_CR_endflank"

             mode = 22
             start_secs = systime()
      }

      (NR % 5000)==0 {
         show_progress()  # progress indicator (do not use the CR anymore)
      }

      {
         if ( ! check_seq("FWD", $seq) )
            check_seq("REV", revcomp($seq))
      }

      END { show_progress(); printf("\n") >"/dev/stderr" }

      function check_seq(prefix, sequence) {
         slen = length(sequence)
         matchInfo = edit_dist(maxdist, begin_flank, beg_len, sequence, slen, mode)
         if (matchInfo == -1) { return 0 }

         split(matchInfo, ar, " ");
         dist = ar[1]; B_begpos = ar[3]; B_endpos = ar[4]; B_cigar = ar[2]

         # found close enough begin flank, look for end flank
         matchInfo = edit_dist(maxdist, end_flank, end_len, sequence, slen, mode)
         if (matchInfo == -1) { return 0 }

         split(matchInfo, ar, " ");
         dist = ar[1]; E_begpos = ar[3]; E_endpos = ar[4]; E_cigar = ar[2]

         # make sure end flank is after beg
         ext_beg = B_endpos+1; ext_end = E_begpos-1
         ext_len = ext_end - ext_beg + 1
         if (ext_len > 50) { # at least 50 bp between them {
            seqs_found++
            if (need_cr) { printf "\n" > "/dev/stderr"; need_cr = 0 }
            if (debug_info) {
               printf("%s rec %d (%s): first flank %d-%d, %s ", prefix, NR, slen, B_begpos, B_endpos, B_cigar)
               printf("end flank %d-%d, %s\tseq btw %d\t%s\n", E_begpos, E_endpos, E_cigar, ext_len, $name)
            }
            else {  # write extraction record example hdr: m64044_201205_132125/23724902/ccs_RC_Pro_CR_Phe 89..4740 4652nt Pro 70 Phe 68 CR 4514
               addtl = (prefix=="REV") ? "_RC" : ""
               printf(">%s%s%s %d..%d %dnt of %d %s %d\n", $name, addtl, flank_cmnt, B_begpos, E_endpos, E_endpos-B_begpos+1, slen, seq_name, ext_len)
               print lc_begin_flank
               print substr(sequence, ext_beg, ext_len)
               print lc_end_flank
            }
            return 1
         }
         return 0
      }

      function show_progress() {
         secs_so_far = systime() - start_secs; secs_per_thousand = secs_so_far / (NR/1000)
         printf "\r%d records searched %d sequences found. %s (%ss per thousand)   ", NR, seqs_found, sec_fmt(secs_so_far), fstr(secs_per_thousand) >"/dev/stderr"
         need_cr = 0
      }
      function fstr(numstr) {
         fnumstr=sprintf("%f", numstr)
         sub("0+$", "", fnumstr) # get rid of trailing zeroes
         sub("\\.$", "", fnumstr) # get rid of decimal point if nothing after it
         return fnumstr
      }
      function sec_fmt(secs) {
         m = ""; s = ""; h = ""
         if (secs >= 3600) {
            h = int(secs/3600) "h"
            secs -= h * 3600
         }
         if (secs >= 60) {
            m = int(secs/60) "m"
            secs -= m * 60
         }
         if (secs > 0)
            s = secs "s"
         return h m s
      }

   ' $reads_to_check
}

################################################################################
#                                 starts here                                  #
################################################################################

set_args $@

do_search_and_extract
