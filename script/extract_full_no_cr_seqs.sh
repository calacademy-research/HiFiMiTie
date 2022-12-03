#!/bin/bash

: ' in some cases the HiFi reads are long enough to incorporate the entirety of
    the mito sequence excluding the control region

    this script takes the one line file (usually the sorted version), the cm table and the mito fasta reads
    then determines from the one line which reads are candidates for having the entire non-cr sequence
    and use the cm_anno table to get the beginning and ending positions of that in the read.

    then given the read name and the begin and end position we can pull out the non-cr sequence from the fasta file
    this output is suitable for using mafft to flush out errors if there are enough reads, that is handled elsewhere

    if no $5 then shows a full_no_cr seq anno format
    if $5 is the /mito_hifi_recs.fasta file then the sequences re extracted that reprsent the full_no_cr portion of the relevant reads
'

one_liner=$1
cm_table=$2

cr_succ=$3  # typical vertebrate case M for Methionine
cr_prev=$4  # typical vertebrate case P for Proline

mito_seqs_fasta=$5

cr_succ_3let=$(awk '{print $2}' <(AAsyms.sh $cr_succ))
cr_prev_3let=$(awk '{print $2}' <(AAsyms.sh $cr_prev))

# [ -z $5 ] && usage

function grep_full_noncr_lines {
   # grep "^m.*M.*V" one_line_per_rec.cm_anno.srt | grep -v "# T" | sed s/~//

   qry="^m.*${cr_succ}.*${cr_prev}"
   grep $qry $one_liner | grep -v "# T" | sed s/~//  # the "# T" excludes those that do not match template
}

function get_prev_succ_anno_lines {
   awk -v cr_succ=$cr_succ -v cr_prev=$cr_prev '
      FNR==NR {
         match($0, cr_succ ".*" cr_prev)
         mtch_template = substr($0, RSTART, RLENGTH)
         ar[$1] = $1 "\t" mtch_template
         next
      }

      $1 in ar && ($7 == cr_succ || $7 == cr_prev) {
         if ( ! (seen[$1]) ) {
            print ar[$1]
            begpos = 0; endpos = 0
         }

         if ($7 == cr_succ && begpos==0) {
            print
            begpos = $2
         } else if ($7 == cr_prev && begpos > 0) {
            endpos = $3
            print
            printf("%s\t%d\t%d\tfull_no_cr_seq\tlen: %d\n", $1, begpos, endpos, endpos-begpos+1)
         }
         seen[$1]++
      }
   ' <(grep_full_noncr_lines) $cm_table
}

# >m64049_220702_011313/179112789/ccs_RC_Met_to_end 1808..5439 3632nt
# m64049_220702_011313/30934586/ccs       29      14043   full_no_cr_seq  len: 14015
function get_mito_seqs {
   if [ ! -s "$mito_seqs_fasta" ]; then
      cat
   else
      grep full_no_cr_seq |
      cawk -t -v succ=$cr_succ_3let -v prev=$cr_prev_3let '
         FNR==NR {
            read_begpos[$1] = $2
            read_endpos[$1] = $3
            next
         }
         $1 in read_begpos {
            beg = read_begpos[$1]; end = read_endpos[$1]
            len = end - beg + 1
            printf(">%s_%s_thru_%s_%d-%d %dnt %s\n", $1, succ, prev, beg, end, len, $4)
            print substr($2, beg, len)
         }
      ' - <(bawk '{print}' $mito_seqs_fasta)
   fi
}

get_prev_succ_anno_lines | get_mito_seqs
