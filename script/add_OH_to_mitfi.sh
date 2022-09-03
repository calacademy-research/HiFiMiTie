#!/bin/bash

# adds an OH blast record found in $wdir/blast_results/OH_blast_to_cand_recs.tsv
# this file is passed in as $1

OH_blast_file=$1
mitfi=$2
min_score_for_OH=250  # kludgy way to throw out the OH look alikes hovering over 12S, we'll try for better method later hopefully

# [ ! -s $OH_blast_file ] && exit 0
# [ ! -s $mitfi ] && echo "No mitfi file found for inserting OH lines" && exit 1

function cat_OH_file { # if OH file is empty output one comment line so we do not confuse mitfi with file onne
   if [ -s $OH_blast_file ]; then
      cat $OH_blast_file
   else
      echo "# no OH blast results"
   fi
}

awk -v min_score_for_OH=$min_score_for_OH '
   BEGIN{FS="\t"; OFS="\t"}
   FNR==NR {
      score = $12

      if ($1 ~ "^#" || $1 in OH_line || score < min_score_for_OH)
         next

      sub(":.*","",$2)
      OH_line[$1] = sprintf("%s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $7, $8, $12, $11, ".", "OH", "blastn", "+")
      OH_startpos[$1] = $7  # we will set this to 0 when we output it so we only print it once
      next
   }

   # process mitfi
   { recid = $1 }  # could be "#" but that is OK

   recid != lst_recid && lst_recid in OH_line && OH_startpos[lst_recid] > 0 { #  it is at end of the sequence
      print OH_line[lst_recid]
      OH_startpos[lst_recid] = 0  # mark it so we do not print again (tho unnecessary here, still mark it)
   }

   { lst_recid = recid }

   recid in OH_line { # see if we print the Srna line before this mitfi line
      rec_start = $2; rec_end = $3
      Sbegin = OH_startpos[recid]
      if (Sbegin > 0 && rec_start > Sbegin) {
         print OH_line[recid]
         OH_startpos[recid] = 0
      }
   }

   { print }

' <(cat_OH_file) $mitfi
