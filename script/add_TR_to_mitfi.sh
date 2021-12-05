#!/bin/bash

# add tandem repeat info from the trf_stats_report.tsv file in $1 into the mitfi file of $2
# consensus_of_19_recs	962	1030	61.65	1.899E-13	GAA	F	Metazoa_F.cm	+
# consensus_of_19_recs	1030	681	723 pos	44 len	61 score	2.0 copies	Consensus pattern 22 bp: AATGCTTGATAGACATTAAGTT

TR_file=$1
mitfi=$2

[ ! -f $TR_file ] && exit 0

function process_TR_file {  # run sed to change 681--723 to 681\t723
   if [ -s $TR_file ]; then
      sed "s/--/\t/" $TR_file
   else # empty file, meaning no repeats. but must put something out as a file or our FNR==NR test will not work
      echo "# dummy line"
   fi
}

awk '
   BEGIN{FS="\t"; OFS="\t"; dot="."}
   FNR==NR {
      if ($1 ~ "^#") next
      if ($1 in TR_line) next

      start = int($3); end = int($4); score=int($6); copies=$7; patt=$8
      evalue = dot

      TR_line[$1] = sprintf("%s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s", $1, start, end, score, evalue, dot, "tandem_repeat", "trf", "+", copies, patt)
      TR_startpos[$1] = start  # we will set this to 0 when we output it so we only print it once
      next
   }

   # process mitfi
   { recid = $1 }  # could be "#" but that is OK

   recid != lst_recid && lst_recid in TR_line && TR_startpos[lst_recid] > 0 { #  it is at end of the sequence
      print TR_line[lst_recid]
      TR_startpos[lst_recid] = 0  # mark it so we do not print again (tho unnecessary here, still mark it)
   }

   { lst_recid = recid }

   recid in TR_line { # see if we print the Srna line before this mitfi line
      rec_start = $2; rec_end = $3
      Sbegin = TR_startpos[recid]
      if (Sbegin > 0 && rec_start > Sbegin) {
         print TR_line[recid]
         TR_startpos[recid] = 0
      }
   }

   { print }

   END {
      if (TR_startpos[$1] > 0) { # tandem repeat is at end
         print TR_line[$1]
         TR_startpos[$1] = 0  # mark it so we do not print again (tho unnecessary here, still mark it)
      }
} ' <(process_TR_file) $mitfi
