#!/bin/bash

: '
m64044_210611_022728/29623447/ccs     -         rrnS                 -          cm        1      899     1829     2803      +    no    1 0.45  16.6  710.0    1e-181 !

m64044_210611_022728/29623447/ccs       1762    1829    62.69   1.098E-12       GAA     F       Metazoa_F.cm    +
m64044_210611_022728/29623447/ccs       1829    2803    710.0   1e-181  complete        !       12S_rRNA.cm     +
'

rrna=$1
[ -z $rrna ] && rrna=rrna_rrnS.tbl

mitfi=$2
[ -z $mitfi ] && mitfi=mito_hifi_recs.mitfi

descrip=$(echo $rrna | awk '{print ($1 ~ "L") ? "16S_rRNA.cm" : "12S_rRNA.cm"}')

awk -v descrip=$descrip '
   FNR==NR && /^#/{next}
   FNR==NR {
      recid = $1
      cmpl = ($11=="no") ? "complete" : $11 "_trunc"
      S_row = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $8, $9, $15, $16, cmpl, $3, descrip, $10)

      S_ar[recid] = S_row
      S_start[recid] = int($8)
      S_end[recid] = int($9)

      next
   }

   { recid = $1 }

   recid != lst_recid && lst_recid in S_ar && S_start[lst_recid] > 0 { #  it is at end of the sequence, likely truncated 3 prime
      print S_ar[lst_recid]
      S_start[lst_recid] = 0
   }

   { lst_recid = recid }

   recid in S_ar { # see if we print the Srna line before this mitfi line
      rec_start = $2; rec_end = $3
      Sbegin = S_start[recid]
      if (Sbegin > 0 && rec_start > Sbegin) {
         print S_ar[recid]
         S_start[recid] = 0
      }
   }

   { print }
' $rrna $mitfi
