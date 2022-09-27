#!/bin/bash

# modify to handle OL.cm as well as the 2 rrna

# 16Sep2022s .tbl file might be empty so put a comment as part of the file strean so we have something to count as the first file

: '
m64044_210611_022728/29623447/ccs     -         rrnS                 -          cm        1      899     1829     2803      +    no    1 0.45  16.6  710.0    1e-181 !

m64044_210611_022728/29623447/ccs       1762    1829    62.69   1.098E-12       GAA     F       Metazoa_F.cm    +
m64044_210611_022728/29623447/ccs       1829    2803    710.0   1e-181  complete        !       12S_rRNA.cm     +
'

rrna=$1
[ -z $rrna ] && rrna=rrna_rrnS.tbl

mitfi=$2
[ -z $mitfi ] && mitfi=mito_hifi_recs.mitfi

descrip=$(echo $rrna | awk '{descr = ($1 ~ "L") ? (($1 ~ "OL.tbl") ? "OL.cm" : "16S_rRNA.cm") : "12S_rRNA.cm"; print descr}')

awk -v descrip=$descrip '
   FNR==NR && /^#/{next}
   FNR==NR && !($1 in S_ar) {
      recid = $1
      cmpl = ($11=="no") ? "complete" : $11 "_trunc"
      s = int($8); e = int($9); if (s > e){t=s; s=e; e=t} # flip them

      S_row = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, s, e, $15, $16, cmpl, $3, descrip, $10)

      S_ar[recid] = S_row
      S_start[recid] = s
      S_end[recid]   = e
   }

   FNR==NR {
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

   END { # 04Sep2022 handle case where it belongs after the last line
      if (recid in S_ar && S_start[recid] > 0) {
         print S_ar[recid]
      }
   }
' <(echo "#so file is not empty" && cat $rrna) $mitfi
