#!/bin/bash

: '
m64044_201011_075919/459824/ccs         ND4L T P ND6 CYTB S1 ND1 rrnL Q I L2 Y ND2 W M N S1 E COX1 C COX2 K D ATP8 ATP6 COX3 G ND3 A OH L1 V rrnS F R ND
'

function make_right_neighbor_list {
   local outfile=$1

   awk -v first=$first '
      {
         for(f=2; f < NF; f++) {
            rtnbr[$f][$(f+1)]++
         }
      }
      END {
         for (f in rtnbr) {
            printf("%s:\t", f)
            PROCINFO["sorted_in"] = "@val_num_desc"
            for (n in rtnbr[f]) {
               printf("%s %s\t", n, rtnbr[f][n])
               if (! (f in chosen_nbr))
                  chosen_nbr[f] = n
            }
            printf("\n")
         }

         # new get template
         cur = first; printf("%s", cur); shown[cur]++
         while ( cur in chosen_nbr ) {
            nbr = chosen_nbr[cur]
            if (nbr in shown) break
            printf(" %s", nbr)
            shown[nbr]++
            cur = nbr
         }
         printf("\n")
      }
   ' one_liners.txt > $outfile
}

function flip_template_to_top {
   local file=$1

   lns=$(wc -l $file | awk '{print $1}')
   ((lns--))

   # flip template line to top of output
   cat <(tail -n 1 $file) <(echo) <(head -n $lns $file)
}

first=$1
[ -z $first ] && first="L1"

tmp=tmp_${RANDOM}

make_right_neighbor_list $tmp
flip_template_to_top $tmp

rm $tmp
