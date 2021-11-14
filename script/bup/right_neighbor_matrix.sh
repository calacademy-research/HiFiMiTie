#!/bin/bash

# using the mifi file, show the counts of the trnas to the right of a trna in each of the mito recs

mitfi=$1
[ ! -s "$mitfi" ] && mitfi=mito_rec_candidates.mitfi

show_first=$2
[ -z "$show_first" ] && show_first="F"

awk -v show_first=$show_first '
   BEGIN{OFS="\t"; show_first = toupper(show_first); sp = " "}
   function most_freq_neighbor(t) {
      biggest_freq = 0; most_freq = ""
      for (n in right_neighbor[t]) {
         if (right_neighbor[t][n] > biggest_freq) {
            biggest_freq = right_neighbor[t][n]
            most_freq = n
         }
      }
      return most_freq
   }

   # skip comment and low quality hits
   { skip = /^#/ || ! ($5 ~ "E") }

   # handle goosehairpin and 12S, 16S rrnas
   $8 == "goose_hairpin" {
      skip = 0
      $7 = "gh"
   }
   $8 ~ "^12S" {
      skip = 0
      $7 = "12S"
   }
   $8 ~ "^16S" {
      skip = 0
      $7 = "16S"
   }

   # skip comment and low quality hits
   skip {next}

   $1 != lst { cur=$7; lst=$1; next }

   {
      right_neighbor[cur][$7]++
      cur = $7
      trna[cur]++
   }

   END {
      # fill array in the order seen most often
      order[++o] = show_first; delete trna[show_first]
      for(cur=show_first; cur in right_neighbor; cur = mf) {
         mf = most_freq_neighbor(cur)
         if (mf == show_first) {
            break
         }
         order[++o] = mf
         delete trna[mf]
      }
      for (t in trna) { # add any we did not hit naturally
         order[++o] = t
         delete trna[t]
      }

      # print header line in the order seen most often
      for (t1=1; t1 <= o; t1++) {
         printf("\t%s", order[t1])
      }
      printf("\n")

      # now for each one show the counts
      for (t1=1; t1 <= o; t1++) {
         cur = order[t1]
         printf ("%s", cur)
         for (t2=1; t2 <= o; t2++) {
            neigh = order[t2]
            val = int(right_neighbor[cur][neigh])
            val = val==0 ? "." : val
            printf("\t%s", val)
         }
         printf("\n")
      }
   }
' $mitfi
