#!/bin/bash

# using the mifi file, show the counts of the trnas to the right of a trna in each of the mito recs
# change to have a minimum for most_freq to use it and also ignore X

mitfi=$1
[ ! -s "$mitfi" ] && mitfi=mito_hifi_recs.mitfi

show_first=$2
[ -z "$show_first" ] && show_first="F"

awk -v show_first=$show_first '
   BEGIN{OFS="\t"; show_first = toupper(show_first); sp = " "}
   function most_freq_neighbor(t) {
      biggest_freq = 0; most_freq = ""
      if (t in right_neighbor) {
         for (n in right_neighbor[t]) {
            if (right_neighbor[t][n] > biggest_freq) {
               biggest_freq = right_neighbor[t][n]
               most_freq = n
            }
         }
      }
      return most_freq
   }
   function prt_order() {
      for(i=1; i<=o; i++)
         printf("%s ", order[i])
      printf("\n")
   }
   function insert_val_before_this_in_order(val, rt_neighbor) {
      num_items=length(order)
      n = 1
      for(i=1; i<=num_items; i++) { # copy items before we see rt_neighbor into new array
         if(order[i] != rt_neighbor)
           nord[n++] = order[i]
         else
            break
      }
      if(order[i] == rt_neighbor) # insert val here
         nord[n++] = val
      for(; i<=num_items; i++) {
         nord[n++] = order[i]
      }
      delete order

      nord_len = length(nord)
      for(n=1; n <= nord_len; n++)
         order[n] = nord[n]
   }

   # skip comment and low quality hits
   { skip = /^#/ || ! ($5 ~ "E") || $7 == "X" }

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
   $8 ~ "^OL" {
      skip = 0
      $7 = "OL"
   }
   $7 == "OH" {
      skip = 0
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
         if (mf == show_first || mf in seen) {
            break
         }
         order[++o] = mf
         seen[mf]++
         delete trna[mf]
      }
      for (t in trna) { # any we did not hit naturally, see what their neighbor would be
         # order[++o] = t
         mf = most_freq_neighbor(t)
         if(mf!="" && biggest_freq >= 8) { # insert it into the order before the mf value
            if(mf == show_first) # append to end
               order[++o] = t
            else # insert t before mf
               insert_val_before_this_in_order(t, mf)
            o = length(order)
         }
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
