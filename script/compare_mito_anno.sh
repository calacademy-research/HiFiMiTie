#!/bin/bash

# given 2 mito anno file report on which looks better
# this will evaluate the trna rrna and genes but not the Control Region. that is doen separately

: '
Anno1 better
Anno2 better
Equal for trna rrna genes. CR not evaluated.
Issues with both. Choosing Anno2
'

anno_1=$1; anno_2=$2

if [[ ! -s $anno_1 || ! -s $anno_2 ]]; then
   echo -e "\n    usage: compare_mito_anno.sh <anno file 1> <anno file 2> \n" >&2 && exit 1
fi

awk -v fname1=$(basename $anno_1) -v fname2=$(basename $anno_2) '
   BEGIN { FS="\t"; OFS="\t"; diff_total = 0 }
   FNR == 1 { FNum++ }

   /^#/{ next }
   $7 == "OL" { next } # ignore OL until we search for it in megahit as we do for msa

   { ele = $7; score = int($4) }
   score == 0 { next }

   ! (ele in seen){ ele_list[++e] = ele; seen[ele] = e}

   { ele_score[ele][FNum] = score }

   END {
      same_list = "Same_score:"
      anno1_Missing = "Anno1_Missing:"; anno2_Missing = "Anno2_Missing:"

      for (i = 1; i <= e; i++) {
         ele = ele_list[i]
         anno1_score = ele_score[ele][1]
         anno2_score = ele_score[ele][2]

         if (anno1_score == 0 || anno2_score == 0) {  # presume one of them is missing the ele
            if (anno1_score == 0) {
               anno1_Missing = anno1_Missing " " ele
               anno1_missing_count++
            } else {
               anno2_Missing = anno2_Missing " " ele
               anno2_missing_count++
            }
         } else if (anno1_score == anno2_score) { # scores for the shared ele are the same
            same_score++
            same_list = same_list " " ele
         } else { # one of them is better than the other
            diff_lines[++d] = ele " " (anno1_score - anno2_score) " " fname1 ": " anno1_score " " fname2 ": " anno2_score
            diff_total += anno1_score - anno2_score  # positive favors anno1, negative favors anno2
         }
      }

      if (same_score == length(ele_list)) {
         print "Equal for trna rrna genes. CR not evaluated.\t" fname2 "\n"
         print same_list
         exit
      }

      call_favored(diff_total, anno1_missing_count, anno2_missing_count)

      if (anno1_missing_count > 0 || anno2_missing_count > 0) {
          if (anno1_missing_count > 0) {
             print anno1_Missing
          }
          if (anno2_missing_count > 0) {
             print anno2_Missing
          }
      }
      print "\nScore difference between same features -- positive favors anno1, negative favors anno2: " diff_total
      if (length(diff_lines) > 0) {
         for (i = 1; i <= d; i++)
            print "    " diff_lines[i]
      }
      print ""
      print same_list
   }

   function call_favored(diff_total, anno1_missing_count, anno2_missing_count) {
      fb = "Anno1 better"; sb = "Anno2 better"
      probs_favor1 = "Issues with both. Choosing Anno1"; probs_favor2 = "Issues with both. Choosing Anno2"
      call = sb
      if (anno1_missing_count==0 && anno2_missing_count==0) {  # none missing go by diff_total
         if (diff_total > 0) call = fb; else call = sb
      } else {  # some missing
         missing = anno2_missing_count - anno1_missing_count # positive favors anno1, negative favors anno2
         if (missing == 0) {  # same number of missing, go with diff_total, 0 favoring anno2
            if (diff_total > 0) call = fb; else call = sb
         } else if (missing < 0) {  # more missing in anno1, either anno2 has none missing or fewer missing
            if (anno2_missing_count > 0) call = probs_favor2; else call = sb
         } else {  # more missing in anno2, either anno1 has none missing or fewer missing
            if (anno1_missing_count > 0) call = probs_favor1; else call = fb
         }
      }
      if (call == fb || call == probs_favor1)
         fname = fname1
      else
         fname = fname2

      print call "\t" fname
   }
' $anno_1 $anno_2
