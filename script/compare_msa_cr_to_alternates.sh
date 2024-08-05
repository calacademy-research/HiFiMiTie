#!/bin/bash

# input: cr1 fasta sequence, cr1 anno, cr2 fasta sequence, cr2 anno
# anno lines that are repeats will have corresponding subsequence lowercased

function compare_cr_seqs {
   [ ! -s "$4" ] && return

   cr1_seq=$1; local cr1_anno=$2
   cr2_seq=$3; local cr2_anno=$4

   edlib_str.py $(prep_seq $cr1_seq $cr1_anno) $(prep_seq $cr2_seq $cr2_anno) | sed "s/'//g" | format_comparison
}

# input: cr fasta sequence, cr anno file -- if cr anno missing just output the sequence
# output: just the sequence with repeats as defined in the anno file lowercased

function prep_seq {
   local cr_seq=$1; local cr_anno=$2

   [ ! -s $cr_anno ] && bawk '{print $seq}' $cr_seq && return

   cawk -t '
      FILENUM==1 { sequence = $1 }
      FILENUM==2 && tolower($7) ~ "repeat" {  # lowercase repeat in sequence
         modstr(sequence, $2, $3-$2+1)
      }
      END { print sequence }
   ' <(bawk '{print $seq}' $cr_seq) $cr_anno
}

function format_comparison {
   cat -
   return

   local show_per_line=$1; [ -z $show_per_line ] && show_per_line=100
   local show_all=$2; [ -z $show_all ] && show_all=1
   local name1=$(basename $cr1_seq)
   local name2=$(basename $cr2_seq)

   bioawk_cas -v show=$show_per_line -v show_all=$show_all -v name1=$name1 -v name2=$name2 '   # we use charcount() from bioawk_cas

      function min(a,b){if(a<=b)return a; return b}
      function dashes(segment,  adash) {
         charcount(segment, adash)
         return int( adash["-"] )
      }
      function update_seq_pos(l, show) { # to keep track of the sequence pos for each we have to discount the dashes in the sequence
         top_seg = substr(ar[1], l, show)
         top_dels = dashes(top_seg)

         bot_seg = substr(ar[3], l, show)
         bot_dels = dashes(bot_seg)

         top_line_beg = top_line_end + 1  # one past wherever we were last time
         top_line_end = top_line_beg + length(top_seg)-1 - top_dels

         bot_line_beg = bot_line_end + 1
         bot_line_end = bot_line_beg + length(bot_seg)-1 - bot_dels

         anno_info()
      }
      function anno_info(   t,tstart,b,bstart) { # update top and bot anno_info based on line and anno coords
         top_anno_str = "";bot_anno_str = ""
         for (t=1; t <= top_anno; t++) {
            tstart = top_start[t]
            if (tstart >= top_line_beg && tstart <= top_line_end) {
               if (top_anno_str=="") top_anno_str = " :"
               top_anno_str = top_anno_str " " top_name[t] " " tstart ".." top_end[t]
            }
         }
         for (b=1; b <= bot_anno; b++) {
            bstart = bot_start[b]
            if (bstart >= bot_line_beg && bstart <= bot_line_end) {
               if (bot_anno_str=="") bot_anno_str = " :"
               bot_anno_str = bot_anno_str " " bot_name[b] " " bstart ".." bot_end[b]
            }
         }
      }

      BEGIN {
         if(!show_all) print "Areas of difference between " name1 " and " name2
         else print "Comparison between " name1 " and " name2
         if (show < 10) show = 100
         top_anno = 0; bot_anno = 0
      }

      FNR==1{ FNum++ }
      FNum < 3 && /^#/ { next }

      FNum==1 {
         top_name[++top_anno] = $7; top_start[top_anno] = $2; top_end[top_anno] = $3
         next
      }
      FNum==2 {
         bot_name[++bot_anno] = $7; bot_start[bot_anno] = $2; bot_end[bot_anno] = $3
         next
      }

      # piped in info is third file
      FNR==1 { gsub("[{}]",""); print $0 "\n" }  # this is the edit distance info

      FNR > 1 {
         ar[FNR-1]=$0
      }

      END {
         len = length(ar[1])
         for (l = 1; l <= len; l += show) {
            mtch_str = substr(ar[2], l, show)
            charcount(mtch_str, chrs)
            mtch_count = int( chrs["|"] )
            indel_count = int( chrs["-"] )
            addtl = (indel_count==0) ? "" : " " indel_count " indels"
            update_seq_pos(l, show)

            to_show = min(l+show-1, len)
            if (show_all || mtch_count < length(mtch_str)) {
               print substr(ar[1], l, show), name1, top_line_beg"-"top_line_end top_anno_str
               print mtch_str, mtch_count " matches" addtl
               print substr(ar[3], l, show), name2, bot_line_beg"-"bot_line_end bot_anno_str
               print ""
            }
         }
      }
   ' $cr1_anno $cr2_anno -
   # reads from piped in info
}

[ -s $4 ] && compare_cr_seqs $@
