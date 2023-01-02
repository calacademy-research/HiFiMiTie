#!/bin/bash

# use the info from mito_hifi_recs.cm_anno and one_line_per_rec.cm_anno.srt
# to pull subseqs from fasta recs in ../mito_hifi_recs.fasta

# using the designations added to one_line_per_rec.cm_anno.srt pull the hiqual template matching sequences
# these just exclude ones that do not match the template or have too few item. keep the wrap around one

# these could be 1 letter or 3 letter AA symbols passed in as args
first_ele=$1
last_ele=$2

# type can be $first_ele, $last_ele or "CR"
type=$3
[ -z $type ] && type=$first_ele

# set file paths based on cm_results dir location
dir=$4
[ ! -d "$dir" ] && dir="."

oneliner=${dir}/one_line_per_rec.cm_anno.srt
cm_anno=${dir}/mito_hifi_recs.cm_anno
recs=${dir}/../mito_hifi_recs.fasta

MINLEN=150

# this lets first_ele and last_ele be either 1 let or 3 let AA syms
function assign_ele_vars {
   beg_ele=$(AAsyms.sh $first_ele | awk '{print $1}')
   beg_3let=$(AAsyms.sh $beg_ele | awk '{print $2}')

   end_ele=$(AAsyms.sh $last_ele | awk '{print $1}')
   end_3let=$(AAsyms.sh $end_ele | awk '{print $2}')
}

function msg {
   echo -e "$@" >&2
}

# use the quality analysis in the the oneliner comments to exlcude some of these
function hiqual_cm_anno {
   awk '
      FNR==NR && /^m/ {
         if (match($0, "# [TN]", ar))
            next

         hiqual[$1] = $0
         if (match($0, "# W"))
            wraps[$1] = $0
      }
      FNR==NR { next }

      $1 in wraps {
         print $0 "\t# W wraps around"
         next
      }
      $1 in hiqual {
         print
      }
   ' $oneliner $cm_anno
}

# typically beg_ele to end of seq but can stop earlier if we have end_ele in the sequence
function subseqs_after_CR {
   awk -v beg_ele=$beg_ele -v end_ele=$end_ele -v beg_3let=$beg_3let -v end_3let=$end_3let -v MINLEN=$MINLEN '

      # cm anno lines per element found. we only care about beg_ele and end_ele
      FNR==NR && $7==beg_ele {
         if (! ($1 in splitseq_start)) # we want the earliest hit of the beg_ele
            splitseq_start[$1] = $2
      }
      FNR==NR && $7==end_ele {
         splitseq_end[$1] = $2 # this could be before or after the start one, this pos in $2 tells us
      }
      FNR==NR{ next }

      # looping through sequence names and sequences in $1 $2
      $1 in splitseq_start {
         start = splitseq_start[$1]
         reclen = length($2)
         seq = $2

         if (! ($1 in splitseq_end) || splitseq_end[$1] < start) { # print from beg_ele to end of sequence
            sublen = reclen - start + 1
            if (sublen > MINLEN) {
               printf(">%s_%s_to_end %d..%d %dnt\n", $1, beg_3let, start, reclen, sublen)
               print substr(seq, start, sublen)
            }
         } else { # print from beg_ele to right before end_ele
            end = splitseq_end[$1] - 1
            sublen = end - start + 1
            if (end > start && sublen > MINLEN) {
               printf(">%s_%s_to_%s %d..%d %dnt\n", $1, beg_3let, end_3let, start, end, sublen)
               print substr(seq, start, sublen)
            }
         }
      }
   ' <(hiqual_cm_anno) <(bawk '{print $name, $seq}' $recs)
}

# from start of seq to end_ele. if the beg_ele is found skip this record
function subseqs_before_CR {
   awk -v beg_ele=$beg_ele -v end_ele=$end_ele -v beg_3let=$beg_3let -v end_3let=$end_3let -v MINLEN=$MINLEN '

      # cm anno lines per element found. we only care about beg_ele and end_ele
      FNR==NR && $7==end_ele {
         if (! ($1 in splitseq_end)) # we want the earliest hit of the end_ele
            splitseq_end[$1] = $3  # end of the element
      }
      FNR==NR && $7==beg_ele {  # if the beg_ele is in the seq, skip this one since it is handled by subseqs_after_CR
         skip_rec[$1] = $2
      }
      FNR==NR{ next }

      # looping through sequence names and sequences in $1 $2
      $1 in splitseq_end {
         end = splitseq_end[$1]
         seq = $2
         sublen = end
         if (sublen > MINLEN) {
            printf(">%s_beg_to_%s 1..%d %dnt\n", $1, end_3let, end, sublen)
            print substr(seq, 1, sublen)
         }
      }
   ' <(hiqual_cm_anno) <(bawk '{print $name, $seq}' $recs)
}

# extract end_ele CR beg_ele. lowercase the elements at begin and end. also show CR len in header (total len minus ele lens)
function CR_with_flanks {
   awk -v beg_ele=$beg_ele -v end_ele=$end_ele -v beg_3let=$beg_3let -v end_3let=$end_3let -v MINLEN=$MINLEN '

      # cm anno lines per element found. we only care about beg_ele and end_ele
      FNR==NR && $7==end_ele {  # end_ele like Pro or Val first
         if (! ($1 in splitseq_beg)) # ignore it if we have seen it before for this rec
            splitseq_beg[$1] = $2
            len_end_ele[$1] = $3 - $2 + 1
      }
      FNR==NR && $7==beg_ele {  # end of beg_ele eg Phe or Met is where we will stop
         # only use this one if we have seen the end_ele for this rec and not seen beg_ele yet
         if ($1 in splitseq_beg && ! ($1 in splitseq_end))
            splitseq_end[$1] = $3 # the ending pos of the beg_ele
            len_beg_ele[$1] = $3 - $2 + 1
      }
      FNR==NR{ next }

     $1 in splitseq_beg && $1 in splitseq_end {
        start = splitseq_beg[$1]; lc_beg_len = len_end_ele[$1]
        end = splitseq_end[$1];   lc_end_len = len_beg_ele[$1]
        seq = $2
        sublen = end - start + 1
        cr_len = sublen - (lc_beg_len + lc_end_len)
        printf(">%s_%s_CR_%s %d..%d %dnt %s %d %s %d CR %d\n", $1, end_3let, beg_3let, start, end, sublen, end_3let, lc_beg_len, beg_3let, lc_end_len, cr_len)

        subseq = substr(seq, start, sublen)
        seq_beg_ele = tolower( substr(subseq, 1, lc_beg_len) )
        seq_end_ele = tolower( substr(subseq, sublen - lc_end_len + 1) )
        seq_cr = toupper( substr(subseq, lc_beg_len+1, sublen-(lc_beg_len+lc_end_len)) )

        printf("%s\n%s\n%s\n", seq_beg_ele, seq_cr, seq_end_ele)
     }
   ' <(hiqual_cm_anno) <(bawk '{print $name, $seq}' $recs)
}

function check_args {
   if [[ -z $last_ele || "$beg_3let" = "Unk" || "$end_3let" = "Unk" ]]; then
      msg "\n    usage: split_seqs_using_cm_results.sh <first_trna> <last_trna> <type of subseqs to retrieve> [<dir for anno files>]\n\n    type can by the first or last trna or CR\n    dir defaults to pwd\n"
      exit 1
   fi
}

# this lets first_ele and last_ele be either 1 let or 3 let AA syms
function assign_ele_vars {
   beg_ele=$(AAsyms.sh $first_ele | awk '{print $1}')
   beg_3let=$(AAsyms.sh $beg_ele | awk '{print $2}')

   end_ele=$(AAsyms.sh $last_ele | awk '{print $1}')
   end_3let=$(AAsyms.sh $end_ele | awk '{print $2}')
}

##################################################################################################
###                           Get whichever subseqs caller asked for                           ###
##################################################################################################

assign_ele_vars
check_args

if [ "$type" = "$first_ele" ]; then
   subseqs_after_CR
elif [ "$type" = "$last_ele" ]; then
   subseqs_before_CR
elif [ "$type" = "CR" ]; then
   CR_with_flanks
else
   msg "Invalid split seq type requested: \"$type\""
fi
