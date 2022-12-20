#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

function get_editdist_result {
   bawk '
      { nm[NR]=$name; lngth[NR]=length($seq); sq[NR]=$seq }
      END {
         if(NR!=2) print "expected 2 sequences to compare, got " NR >"/dev/stderr"

         rslt = edit_dist(-1, sq[1], lngth[1], sq[2], lngth[2], 20)
         split(rslt, parts, " ")  # FS is tab so need to specify space as delim
         printf("%s %s\n", parts[1], parts[2])
    }' <(cat $msa_path $mega_path)
}

function getfld() {
   local space_delim_str="$1"
   local fld_to_retrieve="$2"
   awk -v space_delim_str="$space_delim_str" -v fld_to_retrieve="$fld_to_retrieve" '
      BEGIN {
         split(space_delim_str, ar, " ")
         print ar[fld_to_retrieve]
      }
   '
}

function getseq {
   bawk '{print $seq; exit}' $1
}

function mitolen {
   bawk '{print length($seq); exit}' $1
}

# mito    1       68      62.69   1.938E-12       GAA     F       Metazoa_F.cm    +
function get_msa_anno {
   local msa_anno=${wdir_path}/mito_msa.cm_anno
   [ -s ${wdir_path}/mito_msa.anno ] && msa_anno=${wdir_path}/mito_msa.anno  # 02Dec2022 use the anno with genes if available

   if [ -s $msa_anno ]; then
      cat $msa_anno
   else # placeholder comment so file has at least one line, otherwise not counted as a file by awk
      echo "# no anno file found "
   fi
}
function get_megahit_anno {
   local megahit_anno=${wdir_path}/mito_megahit.cm_anno
   [ -s ${wdir_path}/mito_megahit.anno ] && megahit_anno=${wdir_path}/mito_megahit.anno  # 02Dec2022 use the anno with genes if available

   if [ -s $megahit_anno ]; then
      cat $megahit_anno
   else # placeholder comment so file has at least one line, otherwise not counted as a file by awk
      echo "# no anno file found "
   fi
}

function prt_results {
   local show_per_line=$1; [ -z $show_per_line ] && show_per_line=160
   local show_all=$2; [ -z $show_all ] && show_all=0

   bioawk_cas -v show=$show_per_line -v show_all=$show_all '   # we use charcount() from bioawk_cas

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
         if(!show_all) print "Areas of difference between mito_megahit.fasta and mito_msa.fasta\n"
         else print "Complete sequence comparison between mito_megahit.fasta and mito_msa.fasta\n"
         if (show < 10) show = 160
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
      FNR > 1 {
         ar[FNR-1]=$0
      }

      END {
         len = length(ar[1])
         for (l = 1; l <= len; l += show) {
            mtch_str = substr(ar[2], l, show)
            charcount(mtch_str, chrs)
            mtch_count = chrs["|"]
            indel_count = chrs["-"]
            addtl = (indel_count=="") ? "" : " " indel_count " indels"
            update_seq_pos(l, show)

            to_show = min(l+show-1, len)
            if (show_all || mtch_count < length(mtch_str)) {
               print substr(ar[1], l, show), "mito_megahit.fasta", top_line_beg "-" top_line_end top_anno_str
               print mtch_str, mtch_count " matches" addtl
               print substr(ar[3], l, show), "mito_msa.fasta    ", bot_line_beg "-" bot_line_end bot_anno_str
               print ""
            }
         }
      }
   ' <(get_megahit_anno) <(get_msa_anno) -
   # reads from piped in info
}

function set_anno_vars {
   megahit_anno=${wdir_path}/mito_megahit.cm_anno
   [ -s ${wdir_path}/mito_megahit.anno ] && megahit_anno=${wdir_path}/mito_megahit.anno  # 02Dec2022 use the anno with genes if available

   msa_anno=${wdir_path}/mito_msa.cm_anno
   [ -s ${wdir_path}/mito_msa.anno ] && msa_anno=${wdir_path}/mito_msa.anno  # 02Dec2022 use the anno with genes if available
}

# use the anno info to tell which elements are being shown as different in the overview
function anno_overview {
   # megahit_mitochondrion	14067	14836	280.8	4.3e-72	complete	rrnS	12S_rRNA.cm	-
   set_anno_vars  # set megahit_anno msa_anno
   rm -f $elediffs  # do not want this file if no differing elments

   bioawk_cas -t -v elediffs=$elediffs '
      FILENUM < 3 && /^#/ { next }
      FILENUM == 1 {  # megahit anno
         mh_ele_start[$7] = $2; mh_ele_end[$7] = $3
         long_name[$7] = $8
      }
      FILENUM == 2 {  # msa anno
         msa_ele_start[$7] = $2; msa_ele_end[$7] = $3
      }
      FILENUM != 3 { next }

      /megahit.fasta/ {  # print a line showing the element name that has differences between the 2 assemblies
         cmnt = ""; tmplt = "%s:  mito_megahit.fasta %d-%d  mito_msa.fasta %d-%d"
         if (split($NF, pos, "-") == 2) {  # mito_megahit.fasta 15041-15200
            for (e in mh_ele_start) {
               to_show = (length(e)>2) ? e : long_name[e]
               if ( pos_in_ele(e, pos[1]) || pos_in_ele(e, pos[2]) )
                  cmnt = " " cmnt sprintf(tmplt, to_show, mh_ele_start[e], mh_ele_end[e], msa_ele_start[e], msa_ele_end[e])
            }
         }
         cmnt = substr(cmnt, 2) # skip the initial space
         if (cmnt != "" && cmnt != lst_cmnt) {
            print cmnt
            print cmnt >> elediffs
         }
         lst_cmnt = cmnt
      }

      { print }

      function pos_in_ele(ele, ofs) {
         return ofs >= mh_ele_start[ele] && ofs <= mh_ele_end[ele]
      }
   ' $megahit_anno $msa_anno -
}

# how to grep the triple digit or greater equal runs the -E lets us use the + for one or more
# echo 14962=9D1=7D1=1D2=3D2=1D3=2D1=11D1=5D1321=1I856= | grep -E "[0-9][0-9][0-9]+="
function show_basic_results {

   local cigar_cats=$(echo $cigar | awk '{print gensub("([0-9][0-9][0-9]+=)", "\n\\1\n", "G")}' | awk '!/^$/')
   local bp_in_long_runs=$(echo -e "$cigar_cats" | awk '/^[0-9]+=$/{l+=int($1)}END{print l}')

   local start_matches=$(echo $cigar | grep -E "^[0-9]+=" -o | sed "s/=//")
   [ -z $start_matches ] && start_matches=0
   local end_matches=$(echo $cigar | grep -E "[0-9]+=$" -o | sed "s/=//")
   [ -z $end_matches ] && end_matches=0
   local flank_matches=$(($start_matches+$end_matches))

   echo -e "Comparison between $mega ($mega_len bp) and $msa ($msa_len bp):\n"

   echo edit distance $eddist
   echo $bp_in_long_runs bp in runs of matches 100 or greater
   echo -e $start_matches matches at the beginning and $end_matches at end for $flank_matches contiguous matches"\n"

   echo Following  describes how to transform the $mega_len bp $mega into the $msa_len bp $msa
   echo -e "$cigar_cats\n"
}

function compare_annos {
   set_anno_vars  # set megahit_anno msa_anno

   echo -e "Comparison of the annotation files $(basename $megahit_anno) and $(basename $msa_anno):\n"
   ${script_dir}/compare_mito_anno.sh $megahit_anno $msa_anno
   echo -e "\n"

   # get the first line from the compare script and store info in settings.sh for the Complete analysis to use
   local result=$(${script_dir}/compare_mito_anno.sh $megahit_anno $msa_anno | head -n 1)
   local anno_chosen=$(echo $result | awk '{print $NF; exit}')
   local assembly_chosen=$(replace_ext $anno_chosen fasta)

   update_setting "chosen_anno" "$anno_chosen"
   update_setting "chosen_assembly" "$assembly_chosen"
   update_setting "anno_compare_info" "$result"
}

function update_elediff_setting {
   diffval=$(awk 'BEGIN{FS=":"}NR>1{printf("\t")}{printf("%s", $1)}' $elediffs)
   update_setting "asm_differing_elements" "$diffval"
}

function create_compare_dir {
  [ ! -d $compare_dir ] && mkdir $compare_dir && msglog_module "$(basename $compare_dir) created" && return 0
  [ ! -d $compare_dir ] && msglog_module "Problem creating $compare_dir"
}

#############################################
#      main function for comparison         #
#############################################

function compare_mega_msa_assemblies {
   # ready to start comparison mark it in log
   msglog_module "Comparing $mega ($mega_len bp) and $msa ($msa_len bp)"

   # 02Dec2022 add comparison based on the anno files for megahit and msa
   compare_annos >$overview

   # result has spaces so make sure to quote
   rslt="$(get_editdist_result)"

   # used in show_basic_results and stored in settings
   eddist=$(getfld "$rslt" 1)
   cigar=$(getfld "$rslt" 2)

   show_basic_results >>$overview

   # 03Sep2022 reversing order to mega then msa, since it was backwards from what prt_results expected
   edlib_str.py $(getseq $mega_path) $(getseq $msa_path) | prt_results $show_per_line $show_all | anno_overview >>$overview

   show_all=1
   edlib_str.py $(getseq $mega_path) $(getseq $msa_path) | prt_results $show_per_line $show_all >$fullseq

   msglog ""
   msglog_file $overview

   update_elediff_setting  # 05Dec2022
   update_setting "msa_megahit_edit_dist" $eddist
   update_setting "msa_megahit_edit_cigar" "$cigar"
}

###################################################################################################
#                    comparison of megahit and msa_consensus assemblies                           #
###################################################################################################

mega=mito_megahit.fasta
mega_path=${wdir}/$mega

msa=mito_msa.fasta
msa_path=${wdir}/$msa

# make sure the assemblies are there to compare
[ ! -s $mega_path ] && msglog_module "Could not find $mega for comparison with $msa" && exit 2
[ ! -s $msa_path ]  && msglog_module "Could not find $msa for comparison with $mega" && exit 3

mega_len=$(mitolen $mega_path)
msa_len=$(mitolen $msa_path)

show_all=0
show_per_line=160

#files to hold comparison info
overview=$compare_dir/comparison_overview.txt
fullseq=$compare_dir/fullsequence_comparison.txt
elediffs=$compare_dir/differing_elements.txt

# create a directory to hold the comparison files
create_compare_dir

run_if_no_file compare_mega_msa_assemblies $overview
