#/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

# for each directory in tr_copynum_crs and the cr_recs.fa in it we will run mafft, consensus and make annotation for the consensus
# using mitfi, OH blast, goose hairpin scan, and trf

function show_then_run_cmd {
   local cmd="$@"

   msglog_module "$name: $cmd"
   eval $cmd
}

function blast_OH_db {
   local qry=$1
   local evalue=.00001

   local OH_db=${data_dir}/OrgRepl/OH.fas

   msglog_module blastn -db OH.fas -query $qry -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -max_target_seqs 5 -subject_besthit -evalue $evalue
   blastn -db $OH_db -query $qry -outfmt  "6 std staxid stitle qlen qcovhsp qcovus"  -max_target_seqs 5 -subject_besthit -evalue $evalue
}

function create_cr_tr_copynum_anno {

   local stats=trf_stats_report.tsv
   local mitfi=cr_consensus_${name}.mitfi
   local fa=../cr_consensus_${name}.fasta

   ${script_dir}/add_goosehairpin_to_mitfi.sh $mitfi $fa | \
    ${script_dir}/add_TR_to_mitfi.sh          $stats -   | \
    ${script_dir}/add_OH_to_mitfi.sh          OH.tsv -   | \
    awk '/^#header/{if(NR==1)next; $0="#"} {print}'    > cr_consensus_${name}.anno

   msglog_module "cr_consensus_${name}.anno created with flanking tRNA, OH, tandem repeats, goose hairpin results\n"
}

function process_one_dir {
   local tr_dir=$1
   name=$(basename $tr_dir)  # used in show_then_run_cmd, do not make local

   # need to check to see if the work is already done, that means cr_consensus_${name}.fasta, cr_consensus_${name}.anno and cr_consensus_${name}.compare exist
   unset consensus_made; unset anno_made; unset comp_made; unset m
   [ -s ${tr_dir}/cr_consensus_${name}.fasta ] && m="cr_consensus_${name}.fasta " && consensus_made=yes
   [ -s ${tr_dir}/cr_consensus_${name}.anno ]  && m="${m}cr_consensus_${name}.anno "  && anno_made=yes
   [ -s ${tr_dir}/cr_consensus_${name}.compare ] && m="${m}cr_consensus_${name}.comp " && comp_made=yes

   [ ! -z "$m" ] && msg "${m}already created"

   pushd $tr_dir >/dev/null

   if [ -z $consensus_made ]; then

      show_then_run_cmd "mafft --auto cr_recs.fa 2>cr_recs_mafft.log >cr_recs_mafft.fasta"
      show_then_run_cmd "consensus_from_fasta_alignment.sh cr_recs_mafft.fasta >cr_consensus_${name}.fasta"

   fi

   if [ -z $anno_made ]; then
      mkdir_if_needed anno
      [ ! -d anno ] && msglog_module "error creating $tr_dir/anno" && return 1

      # mitfi for flanking tRNA
      show_then_run_cmd "mitfi.sh -code $code -onlycutoff cr_consensus_${name}.fasta >anno/cr_consensus_${name}.mitfi"

      # blast for OH
      blast_OH_db cr_consensus_${name}.fasta >anno/OH.tsv

      cd anno

      # trf for tandem repeats
      msglog_module "$name: run_trf.sh cr_consensus_${name}.fasta"
      run_trf.sh ../cr_consensus_${name}.fasta
      msg "\n"

      # from these inputs, create an anno file
      create_cr_tr_copynum_anno

      cd ..
   fi

   # make a softlink to the anno file so that the consensus fasta and anno are at same dir level
   [ -s anno/cr_consensus_${name}.anno ] && [ ! -s cr_consensus_${name}.anno ] && ln -s anno/cr_consensus_${name}.anno

   if [ -z $comp_made ]; then
      CR_tr_consensus_path=cr_consensus_${name}.fasta
      CR_tr_consensus_anno=cr_consensus_${name}.anno
      create_comparison  # uses CR_tr_consensus_path and msa_CR_consensus_path file sequences
   fi

   popd >/dev/null
}

#############################
#   compare file creation   #
#############################
function getseq {
   bawk '{print $seq; exit}' $1
}
function get_cr_tr_seq {
   fasta=$1
   anno=$2

   bawk '
      FNR==NR && /tandem_repeat/ && start==0 { start = $5; end = $6 }
      FNR==NR{next}

      {
         if (start > 0 && end > start) {  # convert the tandem repeat portion to lower case
            lngth = end - start + 1
            modstr($seq, start, lngth)
         }
         print $seq
         exit
   } ' <(awk 'BEGIN{print ">dummy line"}/tandem_repeat/{print ">"$0}' $anno) \
       $fasta
}
function get_cr_tr_anno {
   [ ! -s $CR_tr_consensus_anno ] && echo "# no anno" && return  # need to feed it one line so ot looks like a file

   cat $CR_tr_consensus_anno
}

function mitolen {
   bawk '{print length($seq); exit}' $1
}

function prt_results {
   local show_per_line=$1; [ -z $show_per_line ] && show_per_line=100
   local show_all=$2; [ -z $show_all ] && show_all=1
   local name1=$(basename $msa_CR_consensus_path)
   local name2=$(basename $CR_tr_consensus_path)

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
            mtch_count = chrs["|"]
            indel_count = chrs["-"]
            addtl = (indel_count=="") ? "" : " " indel_count " indels"
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
   ' <(echo "# no anno") <(get_cr_tr_anno) -
   # reads from piped in info
}

function create_comparison {
   edlib_str.py $(getseq $msa_CR_consensus_path) $(get_cr_tr_seq $CR_tr_consensus_path $CR_tr_consensus_anno) \
 | sed "s/'//g" | prt_results \
 > cr_consensus_${name}.compare

   [ -s cr_consensus_${name}.compare ] && msglog_file cr_consensus_${name}.compare
}
###############################################################################

###############################################################################
#         set env vars and call the functions for each relevant dir           #
###############################################################################

code=$(get_setting_or_default "code" "5")  # mito code for mitfi, it is usually 5 or 2
msa_CR_consensus_path=$(realpath ${cr_dir}/CR_recs_w_flanks.consensus.fa)

cd $cr_dir

for dir in tr_copynum_crs/*_tr; do
   [ -d $dir ] && [ -s $dir/cr_recs.fa ] && process_one_dir $dir
done
