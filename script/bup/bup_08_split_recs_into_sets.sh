#!/bin/bash

# depending on the mito structure we divide by 2 or more sets
# more than 2 when more than one Control Region is found

# this had been hardwired to use trnF and trnP as the cut points
# however the beginning trna though F in the verebrates could be
# it will be M in Lepidoptera etc.
# also the last trna is often trnE in birds
# so we have put these into the settings file as first_trna and last_trna
# also for multiple control regions we might need another set and we have the trna_order setting as well

: '
split_recs_into_sets.sh: from mito_hifi_recs.fasta using mito_feature_match_to_cand_recs.tsv info, creates:
   recs_linear_split_nocr.fa
   recs_cr_w_trna_flanks.fa
'

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

: ' here are the 2 types of layouts we will process. first type is rarer of the 2, with it we take everything from start of Phe to end of Pro -- any CR portion is after Pro and is not copied
    second type we can get up to 3 recs. from rec begin to Pro end, from Phe begin to rec end. And for CR file from Pro begin to Phe end
Phe_NC_033356   4196521/ccs_RC  94.203  69      4       0       1       69      777     845     Pro_NC_033356   4196521/ccs_RC  94.118  68      3       1       1       68      16146   16212   dist: 15368 trnas: 137 CR: 15231
Pro_NC_033356   4392637/ccs_RC  92.754  69      3       2       1       68      1914    1981    Phe_NC_033356   4392637/ccs_RC  94.203  69      4       0       1       69      2884    2952    dist: 972 trnas: 137 CR: 835
'

fasta=$1
recinfo_tsv=$2

trna_span="trnF_trnP" # todo get this info from the settings
usage="    usage: split_recs_into_sets.sh <mito_hifi_recs.fasta> <${trna_span}_distances.tsv>\n   output: recs_linear_split_nocr.fa recs_cr_w_trna_flanks.fa"
[ -z "$recinfo_tsv" ] && msg -e "\n$usage\n" && exit 1

# make a split_sequences directory in which to store the split fasta files


function log_split_msg {
   splitmsg_file=${wdir}/split_msg
   if [ ! -f $splitmsg_file ]; then
      msglog_module "Problem with reporting split files stats"
   else
      while IFS= read -r line; do
         msglog_module $line
      done < $splitmsg_file
      rm -f  $splitmsg_file
   fi
}

cd $wdir_path # so we write the fa files in the working dir

awk 'BEGIN {
      FS="\t"; OFS="\t"
      min_recsize = 500
      split_recs  = "recs_linear_split_nocr.fa"
      cr_recs     = "recs_cr_w_trna_flanks.fa"
   }

   NR==1 && (has_mach_id=($2 ~ "^m")){}
   FNR==NR && /^Phe/{ type_ones++
      typ[$2] = 1
      Phe_start[$2] = $9; Pro_end[$2] = $20  # we will take from Phe start to Pro end for the record, excluding any portion of CR after Pro
   }
   FNR==NR && /^Pro/{ type_twos++
      typ[$2] = 2
      Pro_end[$2] = $10 # one record is from rec begin to Pro end (if long enough)
      Phe_start[$2] = $19 # another record is from Phe start to end of the rec (if long enough)

      Pro_start[$2] = $9; Phe_end[$2] = $20 # record containing CR with flanking tRNAs Pro and Phe
   }
   FNR==NR{next}

   {rname = $1}
   ! has_mach_id{ sub("m[0-9_]*/", "", rname) }

   typ[rname]==1 {
      start = Phe_start[rname]; end = Pro_end[rname]; len = end - start + 1
      printf(">%s_phe-pro Phe %d Pro %d len %d\n%s\n", $1, start,end,len, substr($2, start, len) ) >split_recs
      written_split_recs++
   }
   typ[rname]==2 {
      start = 1; end = Pro_end[rname]; len = end - start + 1
      if(len>=min_recsize) {
         printf(">%s_bg-pro Pro %d len %d\n%s\n" ,$1, end,len, substr($2, start, len) ) >split_recs
         written_split_recs++
     }
     start = Phe_start[rname]; end = length($2); len = end - start + 1
      if(len>=min_recsize) {
         printf(">%s_phe-end Phe %d len %d\n%s\n" ,$1, end,len, substr($2, start, len) ) >split_recs
         written_split_recs++
     }
     start = Pro_start[rname]; end = Phe_end[rname]; len = end - start + 1
     printf(">%s_pro_CR_phe pro %d phe %d len %d\n%s\n", $1, start,end,len, substr($2, start, len) ) >cr_recs
     written_cr_recs++
   }
   END {printf("%d written to recs_linear_split_nocr.fa\n%d written to recs_cr_w_trna_flanks.fa\n", written_split_recs, written_cr_recs) > "split_msg"}
' $recinfo_tsv <(bawk '{print $name,$seq}' $fasta)

cd ..

log_split_msg
