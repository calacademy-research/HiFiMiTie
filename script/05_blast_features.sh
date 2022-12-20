#!/bin/bash

# we use the -task blastn and -evalue 1-e10 to get the smaller sequences but not too many false matches

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

db_dir=${wdir}/mito_hifi_rec_db
db=${db_dir}/mito_hifi_recs
db_test_file=$db.nin
db_fasta=${wdir}/mito_hifi_recs.fasta

[ ! -d $db_dir ] && msglog_module "$db_dir directory not found" && exit 2
[ ! -f $db_test_file ] && msglog_module "$db does not appear to be intact"

# blast_dir cr_dir set in shared.sh

feat_input_file=${wdir}/top_match_feature_sequences.fasta
feat_output_file=${blast_dir}/mito_feature_match_to_cand_recs.tsv

OH_output_file=${blast_dir}/OH_blast_to_cand_recs.tsv

####### functions ######################################################

function create_output_dirs { # make sure we have the blast_results and cr_analysis directories to place the output

   [ ! -d $blast_dir ] && mkdir $blast_dir && msglog_module "created $(basename $blast_dir) to hold the blast output tsv files"
   [ ! -d $blast_dir ] && msglog_module "directory $blast_dir could not be created" && exit 4

   [ ! -d $cr_dir ] && mkdir $cr_dir && msglog_module "created $(basename $cr_dir) to hold control region(s) analysis files"
   [ ! -d $cr_dir ] && msglog_module "directory $cr_dir could not be created" && exit 5
}

function blast_features_to_cand_recs {

   evalue=1e-10
   max_seqs=$(numrecs $db_fasta)
   query_file=$feat_input_file
   output_file=$feat_output_file

   start=$(date +%s)

   # blast features against candidate rec db and output in position sorted order
   msglog_module blastn -db $db -query $query_file -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -num_threads 8 -task blastn -evalue $evalue -max_target_seqs $max_seqs "|" sort -k2,2V -k9,9n -k12,12nr
   blastn -db $db -query $query_file -outfmt "6 std staxid stitle qlen qcovhsp qcovus" -num_threads 8 -task blastn -evalue $evalue -max_target_seqs $max_seqs \
      | sort -k2,2V -k9,9n -k12,12nr  \
      | bioawk_cas -t '{und=index($1,"_");feat=substr($1,1,und-1)}feat_lst==feat{next}{feat_lst=feat; print}' \
    > $output_file

   msglog_module "blastn completed in $(pprt_till_now $start) results in $output_file"
}

function blast_OH_to_cand_recs {
   qry=$db_fasta
   evalue=.00001
   output_file=$OH_output_file
   OH_db=${data_dir}/OrgRepl/OH.fas

   msglog_module blastn -db OH.fas -query $qry -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads
   blastn -db $OH_db -query $qry -outfmt  "6 std staxid stitle qlen qcovhsp qcovus"  -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads >$output_file
}

###########################################################

function trna_syms {
   AAsyms.sh F V L2 I Q M W A N C Y S2 D K G R H S1 L1 T P E
}

function alt_get_last_trna_count {
   first_trna_3let=$(AAsyms.sh $first_trna | awk '{print $2}')
   remove_non_trna_lines | sed "s/_/ /" | cut -f1 | sort | uniq -c | sort -k1,1nr | awk -v first_trna_3let=$first_trna_3let '
      $2 != first_trna_3let {
         print $2 ":" $1
         exit
   }'
}

function update_last_trna {
   last_trna=$(echo $last_trna_counts | cut -c1)

   if [ ! -z "$last_trna" ]; then
      update_setting_if_changed "last_trna" $last_trna
      update_setting_if_changed "mito_blast_last_trna_counts" "$last_trna_counts"
   fi
}

function get_last_trna {
   last_trna_counts=$( awk '
      FNR==NR {
         ar_let[$2] = $1   # e.g., ar_let["Phe"] = "F"
         next
      }

      {
         sub("_.*", "", $1)
         let1 = ar_let[$1]
         ar_feat_let[let1]++
      }
      END {
         PROCINFO["sorted_in"] = "@val_num_desc"
         for(t in ar_feat_let)
            printf("%s:%s, ", t, ar_feat_let[t])
         printf("\n")
      }
   ' <(trna_syms) <(get_last_trna_info) | sed "s/, *$//" )

   update_last_trna
}

function remove_non_trna_lines {  # remove lines if part before underscore is not one of the 3 letter trna ids
   awk '
      FNR==NR {
         ar_let[$2] = $1   # e.g., ar_let["Phe"] = "F"
         next
      }

      # see if part before underscore is in ar_let, e.g. in Phe_NC_007975 m64044_210611_022728/1639385/ccs_RC
      {
         und_ix = index($0, "_")
         if (und_ix < 1) next
         element = substr($0, 1, und_ix-1)
         if (element in ar_let)
            print
      }
   ' <(trna_syms) $feat_output_file
}

function get_last_trna_info { # create lines where last 2 fields are the first_trna match and the recid and the first 2 fields are the hit right before the first trna hit
   # Glu_NC_007975 m64044_210611_022728/1639385/ccs_RC       Phe_NC_007975 m64044_210611_022728/1639385/ccs_RC
   # Glu_NC_007975 m64044_210611_022728/15534806/ccs_RC      Phe_NC_007975 m64044_210611_022728/15534806/ccs_RC

   awk '
      FNR==NR {  # second field is the 3 letter trna symbol we look for (only a single line in "first file")
         mtch_first_trna= "^" $2
         next
      }

      $1 ~ mtch_first_trna {
         if (lst && lst_rec==$2) {
            print lst_feat,lst_rec "  \t"$1,$2
         }
      }
      {
         lst=$0;lst_rec=$2;lst_feat=$1
      }
   ' <(AAsyms.sh $first_trna) <(remove_non_trna_lines)
}

###########################################################

# do the work

create_output_dirs

run_if_no_file blast_features_to_cand_recs $feat_output_file
run_if_no_file blast_OH_to_cand_recs       $OH_output_file

# this uses the first trna in the feature list of the best match (usually F but can be others; eg in Arthropoda Lepidoptera is M)

source ${script_dir}/first_trna_from_features_or_taxname.sh

get_first_trna_setting
get_last_trna

# 14Dec2022 no longer use the alternate method here. we rely on the cm_results for the trnas, since it is more reliable

exit

if [ -z "$last_trna" ] || [ "$last_trna" == "$first_trna" ] || [ "$last_trna" == ":" ]; then
   msglog_module "Using alternative method to determine last_trna"

   last_trna_counts=$(alt_get_last_trna_count)
   update_last_trna
fi
