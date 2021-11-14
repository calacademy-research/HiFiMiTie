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

blast_dir=${wdir}/blast_results
cr_dir=${wdir}/cr_analysis

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

function create_dist_file {
   set_cr_filenames $1 $2 # sets fiveprime_flank threeprime_flank and dist_file vars

   pos_tsv=$feat_output_file
   [ ! -f $pos_tsv ] && msglog_module "$pos_tsv could not be found" && exit 2

   grepcmd=$(make_grep_command $fiveprime_flank $threeprime_flank)

   # grep -e trnP -e trnF -e "^Phe" -e "^Pro" $pos_tsv | sed  "s|m[0-9]*_[0-9]*_[0-9]*/||" | \
   msglog_module $grepcmd $pos_tsv

   $grepcmd $pos_tsv | sed  "s|m[0-9]*_[0-9]*_[0-9]*/||" | \
   bioawk_cas -t '
        function lrgpos() { return ($9<$10) ? $10:$9}
        lstrd==$2 {
            dist = lrgpos()-lstpos+1
            trnalen = lsttrna + $8
            crlen = dist - trnalen
            print lst, fldcat(1,10), "dist: " dist " trnas: " trnalen " CR: " crlen
        }

        {
           lst = fldcat(1,10)
           lstrd = $2
           lstpos = lrgpos()
           lsttrna = $8
        }
   ' > $dist_file

   [ -s $dist_file ] && msglog_module $(basename $dist_file) created
}

function trna1_trna2_recs { # was pro_phe_recs
   make_grep_command $1 $2 >/dev/null  # sets first_1_let, first_3_let and second_1_let that we use here

   # grep -e ^Pro_ -e "trnP.*trnF" $dist_file
   grep -e "^${first_3_let}_" -e "${first_1_let}.*${second_1_let}" $dist_file
}

function get_stats_from_dist_file {
   trna1=$1; trna2=$2

   mean=$(awk '{sum+=$NF; recs++}END{mean=sum/recs; print mean}' <(trna1_trna2_recs $trna1 $trna2))
   stddev=$(awk -v mean=$mean 'function abs(a){return (a<0)? -a : a}{diffs+=abs($NF-mean); recs++}END{stddev=diffs/recs; print stddev}' <(trna1_trna2_recs $trna1 $trna2))

   mean=$(printf "%.0f\n" $mean) # round to nearest integer
   stddev=$(printf "%.0f\n" $stddev) # round to nearest integer
   [ $stddev = 0 ] && stddev=1

   awk '{print $NF}' <(trna1_trna2_recs $trna1 $trna2) |sort |uniq -c | sort -k1,1nr -k2,2n \
    | awk -v mean=$mean -v stddev=$stddev '
      BEGIN{print "    Num CR_len"; minlen = 1000000; mean=int(mean); stddev=int(stddev); sd_min=mean-stddev; sd_max=mean+stddev}

      NR==1{ mode = $2 }

      $2>maxlen{maxlen=$2}
      $2<minlen{minlen=$2}

      # count how many are with in 1 stddev of the mean
      $2 >=mean && $2<=sd_max{within_1_sd += $1}
      $2 < mean && $2>=sd_min{within_1_sd += $1}

      {numrecs += $1}

      $1 != lstcount{
         if (lstcount>0)
            printf("\n")

         lstcount=$1
         printf("%s", $0)
         next
      }
      {
         printf(",%s", $2)
      }
      $2>=mean && $2<=sd_max{  # count how many are with in 1 stddev of the mean
        within_1_sd++
      }
      $2<mean && $2>=sd_max{  # count how many are with in 1 stddev of the mean
        within_1_sd++
      }
      END{
         pct = (within_1_sd * 100) / numrecs
         printf("\n")
         printf("\nrecs: %d\nmean: %d\nstddev: %d\nmode: %d\nlongest: %d\nshortest: %d\ndiff: %d\nnum within one stddev: %d  %0.2f%%\n", numrecs, mean, stddev, mode, maxlen, minlen, maxlen-minlen, within_1_sd, pct)
      }
   ' > $stat_file

   [ -s $stat_file ] && msglog_module $(basename $stat_file) created
}

function get_first_trna_setting {
   mm_fasta=$wdir_path/top_match_feature_sequences.fasta
   top_mm_first_trna=$(bawk '{print $name; exit}' $mm_fasta)
   [ ! -z $top_mm_first_trna ] && update_setting_if_changed "top_mito_first_trna" "$top_mm_first_trna"

   first_trna=$($script_dir/first_trna_from_features.sh $mm_fasta)
   [ ! -z $first_trna ] && update_setting_if_changed "first_trna" $first_trna
}

function get_last_trna {
   last_trna_counts=$( awk '
      FNR==NR{ ar_let[$2] = $1; next }  # e.g., ar_let["Phe"] = "F"
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
   ' <(AAsyms.sh F V L2 I Q M W A N C Y S2 D K G R H S1 L1 T P E) <(get_last_trna_info) | sed "s/, *$//" )

   last_trna=$(echo $last_trna_counts | cut -c1)
   update_setting_if_changed "mito_blast_last_trna_counts" "$last_trna_counts"
   update_setting_if_changed "last_trna" $last_trna
}
function get_last_trna_info { # create lines where last 2 fields are the first_trna match and the recid amd the first 2 fields are the hit roght before the first trna hit
   # Glu_NC_007975 m64044_210611_022728/1639385/ccs_RC       Phe_NC_007975 m64044_210611_022728/1639385/ccs_RC
   # Glu_NC_007975 m64044_210611_022728/15534806/ccs_RC      Phe_NC_007975 m64044_210611_022728/15534806/ccs_RC

   awk '
      FNR==NR{mtch_first_trna= "^" $2; next}

      $1 ~ mtch_first_trna {
         if (lst && lst_rec==$2)
            print lst_feat,lst_rec"  \t"$1,$2
      }
      {
         lst=$0;lst_rec=$2;lst_feat=$1
      }
   ' <(AAsyms.sh $first_trna) $feat_output_file
}

function make_grep_command {
   trna_1=$1; [ -z $trna_1 ] && trna_1=$last_trna
   trna_2=$2; [ -z $trna_2 ] && trna_2=$first_trna

   # one_letter_AA and three_letter_AA defined in shared.sh
   first_1_let="trn$(one_letter_AA $trna_1)"
   second_1_let="trn$(one_letter_AA $trna_2)"

   first_3_let=$(three_letter_AA $trna_1)
   second_3_let=$(three_letter_AA $trna_2)

   echo grep -e $first_1_let -e $second_1_let -e ^${first_3_let} -e ^${second_3_let}
}

function set_cr_filenames {  # can call with two args for trna or use the first_trna and last_trna currently set
   fiveprime_flank=$1
   threeprime_flank=$2

   [ -z $fiveprime_flank ] && fiveprime_flank=$last_trna
   [ -z $threeprime_flank ] && threeprime_flank=$first_trna

   dist_file=${cr_dir}/trn${fiveprime_flank}_trn${threeprime_flank}_distances.tsv
   stat_file=${cr_dir}/ControlRegion_btw_trn${fiveprime_flank}_trn${threeprime_flank}_length.stats
}


function create_wrap_around_cr_dist_file {
   create_dist_file $last_trna $first_trna
}
function get_wrap_around_cr_stats {
   get_stats_from_dist_file $last_trna $first_trna
}
###########################################################

# do the work

create_output_dirs

run_if_no_file blast_features_to_cand_recs $feat_output_file
run_if_no_file blast_OH_to_cand_recs       $OH_output_file

# this uses the first trna in the feature list of the best match (usually F but can be others; eg in Arthropoda Lepidoptera is M)
get_first_trna_setting
get_last_trna

set_cr_filenames $last_trna $first_trna
run_if_no_file create_wrap_around_cr_dist_file  $dist_file
run_if_no_file get_wrap_around_cr_stats         $stat_file
