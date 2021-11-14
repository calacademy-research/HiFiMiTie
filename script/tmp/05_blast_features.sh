#!/bin/bash

# we use the -task blastn and -evalue 1-e10 to get the smaller sequences but not too many false matches

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things

set_wdir # set the working dir variable wdir
[ -z $wdir ] && msglog "[$this]: The hfmt_<num> working directory not found" && exit 2

db_dir="${wdir}/rec_db"
db=${db_dir}/mito_rec_candidates
db_test_file=$db.nin
db_fasta=${wdir}/mito_rec_candidates.fasta

[ ! -d $db_dir ] && msglog "[$this]: $db_dir directory not found" && exit 2
[ ! -f $db_test_file ] && msglog "[$this]: $db does not appear to be intact"

dist_file=${wdir}/trnF_trnP_distances.tsv
stat_file=${wdir}/ControlRegion_length.stats

####### functions ######################################################
function blast_features_to_cand_recs {
   evalue=1e-10
   max_seqs=$(numrecs $db_fasta)
   query_file=${wdir}/mito_feature_sequences.fasta
   output_file=${wdir}/mito_feature_match_to_cand_recs.tsv

   start=$(date +%s)
   msglog $(date "+%T %d%b%Y")

   # blast features against candidate rec db and output in position sorted order
   msglog blastn -db $db -query $query_file -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -num_threads 8 -task blastn -evalue $evalue -max_target_seqs $max_seqs "|" sort -k2,2V -k9,9n -k12,12nr
   blastn -db $db -query $query_file -outfmt "6 std staxid stitle qlen qcovhsp qcovus" -num_threads 8 -task blastn -evalue $evalue -max_target_seqs $max_seqs \
      | sort -k2,2V -k9,9n -k12,12nr  \
      | bioawk_cas -t '{und=index($1,"_");feat=substr($1,1,und-1)}feat_lst==feat{next}{feat_lst=feat; print}' \
    > $output_file

   echo
   msglog $(date "+%T %d%b%Y")
   msglog "blastn completed in $(pprt_till_now $start) results in $output_file"
}

function create_dist_file {
   pos_tsv=${wdir}/mito_feature_match_to_cand_recs.tsv
   [ ! -f $pos_tsv ] && msglog "[$this]: $pos_tsv could not be found" && exit 2

   grep -e trnP -e trnF -e "^Phe" -e "^Pro" $pos_tsv | sed  "s|m[0-9]*_[0-9]*_[0-9]*/||" | \
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
}

function pro_phe_recs {
   grep -e "^Pro_" -e "trnP.*trnF" $dist_file
}

function get_stats_from_dist_file {
   mean=$(awk '{sum+=$NF; recs++}END{mean=sum/recs; print mean}' <(pro_phe_recs))
   stddev=$(awk -v mean=$mean 'function abs(a){return (a<0)? -a : a}{diffs+=abs($NF-mean); recs++}END{stddev=diffs/recs; print stddev}' <(pro_phe_recs))

   mean=$(printf "%.0f\n" $mean) # round to nearest integer
   stddev=$(printf "%.0f\n" $stddev) # round to nearest integer

   awk '{print $NF}' <(pro_phe_recs) |sort |uniq -c | sort -k1,1nr -k2,2n \
    | awk -v mean=$mean -v stddev=$stddev '
      BEGIN{print "    Num CR_len"; minlen = 1000000; mean=int(mean); stddev=int(stddev); sd_min=mean-stddev; sd_max=mean+stddev}

      $2>maxlen{maxlen=$2}
      $2<minlen{minlen=$2}

      # count how many are with in 1 stddev of the mean
      $2 >=mean && $2<=sd_max{within_1_sd += $1}
      $2 < mean && $2>=sd_max{within_1_sd += $1}

      {numrecs += $1}

      $1 != lstcount{
         if(lstcount>0) printf("\n")
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
         printf("\n")
         printf("\n    recs: %d   mean: %d stddev: %d   longest: %d shortest: %d diff: %d   num within 1 SD: %d\n", numrecs, mean, stddev, maxlen, minlen, maxlen-minlen, within_1_sd)
      }
   ' > $stat_file
}
###########################################################

# do the work
blast_features_to_cand_recs
create_dist_file
get_stats_from_dist_file
