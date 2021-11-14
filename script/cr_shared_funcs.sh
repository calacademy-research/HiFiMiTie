#!/bin/bash

# presume source ${src_dir}/shared.sh has already been run in the script you include this in
# then you can do: source ${script_dir}/cr_shared_funcs.sh

[ ! -d $wdir_path ] && msg "the working directory path in wdir_path could not be found" && exit 2

# blast_dir cr_dir set in shared.sh

# example output in the dist file
# Thr_NC_007975	1639385/ccs_RC	91.429	70	4	2	1	70	4176	4243	Pro_NC_007975	1639385/ccs_RC	100.000	69	0	0	1	69	5674	5742	dist: 1500 trnas: 139 CR: 1361

####### functions ######################################################

function create_dist_file { # called with the 2 trnas that will flank the control region
   set_cr_filenames $1 $2 # sets fiveprime_flank threeprime_flank and dist_file vars

   pos_tsv=${blast_dir}/mito_feature_match_to_cand_recs.tsv
   [ ! -f $pos_tsv ] && msglog_module "$pos_tsv could not be found" && exit 2

   grepcmd=$(make_grep_command $fiveprime_flank $threeprime_flank)

   msglog_module $grepcmd $pos_tsv

   $grepcmd $pos_tsv | sed "s|m[0-9]*_[0-9]*_[0-9]*/||" | \
   bioawk_cas -t '
        function largerpos()  { return ($9 < $10) ? $10 : $9}
        function smallerpos() { return ($9 < $10) ?  $9 : $10}

        { recid = $2 }

        lstrd==recid {
            trnabeg = smallerpos()
            trnaend = largerpos()
            trnalen = trnaend - trnabeg + 1

            dist = trnaend - lstbeg + 1  # dist spans beginning of first trna to end of second trna in the CR sequence flanks
            trnalen = lst_trnalen + trnalen
            crlen = dist - trnalen

            print lst, fldcat(1,10), "dist: " dist " trnas: " trnalen " CR: " crlen
        }

        {
           lst = fldcat(1,10)
           lstrd  = recid
           lstbeg = smallerpos()
           lstend = largerpos()
           lst_trnalen = lstend - lstbeg + 1
        }
   ' > $dist_file

   [ -s $dist_file ] && msglog_module $(basename $dist_file) created
}

function get_stat_file_value {
   awk -v value=$1 '
      BEGIN{ valmatch="^" value }
      $1 ~ valmatch { print $2 }
   ' $stat_file
}

function get_stats_from_dist_file {
   trna1=$1; trna2=$2

   recs=$(awk '{sum+=$NF; recs++}END{mean=sum/recs; print recs}' <(trna1_trna2_recs $trna1 $trna2))
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
         if (maxlen == minlen) stddev = 0
         pct = (within_1_sd * 100) / numrecs
         printf("\n")
         printf("\nrecs: %d\nmean: %d\nstddev: %d\nmode: %d\nlongest: %d\nshortest: %d\ndiff: %d\nnum within one stddev: %d  %0.2f%%\n", numrecs, mean, stddev, mode, maxlen, minlen, maxlen-minlen, within_1_sd, pct)
      }
   ' > $stat_file

   if [ -s $stat_file ]; then # set some values for us to elsewhere if we want
      stddev=$(get_stat_file_value "stddev")
      shortest=$(get_stat_file_value "shortest")
      longest=$(get_stat_file_value "longest")
   fi

   [ -s $stat_file ] && msglog_module $(basename $stat_file) created
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

function trna1_trna2_recs { # was pro_phe_recs
   make_grep_command $1 $2 >/dev/null  # sets first_1_let, first_3_let and second_1_let that we use here

   grep -e "^${first_3_let}_" -e "${first_1_let}.*${second_1_let}" $dist_file
}

function set_cr_filenames {  # can call with two args for trna or use the first_trna and last_trna currently set
   fiveprime_flank=$1
   threeprime_flank=$2

   [ -z $fiveprime_flank ] && fiveprime_flank=$last_trna
   [ -z $threeprime_flank ] && threeprime_flank=$first_trna

   dist_file=${cr_dir}/trn${fiveprime_flank}_trn${threeprime_flank}_distances.tsv
   stat_file=${cr_dir}/ControlRegion_btw_trn${fiveprime_flank}_trn${threeprime_flank}_length.stats
}
####### end defs #######################################################
