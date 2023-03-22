#!/bin/bash

# pass in a file where last field in each line is a CR length
function get_stats_from_CR_lens {
   local CR_lens=$1

   recs=$(awk '{sum+=$NF; recs++}END{mean=sum/recs; print recs}' $CR_lens)
   mean=$(awk '{sum+=$NF; recs++}END{mean=sum/recs; print int(mean + 0.5)}' $CR_lens )
   stddev=$(awk -v mean=$mean 'function abs(a){return (a<0)? -a : a}{diffs+=abs($NF-mean); recs++}END{stddev=diffs/recs; print int(stddev + 0.5)}' $CR_lens )
   [ $stddev = 0 ] && stddev=1

   awk '{print $NF}' $CR_lens |sort |uniq -c | sort -k1,1nr -k2,2n |
   awk -v mean=$mean -v stddev=$stddev '
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
   '
}

function run_get_stats_with_tmp {
   # file is reused so in case we read from process result or stin, make a tmp file
   tmp=cr_lens_${RANDOM}_${RANDOM}
   cp $1 $tmp

   get_stats_from_CR_lens $tmp

   rm $tmp
}

# make this so if we do not have $1 arg we can source it

[ ! -z "$1" ] && run_get_stats_with_tmp $1
