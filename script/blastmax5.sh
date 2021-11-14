#!/bin/bash

# script to limit search to best 5 subjects hit for each query record

# arg1 or arg2 is the fasta query file
# arg2 or arg1 then is the blast db name
# optional arg 3 for number of threads, default 8
# optional arg 3 arg 4 for any more blast options

# output format has qlen qcovhsp qcovus along with staxid stitle
# qcovus looks all lines for same subject with same query and gives tot pct coverage of query

function msg { # write a msg to stderr
   >&2 echo -e $@
}

function usage {
   msg "\n   usage: blastmax5.sh <fasta query file> <blast db> [optional thread count] [optional additonal blastn options]\n"
   exit 1
}

[ -z $2 ] && usage

[ -f $1 ] && qry=$1 && db=$2
[ -z $qry ] && [ -f $2 ] && qry=$2 && db=$1
[ -z $qry ] && usage; [ -z $db ] && usage

is_fasta=1
is_gzipped=$(file -L $qry | grep "gzip compressed") # see if it is gzipped, if so we use fq2fa which handles gzipped or plaintext fasta or fastq files

# if it is a plaintext or gzipped fasta or fastq fastq, then we will use fq2fa
if [ -z "$is_gzipped" ]; then
   first_line=$(head -n1 $qry) && first_char=${first_line:0:1}
   if [ $first_char != ">" ]; then
      if [ $first_char == "@" ]; then
         unset is_fasta
      else
         msg "\n$qry should be a fasta or fastq file. Its first line is:"
         echo $first_line
         exit 2
      fi
   fi
else # looks like a gzipped file
   # msg $is_gzipped $qry
   unset is_fasta
fi

shift; shift

threads=8
re='^[0-9]+$'
[ ! -z $1 ] && [[ $1 =~ $re ]] && threads=$1 && shift

format="6 std staxid stitle qlen qcovhsp qcovus"

msg blastn -db $db -query $qry -outfmt "$format" -max_target_seqs 5 -num_threads $threads $@
if [ ! -z $is_fasta ]; then
   blastn -db $db -query $qry -outfmt "$format" -max_target_seqs 5 -num_threads $threads $@
else
   blastn -db $db -query <(fq2fa $qry) -outfmt "$format" -max_target_seqs 5 -num_threads $threads $@
fi
