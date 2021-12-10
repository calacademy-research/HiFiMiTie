#!/bin/bash

# test was oriented_mito_reads_ge_90_QCov.fa in stubifer/mito

fasta=$1
[ ! -f $fasta ] && fasta=mito_rec_candidates.fasta

threads=$2
[ -z "$threads" ] && threads=32

function msg { # write a msg to stderr
   >&2 echo -e $@
}

[ ! -f $fasta ] && msg "Could not find file: $fasta" && exit 1
prefix=$(echo $fasta | awk '{ sub(".gz$","",$1); sub(".fasta$","",$1); sub(".fa$","",$1); print $1; exit }')

# no round of 22, do not show results where eValue gt .001
addtl_ops="-onlycutoff"

function mitfi_rec_n {
   n=$1
   out=${n}.tmp_mitfi

   bawk -v n=$n 'NR==n{print ">"$name" "$comment;print $seq;exit}' $fasta  >${n}.tmpfa
   mitfi.sh $addtl_ops ${n}.tmpfa > $out

   rm ${n}.tmpfa
}

nrecs=$(numrecs $fasta)
msg "$nrecs records"

atatime=$threads
loops=$(awk -v nrecs=$nrecs -v atatime=$atatime 'BEGIN{print int( (nrecs+atatime-1)/atatime )}')

for (( i=0; i < $loops; ++i)); do
   msg "processing next $atatime"
   for j in $(seq $atatime); do

      ((rec=j+(i*atatime)))

      mitfi_rec_n $rec &

      if [ $rec -ge $nrecs ]; then
         break
      fi
   done
   wait
   msg "$rec processed of $nrecs"
done

cat *.tmp_mitfi >${prefix}.mitfi
rm  *.tmp_mitfi
