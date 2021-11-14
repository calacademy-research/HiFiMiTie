#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

# could use the hfmt thread setting but megahit does fine with 1 on such small input
threads=1

function run_megahit {
    mito_fasta_recs_to_assemble=recs_linear_split_nocr.fa

    msglog_module "running megahit to assemble mito records"
    megahit -t $threads -r $mito_fasta_recs_to_assemble

    cd megahit_out

    # this should be name of best matching contig in the megahit.out directory
    ctg=$(grep "^>" final.contigs.fa | sed "s/multi=//" | sort -k3,3Vr | awk '{print substr($1,2);exit}')

    # we will not fold it since we are going to need to reorient it to start with tPhe
    awk -v ctg="${ctg}" '
       $1 == ">"ctg {
          print $0
          prt=1
          next
       }
       prt && /^>/ {exit}
       prt { print }
    ' final.contigs.fa >megahit_best.fa
    msglog_module "megahit_best.fa mito file created"
}

function run_mitfi {
   # this is to determine where tPhe or tPro is located so we can reorient the sequence to start tPhe
   msglog_module "running MiTFi to determine location of tPhe so that the sequence will start there"
   mitfi.sh megahit_best.fa >megahit_best.mitfi
   msglog_module "MiTFi run completed"
}

# k79_3   15352   15420   62.71   1.718E-12       GAA     F       Metazoa_F.cm    -

# reorient to have it start with tPhe
function reorient_seq {
    mitfi_fail_msg="Could not find megahit_best.mitfi annotation file or the tPhe entry to use for reorienting the mito contig to start at tPhe"
    [ ! -f megahit_best.mitfi ] && msglog_module $mitfi_fail_msg  && exit 3

    tPhe_info=$(awk '
       /^#/{next}
       $7=="F" && start==0{
          if($9!="-") {
             dir = "+"
             start = $2
          } else {
             dir = "-"
             start = $3
         }
       }
       END{printf "%d %s\n", start, dir}
    ' megahit_best.mitfi)

    tPhe_start=$(echo "$tPhe_info" | awk '{print $1}')
    tPhe_dir=$(echo "$tPhe_info" | awk '{print $2}')

    [ "$tPhe_start" -eq 0 ] && msglog_module $mitfi_fail_msg && exit 4

    # use bioawk to make it easier to revcomp if need be
    bawk -v tPhe_start="$tPhe_start" -v tPhe_dir="$tPhe_dir" '
       NR==1 {
          print ">mito megahit " $name " " $comment
          if (tPhe_dir=="-") {
             tPhe_start = length($seq) - tPhe_start + 1
             revcomp($seq)
          }

          if (tPhe_start == 1)
             print $seq
          else
             print substr($seq, tPhe_start) substr($seq, 1, tPhe_start-1)
       }
    ' megahit_best.fa
}

##############################################################################
#                             do the work                                    #
##############################################################################

cd $wdir

run_megahit

run_mitfi
reorient_seq >mito_megahit.fasta

cp mito_megahit.fasta ../mito_megahit.fasta
cd ..

msglog_module "megahit mito fasta starting with tPhe in $wdir/mito_megahit.fasta.  Rerun MiTFi on the reoriented record."
mitfi.sh mito_megahit.fasta >mito_megahit.mitfi
msglog_module "MiTFi on mito_megahit.fasta run completed, results in $wdir/mito_megahit.mitfi"
