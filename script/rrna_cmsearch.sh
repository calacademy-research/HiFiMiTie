#!/bin/bash

# version 2 of this use rrnS.cm and rrnL.cm instead of 12S_s-rna.cm and 16S_l-rna.cm in the models dir
# also use these mitos2 style args on the cmsearch command line: --noF4b --cpu 1 --notextw -E 0.01 --mxsize 80000
# and always do local search not global

# usually called with mito_megahit.fasta or mito_kalign.fasta

: ' # creates 2 output files like:
#target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
mito                 -         12snew_curated.47-1  -          cm        1      952       72     1025      +    no    1 0.38  60.7  845.2   8.5e-59 !   megahit	79_3
...
'
src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir_path ] && wdir="." && wdir_path="."  #msglog_module "The hfmt_<num> working directory not found" && exit 2

mito_fasta=$1
[ -z $mito_fasta ] && msg "    usage: rrna_cmsearch.sh <mito fasta file> [<tbl_prefix>]" && exit 1

tbl_prefix=$2

function check_files {
    [ ! -f $mito_fasta ] && msglog_module "could not find $mito_fasta" && exit 3

    [ -z $(which cmsearch) ] && msglog_module "cmsearch program could not be found" && exit 3

    cm_12S="${src_dir}/models/rrnS.cm"
    cm_16S="${src_dir}/models/rrnL.cm"

    [ ! -f ${cm_12S} ] && msglog_module "could not model ${12S_cm}" && exit 4
    [ ! -f ${cm_16S} ] && msglog_module "could not model ${16S_cm}" && exit 4

    tbl_12S=$(basename $cm_12S .cm).tbl
    tbl_16S=$(basename $cm_16S .cm).tbl

    if [ ! -z $tbl_prefix ]; then
       tbl_12S=${tbl_prefix}_$tbl_12S
       tbl_16S=${tbl_prefix}_$tbl_16S
    fi
}

# simply reports whether line 3 of file in arg1 is a comment or not
function have_tbl_hit {
   tbl=$1
   [ ! -f "$tbl" ] && false && return

   ln3_ch=$(awk 'BEGIN{ch="#"}NR==3{ch=substr($0,1,1)}END{print ch}' $tbl)
   [ ln3_ch == "#" ] && false && return

   true
}

function run_cmsearch {
   tbl=$1
   cmout=$(basename $tbl .tbl).cmout

   # 08Aug2022 could be lots of sequences use what is available
   local cpu_num=$threads
   [ -z $cpu_num ] && cpu_num=1

   cm=$2
   # do global search first
   # cmsearch -g --tblout $tbl $cm $mito_fasta >/dev/null

   # if that does not work do a local search, which we accomplish by just not adding in the -g option
   # ! have_tbl_hit $tbl && cmsearch --tblout $tbl $cm $mito_fasta >/dev/null

   # these rrna settings were determined by looking at the mitos2 cmout output files
   cmsearch --noF4b --cpu $cpu_num --notextw -E 0.01 --mxsize 80000 --tblout $tbl $cm $mito_fasta > $cmout
}

cd $wdir_path
check_files
msglog_module "input files checked"

msglog_module "search for 12S_s-rna"
run_cmsearch $tbl_12S $cm_12S

msglog_module "search for 16S_l-rna"
run_cmsearch $tbl_16S $cm_16S
