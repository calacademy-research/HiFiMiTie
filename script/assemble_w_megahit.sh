#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

# could use the hfmt thread setting but megahit does fine with 1 on such small input
threads=1
# mito_hifi_recs.fasta
function run_megahit {
    local mito_fasta_recs_to_assemble=mito_hifi_recs.fasta

    msglog_module "running megahit to assemble mito records using sequences in $mito_fasta_recs_to_assemble"
    msglog_module "megahit -t $threads -r $mito_fasta_recs_to_assemble"
    megahit -t $threads -r $mito_fasta_recs_to_assemble

    msglog_module "cd megahit_out"
    cd megahit_out

    write_megahit_best
    msglog_module "megahit_best.fa mito file created"
}

# usually the sequence with the highest coverage is the full mito sequence.
# however when there are a lot of duplicated subsequences such as in CR tandem repeats
# this might not be the case. So look for hight coverage and seq length >10,000 bases
# write out the non-mito (ie <10000bp) high coverage sequences into megahit_hicov.fa
# and write the mito seq into megahit_best.fa
function write_megahit_best {
   rm -f megahit_best.fa megahit_hicov.fa
   cawk 'BEGIN{mbest = "megahit_best.fa"; mhicov = "megahit_hicov.fa"; mitomin = 10000}
      length($2) < mitomin { printf(">%s %s\n%s", $1, fldcat(4,NF), $2) >>mhicov; next }
      { printf(">%s %s\n%s", $1, fldcat(4,NF), $2) >mbest ; exit }
   ' <( bawk '{print}' final.contigs.fa | sort -k4,4Vr)
}

function remove_low_evalues {
   awk '
      /^#/{print;next}
      $5 < 0.001
   '
}

function run_mitfi_on_best {
   # expecting to be in the megahit_out dir at this point
   # this is to determine where the first tran,  Phe or Pro usually, is located so we can reorient the sequence to start at the correct tRNA

   msglog_module "running MiTFi to determine location of $first_trna_3let so that the sequence will start there"
   mitfi.sh megahit_best.fa | remove_low_evalues >megahit_best.mitfi
   msglog_module "MiTFi run on megahit_best.fa completed"
}

# k79_3   15352   15420   62.71   1.718E-12       GAA     F       Metazoa_F.cm    -

# reorient to have it start with the first_trna, usually Phe
function reorient_seq {
    local mitfi_fail_msg1="Could not find megahit_best.mitfi annotation file"
    local mitfi_fail_msg2="Could not find $first_trna_3let entry to use for reorienting the mito contig to start at $first_trna_3let"

    [ ! -f megahit_best.mitfi ] && msglog_module $mitfi_fail_msg1 && exit 3

    tPhe_info=$(awk -v first_trna="$first_trna" '
       /^#/{next}
       $7==first_trna && start==0{
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

    [ "$tPhe_start" -eq 0 ] && msglog_module $mitfi_fail_msg2 && exit 4

    # use bioawk to make it easier to revcomp if need be
    bawk -v tPhe_start="$tPhe_start" -v tPhe_dir="$tPhe_dir" '
       NR==1 {
          print ">megahit_mitochondrion " $name " " $comment
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

# search for 12S and 16S rrna sequences in the fasta file of arg 1, directory to store results in arg2
function search_12S_16S {

   local mito_fasta=$(realpath $1)
   local mh_dir=$(realpath $2)
   local pfx=$(get_fastx_basename $mito_fasta)

   tbl_12S=${pfx}_rrnS.tbl
   tbl_16S=${pfx}_rrnL.tbl

   run_if_no_file noop $mh_dir/$tbl_16S
   [ -s $mh_dir/$tbl_16S ] && return 0  # already made the file

   $(get_script rrna_cmsearch) $mito_fasta $pfx

   [ -s ${wdir_path}/$tbl_12S ] && mv ${wdir_path}/${pfx}_rrnS.* $mh_dir/
   [ -s ${wdir_path}/$tbl_16S ] && mv ${wdir_path}/${pfx}_rrnL.* $mh_dir/

   [ -s $mh_dir/$tbl_12S ] && [ -s $mh_dir/$tbl_16S ] && msglog_module "$tbl_12S $tbl_16S created in $(basename $mh_dir) for $(basename $mito_fasta)"
}

function MiTFi_megahit {
   msglog_module "Running MiTFi on the reoriented sequence."
   mitfi.sh mito_megahit.fasta | remove_low_evalues >$mitfi_path
}

function run_mitfi_on_reoriented_sequence {
   local mitfi=$mitfi_path
   cm_anno=$(basename ${mitfi} .mitfi).cm_anno
   cm_anno_path=$(dirname $mitfi)/$cm_anno
   anno_path=$(dirname $mitfi)/$(basename ${mitfi} .mitfi).anno

   run_if_no_file MiTFi_megahit $mitfi
   search_12S_16S $wdir_path/megahit_out/mito_megahit.fasta $wdir_path/megahit_out

   [ -s $mitfi ] && add_goosehairpin_to_mitfi.sh $mitfi mito_megahit.fasta >${mitfi}_w_gh
   [ -s ${mitfi}_w_gh ] && msglog_module "goose hairpin sequence information added, if any found"

   if [ -s megahit_out/$tbl_16S ]; then
      asm_length=$(asmlen megahit_out/mito_megahit.fasta 2>/dev/null)
      settings=$(settings_file)

      ${script_dir}/add_rrna_to_mitfi.sh megahit_out/$tbl_12S ${mitfi}_w_gh | \
      ${script_dir}/add_rrna_to_mitfi.sh megahit_out/$tbl_16S -             | \
      ${script_dir}/add_CR_to_mitfi.sh $settings - $asm_length > $cm_anno_path

      msglog_module "$cm_anno with rrns, rrnL, cr and any goose hairpins added to mitfi results created for mito_megahit.fasta"

      if [ -s $cm_anno_path ]; then  # 01Dec2022 add genes to anno
         [ -s ${mitfi}_w_gh ] && rm ${mitfi}_w_gh
         ${script_dir}/mito_pcg_anno.sh mito_megahit.fasta $cm_anno_path | add_header_to_anno >$anno_path
         [ -s $anno_path ] && msglog_module "$(basename $anno_path) with genes added to $cm_anno created for mito_megahit.fasta"
      fi

   fi
}

function make_softlinks {
   # add softlink to megahit results into the hfmt directory
   [ ! -s mito_megahit.fasta ] && [ -s megahit_out/mito_megahit.fasta ] && ln -s megahit_out/mito_megahit.fasta

   # make softlink to the mitfi and other annotations in the hfmt directory
   [ ! -f $cm_anno ] && [ -s $cm_anno_path ] && ln -s $cm_anno_path $cm_anno
   [ ! -f $anno ] && [ -s $anno_path ] && ln -s $anno_path $anno

   touch -h mito_megahit.fasta
   [ -f $cm_anno ] && touch -h $cm_anno
   [ -f $anno ] && touch -h $anno
}

##############################################################################
#                             do the work                                    #
##############################################################################

cd $wdir_path

# get tRNA that we want to start the mito sequence. one letter and 3 letter versions
first_trna=$(get_setting_or_default "first_trna" "F")
first_trna_3let=$(three_letter_AA $first_trna)

# run megahit on our sequences. this changes us to the megahit_out dir after megahit is run
run_if_no_file run_megahit $wdir_path/megahit_out/megahit_best.fa
cd $wdir_path/megahit_out  # explicitly cd since run_megahit might not be called

# run mitfi so we can use its output to find first_trna
run_if_no_file run_mitfi_on_best $wdir_path/megahit_out/megahit_best.mitfi

# reorient the sequence based on the mitfi output (todo: small chance first_trna is split in megahit_best.fa, handle this)
reorient_seq >$wdir_path/megahit_out/mito_megahit.fasta
msglog_module "mito_megahit.fasta reoriented to start with ${first_trna_3let}."

# pop back to wdir for rest of the script
cd $wdir_path

# make sure to define these before run_mitfi_on_reoriented_sequence called
mitfi_path=megahit_out/mito_megahit.mitfi
cm_anno=$(basename $mitfi_path .mitfi).cm_anno
cm_anno_path=$(dirname $mitfi_path)/$cm_anno
anno=$(basename $mitfi_path .mitfi).anno
anno_path=$(dirname $mitfi_path)/$anno

make_softlinks

# run mitfi on the reoriented fasta to use for QA and comparison purposes
[ -s mito_megahit.fasta ] && run_if_no_file run_mitfi_on_reoriented_sequence megahit_out/mito_megahit.cm_anno

make_softlinks # this call will make the additional cm_anno link and update the other time on the fasta link

msglog_module "MiTFi annotation of mito_megahit.fasta completed, results in $(basename $anno_path)"
