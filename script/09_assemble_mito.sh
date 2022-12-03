#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

function run_assemble_w_megahit {
   msglog_module "Step 9a -- assemble using megahit"
   $(get_script assemble_w_megahit)
}

function run_alignment_assembly {
   msglog ""
   msglog_module "Step 9b -- assemble using multi-sequence alignment (msa) consensus"

   $(get_script assemble_w_mafft_and_consensus.sh)
}

function create_msa_asm_softlink {
   local linkdir=$(basename $alignasm_dir)

   [  ! -s ${wdir}/mito_msa.fasta  ] \
   && [ -s ${alignasm_dir}/mito_msa.fasta ] \
   && ln -s ${linkdir}/mito_msa.fasta ${wdir}/mito_msa.fasta

   [  ! -s ${wdir}/mito_msa.cm_anno  ] \
   && [ -s ${alignasm_dir}/cm_mitfi/mito_msa.cm_anno ] \
   && ln -s ${linkdir}/cm_mitfi/mito_msa.cm_anno ${wdir}/mito_msa.cm_anno

    # 02Dec2022 add anno with genes as well as trna and rrna
   [  ! -s ${wdir}/mito_msa.anno  ] \
   && [ -s ${alignasm_dir}/cm_mitfi/mito_msa.anno ] \
   && ln -s ${linkdir}/cm_mitfi/mito_msa.anno ${wdir}/mito_msa.anno

}

# get the kmer based assembly -- ${megahit_dir}/megahit_best.fa
run_if_no_file run_assemble_w_megahit ${wdir}/mito_megahit.fasta

# get the sequence alignment, msa and consensus assembly
run_if_no_file run_alignment_assembly ${alignasm_dir}/mito_msa.fasta
create_msa_asm_softlink
