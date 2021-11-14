#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

function basic_note {
   local num_reads=$(get_setting HiFi_mito_reads)

   echo -e "
mitochondrion.fasta has the selected mito sequence derived from the $num_reads HiFi reads.
mitochondrion.anno has the basic annotation of tRNAs, genes, and the two rRNAs.

The name of the sequence will likely need to be editted to something more appropriate.

The annotations can be checked by runnng the fasta sequence on the Mitos2 web server.
Mitos2 improves protein coding gene annotations and uses more than just blastn to the
top matching sequence from GenBank, which is what is represented here.

The areas of difference between the 2 assemblies can be seen in the log file and also
in the ${wdir}/compare_megahit_msa/comparison_overview.txt file. In the same directory
is a fullsequence_comparison.txt that has some annotations to orient locations.

The ${wdir}/cm_results/one_line_per_rec.cm_anno.srt file gives an overview of the HiFi read
sequences individual annotations and can show some types of non control region heteroplasmy if it exists.
"
}
