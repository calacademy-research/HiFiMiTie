#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

######################################################################
#        set the path and file name variables we will need           #
######################################################################

mito_fasta=${wdir_path}/mito_rec_candidates.fasta
[ ! -s "$mito_fasta" ] && msglog_module "$mito_fasta missing or empty. could not run cm_search on it." && exit 3

cmdir=cm_results
cmdir_path=$wdir_path/$cmdir
[ ! -d $cmdir_path ] && mkdir $cmdir_path && msglog_module "$cmdir directory created"

pfx=rrna
tbl_12S=${pfx}_rrnS.tbl
tbl_16S=${pfx}_rrnL.tbl

######################################################################
#                       function definitions                         #
######################################################################

# look for 12S and 16S rrna sequences
function search_12S_16S {
   cd $wdir_path # rrna_cmsearch.sh will write to file in the current directory, so keep it here

   $(get_script rrna_cmsearch) $mito_fasta $pfx
   mv ${wdir_path}/${pfx}_* $cmdir_path/

   msglog_module "$tbl_12S $tbl_16S created in $cmdir"
}

# look for all the trna sequences using cm_search
function tRNA_searches {
   cd $cmdir_path
   msglog_module "tRNA search mito reads using: MiTFi - mitochondrial tRNA finder"

   ${script_dir}/multi_mitfi.sh $mito_fasta $threads # creates mito_rec_candidates.mitfi

   [ -s ${wdir_path}/mito_rec_candidates.mitfi ] && mv ${wdir_path}/mito_rec_candidates.mitfi $cmdir_path/
   [ -s ${cmdir_path}/mito_rec_candidates.mitfi ] && msglog_module "mito_rec_candidates.mitfi created in $cmdir"
   [ ! -s ${cmdir_path}/mito_rec_candidates.mitfi ] && msglog_module "Problem creating mito_rec_candidates.mitfi in $cmdir"
}

# add the 12S, 16S and any goose hairpins found into the mitfi contents and output to mito_rec_candidates.ccm_anno file
function create_cm_anno {
   cd $cmdir_path

   ${script_dir}/add_goosehairpin_to_mitfi.sh          | \
    ${script_dir}/add_rrna_to_mitfi.sh rrna_rrnS.tbl - | \
    ${script_dir}/add_rrna_to_mitfi.sh rrna_rrnL.tbl - | \
    awk '/^#header/{if(NR==1)next; $0="#"} {print}' >mito_rec_candidates.cm_anno

   msglog_module "mito_rec_candidates.cm_anno created with mitfi, goose_hairpin, 12S_rna and 16S_rna cm results"
}

# create the right neighbor matrix, ordered by most frequent neighbor seen. this will define our tRNA ordering for additonal analysis (CR etc)
function create_trna_right_neighbor_matrix {
   cd $cmdir_path

   starting_trna=$(get_setting_or_default "first_trna" "F")
   ${script_dir}/right_neighbor_matrix.sh mito_rec_candidates.mitfi $starting_trna >trna_right_neighbor.matrix

   order=$(head -1 trna_right_neighbor.matrix | sed "s/\t/  /g")
   msglog_module "trna_right_neighbor.matrix created using ${BBlue}${starting_trna}${NC} as the first tRNA."
   msglog_module "trna order:" $order

   # remember the last trna in the settings file, for use in defining control region location
   last_trna=$(tail -1 trna_right_neighbor.matrix | cut -c1)
   [ ! -z "$last_trna" ] && update_setting "last_trna" $last_trna

   update_setting "trna_order" "$order"
}

# create the same matrix but including 12S, 16S and potentially any found goose_hairpins, store it in cm_anno_right_neighbor.matrix
function create_cm_anno_right_neighbor_matrix {
   cd $cmdir_path

   starting_trna=$(get_setting_or_default "first_trna" "F")
   ${script_dir}/right_neighbor_matrix.sh mito_rec_candidates.cm_anno $starting_trna >cm_anno_right_neighbor.matrix

   msglog_module "cm_anno_right_neighbor.matrix created using ${BBlue}${starting_trna}${NC} as the first tRNA."

   # create an empty file that shows the cm_anno order
   > $(head -1 cm_anno_right_neighbor.matrix | sed "s/\t/__/g" | sed "s/ //g")__

   order=$(head -1 cm_anno_right_neighbor.matrix | sed "s/\t/  /g")
   update_setting "cm_anno_order" "$order"
   msglog_module "cm_anno rna order:" $order
}

function create_matrix_and_anno_files {
   run_if_no_file create_trna_right_neighbor_matrix    ${cmdir_path}/trna_right_neighbor.matrix
   run_if_no_file create_cm_anno                       ${cmdir_path}/mito_rec_candidates.cm_anno
   run_if_no_file create_cm_anno_right_neighbor_matrix ${cmdir_path}/cm_anno_right_neighbor.matrix
}

##########################################################################
#                                                                        #
#   call the functions if the files they create are not in $cmdir_path   #
#                                                                        #
##########################################################################

run_if_no_file search_12S_16S ${cmdir_path}/$tbl_16S
run_if_no_file tRNA_searches ${cmdir_path}/mito_rec_candidates.mitfi

if [ -s ${cmdir_path}/mito_rec_candidates.mitfi ]; then
   create_matrix_and_anno_files
else
   msglog_module "error: ${cmdir_path}/mito_rec_candidates.mitfi was not found"
fi
