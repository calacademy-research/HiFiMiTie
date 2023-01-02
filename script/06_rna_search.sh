#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

######################################################################
#        set the path and file name variables we will need           #
######################################################################

mito_fasta=${wdir_path}/mito_hifi_recs.fasta
[ ! -s "$mito_fasta" ] && msglog_module "$mito_fasta missing or empty. could not run cm_search on it." && exit 3

cmdir=cm_results
cmdir_path=$wdir_path/$cmdir
[ ! -d $cmdir_path ] && mkdir $cmdir_path && msglog_module "$cmdir directory created"

pfx=rrna
tbl_12S=${pfx}_rrnS.tbl
tbl_16S=${pfx}_rrnL.tbl
cmout_16S=${pfx}_rrnL.cmout
OH_blast_path=$wdir_path/blast_results/OH_blast_to_cand_recs.tsv

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

function do_OL_search {
   cd $exec_path
   search_OL
   [ -s "${cmdir_path}/OL.tbl" ] && msglog_module "OL.tbl created in $cmdir"
}

# look for all the trna sequences using cm_search
function tRNA_searches {
   cd $cmdir_path
   msglog_module "tRNA search mito reads using: MiTFi - mitochondrial tRNA finder"

   ${script_dir}/multi_mitfi.sh $mito_fasta $threads # creates mito_hifi_recs.mitfi

   [ -s ${wdir_path}/mito_hifi_recs.mitfi ] && mv ${wdir_path}/mito_hifi_recs.mitfi $cmdir_path/
   [ -s ${cmdir_path}/mito_hifi_recs.mitfi ] && msglog_module "mito_hifi_recs.mitfi created in $cmdir"
   [ ! -s ${cmdir_path}/mito_hifi_recs.mitfi ] && msglog_module "Problem creating mito_hifi_recs.mitfi in $cmdir"
}

function add_dist_from_last { # also cleans up headers
   awk '
      BEGIN{FS="\t"; OFS="\t"}

      /^#header/ {
         if (NR==1) {print} else if (! lst_hdr) {print "#"}
         lst_hdr = 1
         next
      }
      { lst_hdr = 0 }

      lst==$1 { $0 = $0 "\t" $2-epos }
      lst!=$1 && NF==9 { $0 = $0 "\t." }
      { print }
      { lst = $1; epos = $3 }
   '
}

function add_counts_to_settings {
   one_liner_annos=$1
   [ ! -s $one_liner_annos ] && return

   tot_recs=$(tail -n +2 $one_liner_annos | wc -l)
   prob_annos=$(grep -c "#" $one_liner_annos)

   # count the reads that are candidates for having the whole non-cr sequence in them
   local starting_trna=$(get_setting_or_default "first_trna" "F")
   local last_trna=$(get_setting_or_default "last_trna" "P")
   qry="^m.*${starting_trna}.*${last_trna}"
   complete_non_cr_recs=$(grep $qry $one_liner_annos | grep -c -v "# T")  # the "# T" excludes those that do not match template

   update_setting_if_changed "cm_tot_anno_recs" "$tot_recs"
   update_setting_if_changed "cm_low_qual_annos" "$prob_annos"
   update_setting_if_changed "cm_complete_no_cr_anno" "$complete_non_cr_recs"
}

function create_cm_anno_one_liners {
   one_liner=$cmdir_path/one_line_per_rec.cm_anno
   ${script_dir}/one_liner_per_cm_anno_rec.sh > $one_liner

   # now make one sorted and also a file showing distribution counts of the starting trna or other elements per line, useful if subsampling warranted
   if [ -s $one_liner ]; then
      ${script_dir}/one_liner_per_cm_anno_rec.sh "sort" > ${one_liner}.srt

      if [ -s ${one_liner}.srt ]; then
         # create starting element distribution file, exclude suspect lines, i.e., those with a comment char
         dist_file=$(anno_start_item_distribution_file)
         grep -v "#" ${one_liner}.srt | awk 'NR>1{sub("~","",$2); print $2}' | uniq -c > $dist_file
         add_counts_to_settings ${one_liner}.srt
      fi
   fi

   [ -s ${one_liner}.srt ] && msglog_module "$(basename ${one_liner}.srt) created"
}

# add the 12S, 16S and any goose hairpins found into the mitfi contents and output to mito_hifi_recs.cm_anno file
function create_cm_anno {
   cd $cmdir_path

   ${script_dir}/add_goosehairpin_to_mitfi.sh          |
   ${script_dir}/add_rrna_to_mitfi.sh rrna_rrnS.tbl -  |
   ${script_dir}/add_rrna_to_mitfi.sh rrna_rrnL.tbl -  |
   ${script_dir}/add_rrna_to_mitfi.sh OL.tbl        -  |
   ${script_dir}/add_OH_to_mitfi.sh $OH_blast_path  -  |
   add_dist_from_last > mito_hifi_recs.cm_anno

   msglog_module "mito_hifi_recs.cm_anno created with mitfi, goose_hairpin, 12S_rna and 16S_rna cm results"
}

# create the right neighbor matrix, ordered by most frequent neighbor seen. this will define our tRNA ordering for additonal analysis (CR etc)
function create_trna_right_neighbor_matrix {
   cd $cmdir_path

   starting_trna=$(get_setting_or_default "first_trna" "F")
   ${script_dir}/right_neighbor_matrix.sh mito_hifi_recs.mitfi $starting_trna >trna_right_neighbor.matrix

   order=$(head -1 trna_right_neighbor.matrix | sed "s/\t/ /g" | sed "s/^ *//")
   msglog_module "trna_right_neighbor.matrix created using ${BBlue}${starting_trna}${NC} as the first tRNA."
   msglog_module "trna order:" $order

   # remember the last trna in the settings file, for use in defining control region location
   last_trna=$(tail -1 trna_right_neighbor.matrix | cut -c1)
   [ ! -z "$last_trna" ] && update_setting "last_trna" $last_trna

   update_setting "trna_order" "$order"
}

function goosehairpin_trnaset {
   awk '
      function trna_prev(p) {  # passes in i-1 and we return $p unless it is OH where we will return $(p-1)
         return ($p != "OH") ? $p : $(p-1)
      }

      NR==1 {
      for (i=1;i<=NF;i++) {
         if ($i=="gh") {
            if (i==NF)  # P F case
               print trna_prev(i-1), $1
            else if (i==1) # should not happen but just to cover edge case
               print trna_prev(NF), $2
            else # T P case
               print trna_prev(i-1), $(i+1)
            exit
         }
      }
   }' cm_anno_right_neighbor.matrix
}
function find_goosehairpin_bounding_trnas {
   trna_set=$(goosehairpin_trnaset)

   first_let=$(echo $trna_set | cut -c1)
   last_let=$(echo $trna_set  | cut -c3)
   [ ! -z $first_let ] && update_setting "gh_prev_trna" $first_let
   [ ! -z $last_let ]  && update_setting "gh_succ_trna" $last_let
}

# 06Sep2022 normally existing mito_hifi_recs.cm_anno will work to create_cm_anno_right_neighbor_matrix
# but if there are multiple gh in the there it will not work yet
# so for now remove the gh is the first time did not work
function get_mito_hifi_recs_cm_anno {
   if [ -z $filter_gh ] && [ -z $filter_OH ]; then
      cat mito_hifi_recs.cm_anno
   elif [ ! -z $filter_OH ]; then
      grep -v -e gh -e OH mito_hifi_recs.cm_anno
   else
      grep -v gh mito_hifi_recs.cm_anno
   fi
}
function cm_matrix_ok { # check to see if we had a problem by counting items in header. has to be at least what we had for the trna only version
   ok=$(awk 'NR==1{trna_items=NF; next}
             FNR!=NR && FNR==1{cm_anno_items=NF; exit}
             END{if(cm_anno_items >= trna_items) print "ok"; else print "no"}
        ' trna_right_neighbor.matrix cm_anno_right_neighbor.matrix)

   if [ $ok = "no" ]; then
      false
   else
      true
   fi
}

function create_cm_anno_right_neighbor_matrix {
   unset filter_gh; unset filter_OH

   run_create_cm_anno_right_neighbor_matrix

   if ! cm_matrix_ok ; then
      msglog_module "too few items found in cm_anno_right_neighbor_matrix. rerun removing goosehairpin lines"
      filter_gh="filter_gh"
      run_create_cm_anno_right_neighbor_matrix
      if ! cm_matrix_ok; then
         msglog_module "too few items found in cm_anno_right_neighbor_matrix. rerun removing goosehairpin and OH lines"
         filter_gh="filter_OH"
         run_create_cm_anno_right_neighbor_matrix
      fi
   fi
}

# create the same matrix but including 12S, 16S and potentially any found goose_hairpins, store it in cm_anno_right_neighbor.matrix
function run_create_cm_anno_right_neighbor_matrix {
   cd $cmdir_path

   starting_trna=$(get_setting_or_default "first_trna" "F")
   get_mito_hifi_recs_cm_anno >tmp_cm_anno
   ${script_dir}/right_neighbor_matrix.sh tmp_cm_anno $starting_trna >cm_anno_right_neighbor.matrix
   rm tmp_cm_anno

   msglog_module "cm_anno_right_neighbor.matrix created using ${BBlue}${starting_trna}${NC} as the first tRNA."

   # create an empty file that shows the cm_anno order
   > $(head -1 cm_anno_right_neighbor.matrix | sed "s/\t/__/g" | sed "s/ //g")__

   order=$(head -1 cm_anno_right_neighbor.matrix | sed "s/\t/ /g" | sed "s/^ *//" | sed "s/ *$//")
   update_setting "cm_anno_order" "$order"
   msglog_module "cm_anno rna order:" "$order"

   find_goosehairpin_bounding_trnas
}

function create_full_no_cr_seq_info {
   cd $cmdir_path

   local starting_trna=$(get_setting_or_default "first_trna" "F")
   local last_trna=$(get_setting_or_default "last_trna" "P")

   ${script_dir}/extract_full_no_cr_seqs.sh one_line_per_rec.cm_anno.srt mito_hifi_recs.cm_anno $starting_trna $last_trna > $cmdir_path/full_no_cr_seq.info
   if [ -s $cmdir_path/full_no_cr_seq.info ]; then
      msglog_module "full_no_cr_seq.info created created"
   elif [ -f $cmdir_path/full_no_cr_seq.info ]; then  # empty file delete it
      rm $cmdir_path/full_no_cr_seq.info
   fi
}

function create_matrix_and_anno_files {
   run_if_no_file create_trna_right_neighbor_matrix    $cmdir_path/trna_right_neighbor.matrix
   run_if_no_file create_cm_anno                       $cmdir_path/mito_hifi_recs.cm_anno
   run_if_no_file create_cm_anno_right_neighbor_matrix $cmdir_path/cm_anno_right_neighbor.matrix
   run_if_no_file create_cm_anno_one_liners            $cmdir_path/one_line_per_rec.cm_anno.srt
   run_if_no_file create_full_no_cr_seq_info           $cmdir_path/full_no_cr_seq.info
}

function count_and_record_gh_recs {
   local cm_anno=${cmdir_path}/mito_hifi_recs.cm_anno
   [ ! -s $cm_anno ] && return 1

   gh_recs=$(awk 'BEGIN{delete gh; rec=0}/^#/{rec++}NR==1{rec=1}/goose_hairpin/{gh[rec]++}END{print length(gh)}' $cm_anno)
   num_gh=$(awk 'BEGIN{num=0}/goose_hairpin/{nm++}END{print int(nm)}' $cm_anno)

   update_setting "gh_rec_count"     $gh_recs
   update_setting "gh_found_in_recs" $num_gh

   if [ ! $num_gh = 0 ]; then
      awk '/goose_hairpin/ {
            seqname[$1]++
            if (seqname[$1]==1) print $1
      }' $cm_anno > ${cmdir_path}/mito_hifi_seq_names_w_gh.txt
   fi
}

##########################################################################
#                                                                        #
#   call the functions if the files they create are not in $cmdir_path   #
#                                                                        #
##########################################################################

run_if_no_file search_12S_16S ${cmdir_path}/$cmout_16S
run_if_no_file do_OL_search ${cmdir_path}/OL.tbl
run_if_no_file tRNA_searches ${cmdir_path}/mito_hifi_recs.mitfi

if [ -s ${cmdir_path}/mito_hifi_recs.mitfi ]; then
   create_matrix_and_anno_files
   run_if_no_file count_and_record_gh_recs ${cmdir_path}/mito_hifi_seq_names_w_gh.txt
else
   msglog_module "error: ${cmdir_path}/mito_hifi_recs.mitfi was not found"
fi

#########################################################################

