#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

###############################################################################
#                                functions                                    #
###############################################################################

function make_chosen_mito_blastdb {
   local name=mitochondrion.fasta

   mito_asm_db=$complete_dir/blastdb/$name

   pushd $complete_dir >/dev/null

   mkdir_if_needed blastdb; cd blastdb
   [ -s ${name}.nin ] && popd >/dev/null && return 0  ## blastdb already created

   makeblastdb -in ../$name -dbtype nucl -parse_seqids -out ./$name

   popd >/dev/null
}

function anno_file_name {
   echo $complete_dir/mitochondrion.anno
}
function gff_file_name {
   echo $complete_dir/mitochondrion.gff
}

function feat_output_file {
   echo $complete_dir/mito_feature_matches.tsv
}

function make_anno_symbol_file {
   local syms=$( awk '/^#/{next}BEGIN{printf("_")}{printf("_%s",$7)}END{printf("__\n")}' $(anno_file_name) )
   > ${complete_dir}/$syms

   cleaned_syms=$(echo $syms | sed "s/_/ /g" | sed -e "s/^ *//" -e "s/ *$//")
   update_setting_if_changed anno_order "$cleaned_syms"
}

function get_pcg_ids {
    AAsyms.sh $( awk '{printf(" %s", substr($1,1,index($1, "_")-1) )}END{print ""}' $(feat_output_file) ) \
  | grep Unknown | grep -v ^1[26]S | awk '{print "^"$1}'
}

function get_pcg_tsv_lines {
   grep -f <( get_pcg_ids ) $(feat_output_file)
}

function pcg_lines_in_mitfi_format {
   awk 'BEGIN{FS="\t"; OFS="\t"}
      {print $2, $9, $10, $12, $11, ".", substr($1,1,index($1, "_")-1), $1, "."}
   ' <(get_pcg_tsv_lines)
}

function insert_pcg_into_cm_anno {
   awk 'BEGIN{FS="\t"; OFS="\t"}
      FNR==NR {
         ar[++pcg] = $0
         start[pcg] = $2
         next
      }

      FNR==1{ix=1; next_pcg_start = start[ix]}
      /^#/{ print; next }

      {
         begpos = $2; endpos = $3
         if ( begpos < next_pcg_start )
            print
         else {
            while ( begpos > next_pcg_start ) {
               print ar[ix]
               ix++
               next_pcg_start = (ix <= pcg) ? start[ix]: 1000000
            }
            print $0
         }
      }
   ' <(pcg_lines_in_mitfi_format) $cm_anno
}

function add_anno_link {
   local link_name=${wdir}/$(basename $(anno_file_name))
   [ ! -s $link_name ] && [ -s $(anno_file_name) ] && ln -s complete/$(basename $(anno_file_name)) $link_name
}

function add_gff_link {
   local gff_base=$(basename $(gff_file_name))
   local link_name=${wdir}/$gff_base

   [ ! -s $link_name ] && [ -s $(gff_file_name) ] && ln -s complete/$gff_base $link_name
}

function make_gff_file {
   [ ! -s $(anno_file_name) ] && return

   local anno=$(anno_file_name)
   local gff=$(gff_file_name)

   convert_hfmt_anno_to_gff.sh $anno >$gff
}

function create_anno_file {
   cm_anno=${wdir}/mito_msa.cm_anno
   [ ! -s $cm_anno ] && msglog_module "could not find $cm_anno" && return 1

   awk ' # add lengths and pos from last as extra fields
      BEGIN{FS="\t"; OFS="\t"}
      /^#/{ if (FNR==1) { print $0, "length", "dist_from_last" } else print; next }

      { n++ }
      { gh = $0 ~ "goose_hairpin" }

      {
         len = $3 - $2 + 1
         dist_from_lst = (n==1) ? (0) : ($2 - lst_end)
         if(gh) dist_from_lst = $2 - lst_beg
         print $0, len, dist_from_lst

         lst_beg = $2
         if (!gh) lst_end = $3
   } ' <(insert_pcg_into_cm_anno) > $(anno_file_name)

   add_anno_link

   make_gff_file
   add_gff_link

   make_anno_symbol_file

   msglog_module "$(basename $(anno_file_name)) created"
}

function blast_features_to_chosen_mito {

   evalue=1e-10
   query_file=$(feature_sequence_file)
   output_file=$(feat_output_file)
   db=$mito_asm_db

   # blast features against candidate rec db and output in position sorted order
   msglog_module "blastn of features against $(basename $mito_asm_db) for protein coding gene annotation"
   msglog_module blastn -db $db -query $(basename $query_file) -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -task blastn -evalue $evalue "|" sort -k2,2V -k9,9n -k12,12nr
   blastn -db $(realpath $db) -query $query_file -outfmt "6 std staxid stitle qlen qcovhsp qcovus" -task blastn -evalue $evalue \
      | sort -k2,2V -k9,9n -k12,12nr  \
      | bioawk_cas -t '{und=index($1,"_");feat=substr($1,1,und-1)}feat_lst==feat{next}{feat_lst=feat; print}' \
    > $output_file
}

function handle_msa_choice {
   local msa=$wdir/mito_msa.fasta
   local megahit=$wdir/mito_megahit.fasta
   local name=mitochondrion.fasta
   local mito=$complete_dir/$name

   if [ ! -s $mito ]; then

      if [ ! -s $msa ] && [ -s $megahit ]; then  # missing msa but we have a megahit: this needs more checks
         chosen=$megahit
         msglog_module "megahit consensus mitochondrion chosen as best representative"
      else
         chosen=$msa
         msglog_module "msa consensus mitochondrion chosen as best representative"
      fi

      [ ! -s $chosen ] && msglog_module "$chosen could not be found" && return 1

      mkdir_if_needed $complete_dir
      cp -a -H $chosen $mito && msglog_module "$name created"
   else
      msg "$name already created"
   fi

   [ ! -s $wdir_path/$name ] && [ -s $mito ] && ln -s complete/$name $wdir_path/$name

   produce_annotations
}

function produce_annotations {
   make_chosen_mito_blastdb
   run_if_no_file blast_features_to_chosen_mito $(feat_output_file)
   run_if_no_file create_anno_file $(anno_file_name)

   make_source_dir_links_and_file
   touch_wdir_links # move links underneath complete dir for lt listing
}

# everything else is in our hfmt working dir, but we will put softlinks
# and one mito_note file in the parent dir of it
function make_source_dir_links_and_file {
   local name=mitochondrion.fasta
   local mito=$complete_dir/$name
   local anno=mitochondrion.anno
   local anno_path=$complete_dir/$anno
   local gff=mitochondrion.gff
   local gff_path=$complete_dir/$gff
   local note_file=mito_note.txt

  [ ! -s $name ] && [ -s $mito ]      && ln -s $mito
  [ ! -s $anno ] && [ -s $anno_path ] && ln -s $anno_path
  [ ! -s $gff ]  && [ -s $gff_path ]  && ln -s $gff_path

  [ ! -s $note_file ] && basic_note > $note_file
}

function touch_wdir_links { # touch them unless we have finished the run, we do not want to keep changing them
   if [ -z "$(get_setting run_time)" ]; then
      touch -c -h ${wdir_path}/mitochondrion.fasta
      touch -c -h ${wdir_path}/mitochondrion.anno
   fi
}

# set assemble_CR_result based on the annos, preferring msa if it is within 60% megahit
function compare_annos {
   local msa=$wdir/mito_msa.cm_anno
   local megahit=$wdir/mito_megahit.cm_anno

   local msa_score=0
   [ -s $msa ] && msa_score=$(awk 'BEGIN{score=0.0}NR>1{score += $4}END{print score}' $msa)

   local mega_score=0
   [ -s $mega ] && mega_score=$(awk 'BEGIN{score=0.0}NR>1{score += $4}END{print score}' $megahit)


   # get a phrase repsenting the score comparison
   cmp_to_mega=$(echo $msa_score $mega_score | awk '
            function abs(a){return (a<0) ? -a : a}

            $1 == $2{print "equal to"; exit}
            $1 > $2{print "greater than"; exit}
            abs($1-$2) <= 60{print "near"; exit}
            {print "lt"}
        ')

   if [ "$cmp_to_mega" = "lt" ]; then
      msglog_module megahit assembly has cm_anno score of $mega_score compared to msa assembly score $msa_score.
      update_setting_if_changed assemble_CR_result "megahit"
   else
      msglog_module msa assembly cm_anno score of $msa_score is $cmp_to_mega megahit score $mega_score.
      update_setting_if_changed assemble_CR_result "msa"
   fi
}

###############################################################################

source ${script_dir}/mito_note_file.sh

mkdir_if_needed $complete_dir

result=$(get_setting assemble_CR_result)
[ -z $result ] && compare_annos && result=$(get_setting assemble_CR_result)

if [ -z "$result" ]; then
   msglog_module "No result found from last step, you'll need to look at the 2 assembly comparisons and decide."
elif [ "$result" = "msa" ] || [ "$result" = "megahit" ]; then
   handle_msa_choice
else
   msglog_module "Result \"$result\" found from last step is not recognized, you'll need to look at the 2 assembly comparisons and decide."
fi
