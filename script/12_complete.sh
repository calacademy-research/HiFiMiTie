#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things

[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

###############################################################################
#                                functions                                    #
###############################################################################

function anno_file_name {
   echo $complete_dir/mitochondrion.anno
}
function gff_file_name {
   echo $complete_dir/mitochondrion.gff
}

function make_anno_symbol_file {
   local syms=$( awk '/^#/{next}BEGIN{printf("_")}{printf("_%s",$7)}END{printf("__\n")}' $(anno_file_name) )
   > ${complete_dir}/$syms

   cleaned_syms=$(echo $syms | sed "s/_/ /g" | sed -e "s/^ *//" -e "s/ *$//")
   update_setting_if_changed anno_order "$cleaned_syms"
}

function make_gff_file {
   [ ! -s $(anno_file_name) ] && return

   local anno=$(anno_file_name)
   local gff=$(gff_file_name)

   $script_dir/convert_hfmt_anno_to_gff.sh $anno >$gff
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

###############################################################################
#                            script starts here                               #
###############################################################################

source ${script_dir}/mito_note_file.sh

mkdir_if_needed $complete_dir

###############################################################################
# full megahit and msa anno files created: use one of them, no need for pcg   #
# anno done here previously since we have those -- 05Dec2022                  #
###############################################################################

chosen_asm=$(get_setting chosen_assembly)
chosen_anno=$(get_setting chosen_anno)

# remove link files if they exist
rm -f $wdir/mitochondrion.fasta; rm -f $wdir/mitochondrion.anno; rm -f $wdir/mitochondrion.gff

# copy the fasta assembly file
if [ -s "${wdir}/$chosen_asm" ]; then
   cp -aL ${wdir}/$chosen_asm $complete_dir/mitochondrion.fasta
   ln -s complete/mitochondrion.fasta $wdir_path/mitochondrion.fasta
else
   msglog_module "No $chosen_asm found, you'll need to look at the 2 assembly comparisons and decide."
   exit 1
fi

# copy the anno file previously made in step 9 and softlinked in $wdir
if [ -s "${wdir}/$chosen_anno" ]; then
   cp -aL ${wdir}/$chosen_anno $complete_dir/mitochondrion.anno
   ln -s complete/mitochondrion.anno $wdir_path/mitochondrion.anno

   # make a gff file from the anno file
   make_gff_file
   ln -s complete/mitochondrion.gff $wdir_path/mitochondrion.gff

   # at an empty file whose name shows the mito symbols
   make_anno_symbol_file
else
   msglog_module "No $chosen_anno found."
fi

# link to these from the original directory
make_source_dir_links_and_file

# final messages for pipeline
msglog_module "links to mitochondrion.fasta, mitochondrion.anno and mitochondrion.gff created."
msglog_module "log for the pipeline is in $wdir/hfmt.log and settings used and found in $wdir/settings.tsv"
msglog_module "please read mito_notes.txt"
