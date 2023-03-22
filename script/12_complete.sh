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
      touch -c -h ${wdir_path}/mitochondrion.gff
      touch -c -h ${wdir_path}/alternate_control_region_assemblies
      touch -c -h ${wdir_path}/mito_note.txt
   fi
}

# 18Mar2023
function Create_alternate_mitos {
   # see if there are any alternates
   num_alts=$(grep -m1 "^" $cr_dir/CR_alternates/CR_alt*.fa 2>/dev/null | wc -l)

   function seqs_same {  # takes 2 fasta files as args and is true if first seqs in each are the same
      local fa1=$1; local fa2=$2

      if [ ! -s $fa1 ] || [ ! -s $fa2 ]; then false; return; fi

      rslt=$(bawk '
         BEGIN{ rslt = "different" }
         FILENUM==1 && FNR==1 { len1=length($seq); seq1=$seq }
         FILENUM==2 && FNR==1 && len1==length($seq) && seq1==$seq { rslt = "same" }
         END { print rslt }
      ' $fa1 $fa2)

      [ $rslt == "same" ]  # set exit code based on the rslt string value
   }

   CR_in_chosen_mito=$cr_dir/CR_consensus_no_flanks.fasta
   [ ! -s $CR_in_chosen_mito ] && msglog_module $CR_in_chosen_mito "could not be found" && return

   if [ $num_alts -gt 0 ]; then
      for alt_cr in $cr_dir/CR_alternates/CR_alt*.fa; do
         if ! seqs_same $CR_in_chosen_mito $alt_cr; then
            Make_mito_with_alternate_CR $alt_cr
         else
            msglog_module "$(basename $alt_cr .fa) is identical to the Control Region used in the primary assembly. Skipped."
         fi
      done

      if [ -d $alternate_mito_path ]; then  # make link to alt dir in $wdir. alt dir only created if we make at least one cr_alt fasta
         local alternate_dir_name=alternate_control_region_assemblies
         local alternate_mito_path=$complete_dir/alternate_mitos
         [ ! -s $alternate_dir_name ] && ln -s $alternate_mito_path $alternate_dir_name
         [ ! -s $wdir/$alternate_dir_name ] && ln -s complete/alternate_mitos $wdir/$alternate_dir_name

         echo "alternates created" > $alternate_mito_path/alternates.created
         msglog ""
         touch -c -h mito_note.txt # we want notes file after the softlinks
      fi
   fi
}

function make_mito_CR_pre_post_seqs { # print seqs we need to glue around our new CR to files. do not use ending CR
   pre=$alt_dir/sequence_before_CR
   post=$alt_dir/sequence_after_CR
   bawk -v pre=$pre -v post=$post '
      NR==1 { CR_seq = $seq  }
      FILENUM==2 {
         cr_begin = index($seq, CR_seq)
         printf("%s", substr($seq, 1, cr_begin-1)) > pre
         printf("%s", substr($seq, cr_begin+length(CR_seq))) > post
      }
   ' $CR_in_chosen_mito $complete_dir/mitochondrion.fasta
}

function Make_mito_with_alternate_CR {
   alt_cr=$1; bn=$(basename $alt_cr)
   alt_dir=$complete_dir/alternate_mitos
   pre=$alt_dir/sequence_before_CR
   post=$alt_dir/sequence_after_CR

   mkdir_if_needed $alt_dir
   [ ! -f $alt_dir/sequence_before_CR ] && make_mito_CR_pre_post_seqs

   # ready to glue new CR into mito to create alternate
   new_name=$(bawk '{print "mito_w_" $name "_" $comment}' $alt_cr)

   (echo ">"$new_name
    fold -w 120 $pre
    bawk '{print "\n" $seq}' $alt_cr | fold -w 120
    fold -w 120 $post
   ) > ${alt_dir}/${new_name}.fasta

   if [ -s ${alt_dir}/${new_name}.fasta ]; then
      msglog_module mito assembly with alternate Control Region $bn created.
      Anno_mito_with_alternate_CR ${alt_dir}/${new_name}.fasta $alt_cr
   fi
}

function Anno_mito_with_alternate_CR { # put the alternate mito CR's annotation in place of the primary mito CR anno
   new_mito=$1; new_anno=$(replace_ext.sh $new_mito anno)
   mito_name=$(basename $new_anno .anno)
   alt_cr_anno=$(replace_ext.sh $2 anno) # where we get the CR anno

   cr_add=$(cawk '{print length($0); exit}' $complete_dir/alternate_mitos/sequence_before_CR)

   # replace the CR lines in $complete_dir/mitochondrion.anno with $alt_cr_anno
   cawk -t -v mito_name=$mito_name -v cr_add=$cr_add '
      FILENUM==1 {
         if (FNR==1) new_cr_len = $3
         $1 = mito_name; $2 += cr_add; $3 += cr_add
         new_cr_anno[++anno] = $0
      }

      FILENUM != 2 { next }
      $8 == "Control Region" {
         old_cr_len = $3 - $2 + 1
         loc_mod = new_cr_len - old_cr_len
         skip_till_gt = $3; skip = 1
         for (a=1; a <= anno; a++)
            print new_cr_anno[a]
      }

      /^#/ { print; next }
      { $1 = mito_name }
      skip == 0 { print; next }
      $2 > skip_till_gt { $2 += loc_mod; $3 += loc_mod; print }

   ' $alt_cr_anno $complete_dir/mitochondrion.anno > $new_anno
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

# 18Mar2023 create a directory for alternate mito fasats using any CR_alterantes not same as CR used in chosen assembly
run_if_no_file Create_alternate_mitos $complete_dir/alternate_mitos/alternates.created

# final messages for pipeline
msglog_module "links to mitochondrion.fasta, mitochondrion.anno and mitochondrion.gff created."
msglog_module "log for the pipeline is in $wdir/hfmt.log and settings used and found in $wdir/settings.tsv"
msglog_module "please read mito_notes.txt"
