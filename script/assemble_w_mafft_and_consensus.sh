#!/bin/bash

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, wdir_path, script_dir, msglog, msglog_module, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

function run_mafft_and_consensus {
   fasta=$1
   names_from_fasta $fasta

   msglog_module "running mafft --auto $(basename $fasta) >$msa_out 2>$msa_log"
   mafft --auto $fasta > $msa_out_path 2> $msa_log_path

   [ ! -s $msa_out_path ] && msglog_module "error: empty alignment file $msa_out_path" && return 1

   msglog_module "consensus_from_fasta_alignment.sh $msa_out >$consensus_fa"
   consensus_from_fasta_alignment.sh $msa_out_path >$consensus_fa_path
}

function make_reverse {
   fasta=$1
   set_fastx_basename $fasta # sets fastx_basename var

   reversed_fasta=reversed_${fastx_basename}.fa
   reversed_fasta_path=$msa_dir/$reversed_fasta

   [ -s $reversed_fasta_path ] && msglog_module "$reversed_fasta already created" && return  # already exists

   bawk ' {
      print ">"$name"_reversed " $comment
      print reverse($seq)
   } ' $fasta >$reversed_fasta_path
}

# made for beg_to_Pro.fasta when things line up at the end not the beginning of each sequence
# reverse the fasta sequences, run mafft on this,
# then call consensus on the reversed sequence msa from mafft
# finally reverse the consensus so the everything is in the original order
function run_reverse_mafft_and_consensus {
   fasta=$1
   msglog_module running multiple sequence alignment on reversed sequences of $(basename $fasta)

   make_reverse $fasta

   # make names with reversed in them
   names_from_fasta $reversed_fasta

   msglog_module "running mafft --auto $reversed_fasta >$msa_out 2>$msa_log"
   [ ! -s $msa_out_path ] && mafft --auto $reversed_fasta_path >$msa_out_path 2>$msa_log_path

   [ ! -s $msa_out_path ] && msglog_module "error: empty alignment file $msa_out_path" && return 1

   consensus_reverse_out=$consensus_fa_path
   msglog_module "consensus_from_fasta_alignment.sh $msa_out >$consensus_fa"
   consensus_from_fasta_alignment.sh $msa_out_path >$consensus_reverse_out

   names_from_fasta $1  # $fasta may have changed do not use it
   msglog_module "re-reversing consensus file to make forward version $consensus_fa"
   bawk '{print ">"$name" "length($seq); print reverse($seq)}' $consensus_reverse_out >$consensus_fa_path
}

function make_dirs {
  mkdir_if_needed $alignasm_dir
  mkdir_if_needed $msa_dir
  mkdir_if_needed $consensus_dir
  mkdir_if_needed $aligncm_dir
}

function names_from_fasta {
   fasta=$1
   set_fastx_basename $fasta # sets fastx_basename var

   msa_out=${fastx_basename}.mafft.fa
   msa_out_path=$msa_dir/$msa_out

   msa_log=${fastx_basename}.mafft.log
   msa_log_path=$msa_dir/$msa_log

   consensus_fa=${fastx_basename}.consensus.fa
   consensus_fa_path=${consensus_dir}/$consensus_fa
}

function set_seq_filenames {
   # get the 3 sequence file paths
   fasta_trna_to_end=$(ls $splitseq_dir/*_end.fasta)
   fasta_beg_to_trna=$(ls $splitseq_dir/beg_*.fasta)
   fasta_CR_btw_trnas=$(ls $splitseq_dir/[A-Z]*_CR_*.fasta)

   if [ -z "$fasta_trna_to_end" ] || [ -z "$fasta_beg_to_trna" ] || [ -z "$fasta_CR_btw_trnas" ]; then
      false
   else
      true
   fi

   consensus_trna_to_end=$(basename $fasta_trna_to_end .fasta).consensus.fa
   path_consensus_trna_to_end=${consensus_dir}/$consensus_trna_to_end

   consensus_beg_to_trna=$(basename $fasta_beg_to_trna .fasta).consensus.fa
   path_consensus_beg_to_trna=${consensus_dir}/$consensus_beg_to_trna

   consensus_CR_btw_trnas=$(basename $fasta_CR_btw_trnas .fasta).consensus.fa
   path_consensus_CR_btw_trnas=${consensus_dir}/$consensus_CR_btw_trnas

   non_cr_consensus=non_cr_consensus.fasta
   path_non_cr_consensus=${consensus_dir}/$non_cr_consensus

   msa_no_cr_mito_fa=mito_msa_no_cr.fasta
   msa_no_cr_mito_fa_path=${alignasm_dir}/$msa_no_cr_mito_fa

   msa_mito_fa=mito_msa.fasta
   msa_mito_fa_path=${alignasm_dir}/$msa_mito_fa
}

# DEPRECATE: this can have issues is overlaps replace with function following
function make_consensus_of_two_non-cr_consensus_outputs {
   [ -s $path_non_cr_consensus ] && run_if_no_file noop $path_non_cr_consensus && return

   msglog_module "$non_cr_consensus from $consensus_trna_to_end & $consensus_beg_to_trna"
   [ ! -s $path_consensus_trna_to_end ] && msglog_module "$consensus_trna_to_end could not be found" && return 1
   [ ! -s $path_consensus_beg_to_trna ] && msglog_module "$consensus_beg_to_trna could not be found" && return 1

   names_from_fasta $path_non_cr_consensus
   msglog_module "running mafft --localpair <(cat $consensus_trna_to_end $consensus_beg_to_trna) >$msa_out 2>$msa_log"
   mafft --localpair --maxiterate 1000 <(cat $path_consensus_trna_to_end $path_consensus_beg_to_trna) > $msa_out_path 2> $msa_log_path

   [ ! -s $msa_out_path ] && msglog_module "$msa_out_path empty or not found. no $non_cr_consensus created" && return 1

   msglog_module "consensus_from_fasta_alignment.sh $msa_out >$non_cr_consensus"
   consensus_from_fasta_alignment.sh $msa_out_path >$path_non_cr_consensus
}

# 05Sep2022 this replaces above and makes mitfi output for each and blends the 2 based on those settings
function merge_two_non-cr_consensus_fa {
   [ -s $path_non_cr_consensus ] && run_if_no_file noop $path_non_cr_consensus && return

   msglog_module "$non_cr_consensus from $consensus_trna_to_end & $consensus_beg_to_trna"
   [ ! -s $path_consensus_trna_to_end ] && msglog_module "$consensus_trna_to_end could not be found" && return 1
   [ ! -s $path_consensus_beg_to_trna ] && msglog_module "$consensus_beg_to_trna could not be found" && return 1

   names_from_fasta $path_non_cr_consensus
   msglog_module "Merging $consensus_trna_to_end $consensus_beg_to_trna"

   local code_setting=$(get_setting_or_default "code" "5")
   msglog_module "merge_two_consensus_fa.sh $path_consensus_trna_to_end $path_consensus_beg_to_trna $code_setting >$path_non_cr_consensus"

   $script_dir/merge_two_consensus_fa.sh $path_consensus_trna_to_end $path_consensus_beg_to_trna $code_setting >$path_non_cr_consensus
}

function mitfi_fasta {
   [ ! -s "$1" ] && msglog_module $1 not found or empty && return 1

   set_fastx_basename $1 # sets fastx_basename var
   mitfi_file=${aligncm_dir}/${fastx_basename}.mitfi

   run_if_no_file noop $mitfi_file
   [ -s $mitfi_file ] && return 0  # already made the file

   local code_setting=$(get_setting_or_default "code" "5")
   local code_option="-code $code_setting"
   local any_addtl_options="$2"

   msglog_module running mitfi $code_option $any_addtl_options analysis on $(basename $1)
   mitfi.sh $code_option $any_addtl_options $1 >$mitfi_file

   [ -s $mitfi_file ] && msglog_module $(basename $mitfi_file) created with $(nocomment $mitfi_file | wc -l) entries
}

function do_OL_search {
   local mito_fasta=$(realpath $1)
   set_fastx_basename $mito_fasta
   local pfx=${fastx_basename}
   local OL_tbl=$aligncm_dir/${pfx}_OL.tbl

   run_if_no_file noop $OL_tbl
   [ -s $OL_tbl ] && return 0 # already made the file

   search_OL $mito_fasta $aligncm_dir
   [ -s $aligncm_dir/OL.tbl ] && mv $aligncm_dir/OL.tbl $OL_tbl
   [ -s $OL_tbl ] && msglog_module "OL.tbl created in $aligncm_dir for $(basename $mito_fasta)"
}

# search for 12S and 16S rrna sequences in the fasta file of arg 1
function search_12S_16S {
   local mito_fasta=$(realpath $1)

   set_fastx_basename $mito_fasta
   local pfx=${fastx_basename}

   tbl_12S=${pfx}_rrnS.tbl
   tbl_16S=${pfx}_rrnL.tbl

   run_if_no_file noop $aligncm_dir/$tbl_16S
   [ -s $aligncm_dir/$tbl_16S ] && return 0  # already made the file

   $(get_script rrna_cmsearch) $mito_fasta $pfx

   [ -s ${wdir}/$tbl_12S ] && mv ${wdir}/${pfx}_rrnS.* $aligncm_dir/
   [ -s ${wdir}/$tbl_16S ] && mv ${wdir}/${pfx}_rrnL.* $aligncm_dir/

   [ -s $aligncm_dir/$tbl_12S ] && [ -s $aligncm_dir/$tbl_16S ] && msglog_module "$tbl_12S $tbl_16S created in $aligncm_dir for $(basename $mito_fasta)"
}

function blast_OH_to_fasta {
   qry=$1
   evalue=.00001
   output_file=$2
   OH_db=${data_dir}/OrgRepl/OH.fas

   msglog_module blastn -db OH.fas -query $qry -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads
   blastn -db $OH_db -query $qry -outfmt  "6 std staxid stitle qlen qcovhsp qcovus"  -max_target_seqs 5 -subject_besthit -evalue $evalue -num_threads $threads >$output_file
}

function anno_fasta { # mitfi file and one with gh and 12S and 16S file added created
   local fasta=$1
   [ ! -s "$fasta" ] && msglog_module $fasta not found or empty && return 1

   mitfi_fasta $fasta
   do_OL_search $fasta
   search_12S_16S $fasta

   set_fastx_basename $fasta
   local pfx=$aligncm_dir/${fastx_basename}

   local OH_blast_path=${pfx}_OH.tsv
   blast_OH_to_fasta $fasta $OH_blast_path

   ${script_dir}/add_goosehairpin_to_mitfi.sh ${pfx}.mitfi $fasta |
   ${script_dir}/add_rrna_to_mitfi.sh ${pfx}_rrnS.tbl - |
   ${script_dir}/add_rrna_to_mitfi.sh ${pfx}_rrnL.tbl - |
   ${script_dir}/add_rrna_to_mitfi.sh ${pfx}_OL.tbl   - |
   ${script_dir}/add_OH_to_mitfi.sh $OH_blast_path    - |
   ${script_dir}/add_CR_to_mitfi.sh ${wdir_path}/settings.tsv - $(asmlen $fasta) > ${pfx}.cm_anno

   [ -s ${pfx}.cm_anno ] && msglog_module $(basename ${pfx}.cm_anno) created with rrnS, rrnL, cr and any goose hairpin added to mitfi results

   if [ -s ${pfx}.cm_anno ]; then  # 01Dec2022 add genes to anno
         if [ ! -s ${pfx}.anno ]; then
            ${script_dir}/mito_pcg_anno.sh $fasta ${pfx}.cm_anno | add_header_to_anno > ${pfx}.anno
            [ -s ${pfx}.anno ] && msglog_module "$(basename ${pfx}.anno ) with genes added to ${pfx}.cm_anno created for $(basename $fasta)"
         else
            msg $(basename ${pfx}.anno) already created
         fi
      fi
}

function make_no_cr_mito_fasta  {
   # make the msa_no_cr_mito.fasta, which is cr succ thru cr prev_trna
   # by clipping off any sequence before the cr_succ_trna or after the cr_prev_trna

   cr_succ_trna=$(get_setting_or_default gh_succ_trna F); tb=$(three_letter_AA $cr_succ_trna)
   cr_prev_trna=$(get_setting_or_default gh_prev_trna E); te=$(three_letter_AA $cr_prev_trna)

   get_mitfi_entry_info $path_no_cr_mitfi $cr_succ_trna
   beg_succ=$tbeg; end_succ=$tend; entry_succ=$tentry

   get_mitfi_entry_info $path_no_cr_mitfi $cr_prev_trna
   beg_prev=$tbeg; end_prev=$tend; entry_prev=$tentry

   seqlen=$(bawk '{print length($seq);exit}' $path_non_cr_consensus)

   [ ! $entry_succ = "1" ]  && msglog_module "WARNING: expected $cr_succ_trna to be first trna. reported as number $entry_succ"
   [ ! $entry_prev = "22" ] && msglog_module "WARNING: expected $cr_prev_trna to be 22nd trna. reported as number $entry_prev"

   let remove_at_beg=beg_succ-1
   let remove_at_end=seqlen-end_prev
   let newlen=seqlen-remove_at_beg-remove_at_end
   msglog_module removing $remove_at_beg from beginning and removing $remove_at_end from the end to create ${newlen}bp $msa_no_cr_mito_fa

   bawk -v begremove=$remove_at_beg -v endremove=$remove_at_end -v tb=$tb -v te=$te '{
      newseq = substr($seq, begremove+1)
      slen = length(newseq)
      if(endremove) newseq = substr(newseq, 1, slen-endremove)
      printf(">mitochondrial_sequence_from_%s_to_%s_excluding_CR length %d\n", tb, te, length(newseq))
      print newseq
   }' $path_non_cr_consensus >$msa_no_cr_mito_fa_path
   [ -s $msa_no_cr_mito_fa_path ] && msglog_module $msa_no_cr_mito_fa created
}

function make_consensus_mito {
   cr_succ_trna=$(get_setting_or_default gh_succ_trna F); tb=$(three_letter_AA $cr_succ_trna)
   cr_prev_trna=$(get_setting_or_default gh_prev_trna E); te=$(three_letter_AA $cr_prev_trna)

   get_mitfi_entry_info $path_mitfi_msa_no_cr_mito $cr_prev_trna
   no_cr_prev_beg=$tbeg; no_cr_prev_end=$tend; no_cr_prev_entry=$tentry
   let no_cr_prev_len=no_cr_prev_end-no_cr_prev_beg+1

   get_mitfi_entry_info $path_mitfi_CR_btw_trnas $cr_prev_trna
   cr_prev_beg=$tbeg; cr_prev_end=$tend; cr_prev_entry=$tentry
   let cr_prev_len=cr_prev_end-cr_prev_beg+1

   get_mitfi_entry_info $path_mitfi_CR_btw_trnas $cr_succ_trna
   cr_succ_beg=$tbeg

   [ ! $cr_prev_entry = "1" ]  && msglog_module "WARNING: expected $cr_prev_trna to be first trna of the CR with flanking consesnus sequence. reported as number $entry_succ"
   [ ! $no_cr_prev_entry = "22" ] && msglog_module "WARNING: expected $cr_prev_trna to be 22nd trna of $msa_no_cr_mito_fa. reported as number $entry_prev"

   # extract just CR from the consensus by removing the trna flank sequence
   # if cr_succ_trna is not the same as the first_trna setting we will need to reflow to put this at the beginning
   bawk -v prev_end=$cr_prev_end -v succ_beg=$cr_succ_beg '
    FNR==NR { print ">msa_consensus_mitochondrion\n" $seq }
    FNR!=NR {
      newseq=substr($seq, 1, succ_beg-1)
      newseq=substr(newseq, prev_end+1)
      print newseq
    }' $msa_no_cr_mito_fa_path $path_consensus_CR_btw_trnas >$msa_mito_fa_path
    [ -s $msa_mito_fa_path ] && msglog_module $msa_mito_fa created
}

function reflow_to_first_trna {
   cr_succ_trna=$(get_setting_or_default gh_succ_trna F); tb=$(three_letter_AA $cr_succ_trna)
   first_trna=$(get_setting_or_default first_trna F); tf=$(three_letter_AA $first_trna)

   [ $cr_succ_trna != $first_trna ] && msglog_module "Control Region followed by $tb ($cr_succ_trna) but first_trna set to $tf, file will be reflowed to begin with $tf ($first_trna)"

   # get location of first_trna from the mito_no_cr mitfi file
   # use this location to do the reflow into a temporary file if the location is not 1
   # rename the temp file to the full mito msa consensus fasta
   get_mitfi_entry_info $path_mitfi_msa_no_cr_mito $first_trna
   first_trna_beg_pos="$tbeg"; first_trna_end_pos="$tend"; first_trna_entry="$tentry"

   [ -z $first_trna_beg_pos ] && msglog_module "Error: could not find first_trna $first_trna in $path_mitfi_msa_no_cr_mito" && return 1
   [ $first_trna_beg_pos = 1 ] && return 0

   reflow_temp=${msa_mito_fa_path}.reflow
   bawk -v first_trna_beg_pos=$first_trna_beg_pos '{
      first_trna_to_end = substr($seq, first_trna_beg_pos)
      bases_before_first_trna = substr($seq, 1, first_trna_beg_pos-1)
      newseq = first_trna_to_end bases_before_first_trna
      print ">"$name " length " length(newseq)
      print newseq
   }' $msa_mito_fa_path > $reflow_temp
   mv $reflow_temp $msa_mito_fa_path
   msglog_module "$msa_mito_fa reflowed to begin with $tf ($first_trna)"
}

function check_for_tandem_repeats {
   local trf_dir=$msa_trf_dir
   local fa_file=../$msa_mito_fa

   [ ! -d $trf_dir ] && mkdir $trf_dir; local retcode=$?
   [ ! -d $trf_dir ] && msglog_module "Problem creating $trf_dir [$retcode]" && return $retcode

   pushd $trf_dir >/dev/null
   run_trf.sh $fa_file
   echo "" ## trf needs final newline
   popd >/dev/null

   add_repeats_to_anno
}
function get_formatted_trf_line {  # depends on trf_stats being set and a file. check before calling
   head -n 1 $trf_stats |
   awk  'BEGIN{OFS="\t"}
      {sub("--","\t",$3)
      printf("%s\t%s\t%s\t.\t.\trepeat\tRepeat Region\t+\t%s copies %s bp consensus: %s\n", $1, $3, $7, $9,$13,$NF)}
   '
}
function add_repeats_to_anno {
   anno=$aligncm_dir/mito_msa.anno
   trf_stats=$msa_trf_dir/trf_stats_report.tsv

   trf_anno=${anno}_w_trf

   [ ! -s $anno ] && return 3; [ ! -s $trf_stats ] && return 4

   if [ ! -z "$(get_formatted_trf_line)" ]; then
      cawk -t '
         FILENUM==1 && FNR==1{trf_start = $2; trf_end = $3; trf_line = $0}

         FILENUM==2 && /^#/ {print; next}
         FILENUM==2 {
            if ($2 > trf_start && ! trf_prtd) {
               print trf_line
               trf_prtd = 1
            }
            print
         }
         END { if (! trf_prtd) print trf_line }
      ' <(get_formatted_trf_line) $anno > $trf_anno

      anno_size=$(stat -c %s $anno)
      trf_anno_size=$(stat -c %s $trf_anno)
      if (( trf_anno_size > anno_size )); then
         mv $trf_anno $anno
         msglog_module "tandem repeat info added to mito_msa.anno"
      else
         rm $trf_anno
      fi
   fi
}

# 11Aug2022 this currently uses gh_succ_trna and gh_prev_trna when if missing defaults as F P
# gh is often not found so this is flawed. what we will do here is set these from more rational presumptions if they
# have not already been set. also later on we may use OH setting alone or in tandem with gh
# but always need to handle possibility than only first_trna and last_trna is set (and we may expand to all 12S for last_trna later)
function clean_up_cr_settings {
   gh_succ=$(get_setting gh_succ_trna)
   gh_prev=$(get_setting gh_prev_trna)

   [ -z $gh_succ ] && update_setting gh_succ_trna $(get_setting_or_default first_trna F)
   [ -z $gh_prev ] && update_setting gh_prev_trna $(get_setting_or_default last_trna P)
}

###############################################################################
#                                do the work                                  #
###############################################################################

make_dirs
set_seq_filenames
clean_up_cr_settings # 11Aug2022

run_if_no_file noop $path_consensus_trna_to_end
[ ! -s $path_consensus_trna_to_end ] && run_mafft_and_consensus $fasta_trna_to_end

run_if_no_file noop $path_consensus_beg_to_trna
[ ! -s  $path_consensus_beg_to_trna ] && run_reverse_mafft_and_consensus $fasta_beg_to_trna

run_if_no_file noop $path_consensus_CR_btw_trnas
[ ! -s $path_consensus_CR_btw_trnas ] && run_mafft_and_consensus $fasta_CR_btw_trnas

# take the 2 consensus sequences and run mafft in sensitive mode
## make_consensus_of_two_non-cr_consensus_outputs
# 05Sep2022 replace above with merge which relies on mitfi of the 2 consensus sequences
merge_two_non-cr_consensus_fa

# clean up the non-cr_consensus file to create the final msa_no_cr_mito.fasta, this requires the mitfi file for the consensus
mitfi_fasta $path_non_cr_consensus "-onlycutoff"
path_no_cr_mitfi=$mitfi_file

run_if_no_file make_no_cr_mito_fasta $msa_no_cr_mito_fa_path
anno_fasta $msa_no_cr_mito_fa_path
path_mitfi_msa_no_cr_mito=$mitfi_file

# get the trna info for the CR consensus flanked by the trna
mitfi_fasta $path_consensus_CR_btw_trnas "-onlycutoff"
path_mitfi_CR_btw_trnas=$mitfi_file

# take the msa_no_cr_mito_fa and append the consensus CR
run_if_no_file make_consensus_mito $msa_mito_fa_path
set_fastx_basename $msa_mito_fa_path
path_mitfi_msa_mito=${aligncm_dir}/${fastx_basename}.mitfi

# if cr_succ_trna is not the same as the first_trna setting we will need to reflow to put this at the beginning
[ ! -s $path_mitfi_msa_mito ] && reflow_to_first_trna
anno_fasta $msa_mito_fa_path

# search for tandem repeats in the mito, specifically expect them if anywhere in the Control Region(s)
msa_trf_dir=${alignasm_dir}/trf_output
run_if_no_file check_for_tandem_repeats ${msa_trf_dir}/trf_stats.txt
