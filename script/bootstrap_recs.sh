#!/bin/bash

: ' bootstrap from the cand recs to search for more recs from the input using a preliminary build of the mitochondrion
   __1_assemble_w_megahit_using_mito_hifi_candidates_recs__
   __2_blast_megahit_best_rec_to_db_mito_to_orient_seq__
   __3_make_blastdb_of_seq_and_use_input_as_query_to_blast_this_db__
   __4_add_additional_found_recs_to_candidates__
   megahit_out
'

source $(dirname $(realpath $0))/shared.sh
[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

################################################################################
#                              utility functions                               #
################################################################################

function create_info_files {
   cd $bootstrap_dir
   > "__1_assemble_w_megahit_using_mito_hifi_candidates_recs__"; sleep 1
   > "__2_blast_megahit_best_rec_to_db_mito_to_orient_seq__"; sleep 1
   > "__3_make_blastdb_of_seq_and_use_input_as_query_to_blast_this_db__"; sleep 1
   > "__4_add_additional_found_recs_to_bootstrap_recs__"
   cd $curdir
}

function cat_input_files {  # use fq2fa which handles gzipped file
   fofn=$wdir_path/files_to_search.fofn
   msglog_module "Reading files listed in $fofn"

   while IFS= read -r hifi_input_file; do
      [ -z "$hifi_input_file" ] && continue; [ ! -s "$hifi_input_file" ] && continue

      msglog_module $hifi_input_file
      fq2fa $hifi_input_file
   done < $fofn
}

function mitodb_recs_not_in_bootstrap_recs {
   bawk '
      { nm = $name; sub("_RC", "", nm) }
      FILENUM==1 { ar[nm] = $name; next }
      ! (nm in ar) { print ">" $name " " $comment; print $seq }
   ' bootstrap_matching_recs.fasta ../mito_hifi_candidates.fasta | fold -w 120
}

function get_comparison_str {
   cawk  '
      ! /^>/ { next }
      {
         full_nm = substr($1,2)
         nm = full_nm; sub("_RC", "", nm)
         is_rc = !(full_nm==nm)
         qpct = $2; gsub("[^0-9]+", "", qpct); qpct = int(qpct)
      }
      FILENUM==1 { mitodb_rec_qpct[nm]    = qpct; mitodb_rec_rc[nm]    = is_rc }
      FILENUM==2 { bootstrap_rec_qpct[nm] = qpct; bootstrap_rec_rc[nm] = is_rc }

      END {
         for (r in mitodb_rec_qpct) {
            if (r in bootstrap_rec_qpct) {
               in_bootstrap++
               bootstrap_same_orientation_as_mitodb_rec += (bootstrap_rec_rc[r]==mitodb_rec_rc[r])
            }
            else not_in_boostrap++
         }
         printf("mitodb_recs_in_bootstrap_recs:\t%d\t%d\t%d\t%d",  in_bootstrap, not_in_boostrap, in_bootstrap - bootstrap_same_orientation_as_mitodb_rec, length(mitodb_rec_qpct))

         for (r in bootstrap_rec_qpct) {
            if (r in mitodb_rec_qpct) {
               in_mitodb++
               mitodb_same_orientation_as_bootstrap += (bootstrap_rec_rc[r]==mitodb_rec_rc[r])
            }
            else not_in_mitodb++
         }
         printf("\tbootstrap_recs_in_mitodb_recs\t%d\t%d\t%d\t%d\n", in_mitodb, not_in_mitodb, in_mitodb - mitodb_same_orientation_as_bootstrap, length(bootstrap_rec_qpct))
      }
   ' ../mito_hifi_candidates.fasta bootstrap_matching_recs.fasta
}

function set_comparison_vars { # expecting to be in bootstrap_dir
   # returns for example: mitodb_recs_in_bootstrap_recs:   211     0       0       211     bootstrap_recs_in_mitodb_recs   211     156     0       367
   # 1st number in set is number matching recs, 2nd is number missinf, 3rd is revcomp mismatches, 4th is total recs. then same four for the bootstrap fasta file

   cmpstr=$(get_comparison_str)

   mitodb_cand_recs=$(get_cmp_field 5)
   update_setting_if_changed "mitodb_cand_recs" $mitodb_cand_recs

   bootstrap_recs=$(get_cmp_field 10)
   update_setting_if_changed "bootstrap_recs" $bootstrap_recs

   mitodb_cand_recs_not_in_bootstrap_recs=$(get_cmp_field 3)
   update_setting_if_changed "mitodb_cand_recs_not_in_bootstrap_recs" $mitodb_cand_recs_not_in_bootstrap_recs
}

function get_cmp_field {
   local fnum=$1
   echo $cmpstr | awk -v fnum=$fnum '{print $fnum}'
}

################################################################################
#                                main functions                                #
################################################################################

function make_bootstrap_dir {
   if [ ! -d $bootstrap_dir ]; then
      mkdir $bootstrap_dir

      msglog ""
      msglog_module "Creating $(basename $bootstrap_dir) directory to search for additional mitochondrial reads in the input files"

      create_info_files
   else
      msg "$(basename $bootstrap_dir) directory already created"
   fi
}

function run_megahit {
   local threads=1
   local mito_fasta_recs_to_assemble=$wdir_path/mito_hifi_candidates.fasta

   cd $bootstrap_dir

   msglog_module "creating boostrap assembly with megahit from $mito_fasta_recs_to_assemble"
   msglog_module "megahit -t $threads -r $mito_fasta_recs_to_assemble"
   megahit -t $threads -r $mito_fasta_recs_to_assemble

   cd megahit_out

   # this should be name of best matching contig in the megahit.out directory
   ctg=$(grep "^>" final.contigs.fa | sed "s/multi=//" | sort -k3,3Vr | awk '{print substr($1,2);exit}')

   awk -v ctg="${ctg}" '
      $1 == ">"ctg {
         print $0
         prt=1
         next
      }
      prt && /^>/ {exit}
      prt { print }
   ' final.contigs.fa > megahit_best.fa
   msglog_module "megahit_best.fa mito file created"

   cd $curdir
}

function orient_megahit_rec { # leave as is or revcomp depending on best mito match
   cd $bootstrap_dir/megahit_out

   blastn -db mito -query megahit_best.fa -outfmt "6 std staxid stitle qlen qcovhsp qcovus" > blast_hits_for_revcomp_check.tsv

   # check $9 $10 to see if $9 > $10 and revcomp if so
   do_rc=$(awk '$9 > $10{print "RC"}{exit}' blast_hits_for_revcomp_check.tsv)

   if [ ! -z "$do_rc" ]; then # need to revcomp megahit_best.fa
      msglog_module "revcomp megahit_best.fa"
      bawk '{print ">"$name"_RC " $comment; print revcomp($seq)}' megahit_best.fa > bootstrap_rec.fa
   else # create softlink to megahit_best.fa
      ln -s megahit_best.fa bootstrap_rec.fa
   fi

   cd $curdir
}

function blast_bootstrap_rec {
   cd $bootstrap_dir/megahit_out

   # make blastdb for the fasta mito bootstrap record
   makeblastdb -dbtype nucl -parse_seqids -in bootstrap_rec.fa

   # blast against files named in files_to_search.fofn
   msglog_module "blastn -db bootstrap_rec.fa -query <(cat_input_files) -outfmt \"6 std staxid stitle qlen qcovhsp qcovus\" -num_threads $threads -mt_mode 1 > ../bootstrap_blast_hits_to_input.tsv"
   blastn -db bootstrap_rec.fa -query <(cat_input_files) -outfmt "6 std staxid stitle qlen qcovhsp qcovus" -num_threads $threads -mt_mode 1 > ../bootstrap_blast_hits_to_input.tsv

   cd $curdir
}

function extract_recs {
   cd $bootstrap_dir

   [ ! -s bootstrap_blast_hits_to_input.tsv ] && msglog_module "Could not find $bootstrap_dir/bootstrap_blast_hits_to_input.tsv" && cd $curdir && return 1

   msglog_module "Extracting records from input that have at least ${MIN_COV}% coverage to the bootstrap assembly"

   bawk -v MIN_COV=$MIN_COV  -v wdir=$wdir '
      FILENUM==1 {
         if ($NF >= MIN_COV) {
           extract_ar[$1] = int($NF)
           if ($11 > $12) rc_ar[$1]++
         }
         next
      }

      ! ($name in extract_ar) { next }  # if this is not one we want go to the next record

      # set qpct, suffix if RC, and gh_pos str addtl
      { qpct = extract_ar[$name]; suffix = "" }
      $name in rc_ar { revcomp($seq); suffix = "_RC" }
      { gh_pos = goosehairpin_pos(); addtl = (gh_pos) ? " gh:" gh_pos : "" }

      {print ">" $name suffix " QPct:" qpct addtl; print $seq; pulled++}

      function goosehairpin_pos() {
         return match($seq, /CCCCCCC[AGT][AGT]?[AGT]?CCCCCCC/))   # nawk does not support quantification notation /CCCCCCC[AGT]{1,3}CCCCCCC/) so use ? twice
      }
   ' <(prefix_gt bootstrap_blast_hits_to_input.tsv) <(cat_input_files) | fold -w 120 > bootstrap_matching_recs.fasta

   cd $curdir
}

function make_mito_hifi_recs_fasta {
   cd $bootstrap_dir

   set_comparison_vars  # sets: cmpstr mitodb_cand_recs bootstrap_recs mitodb_cand_recs_not_in_bootstrap_recs
   msglog_module $cmpstr

   if [ "$mitodb_cand_recs_not_in_bootstrap_recs" = "0" ]; then
      msglog_module "Using the bootstrap records as the mito_hifi_recs.fasta dataset"
      to_link=bootstrap_matching_recs.fasta
   else
      to_link=bootstrap_mitodb_rec_union.fasta
      msglog_module "adding $mitodb_cand_recs_not_in_bootstrap_recs recs in the mitodb set not found in the bootstrap recs to create $to_link"
      cat bootstrap_matching_recs.fasta <(mitodb_recs_not_in_bootstrap_recs) > $to_link
   fi

   cd ..
   [ -f mito_hifi_recs.fasta ] && rm mito_hifi_recs.fasta
   ln -s $(basename $bootstrap_dir)/$to_link mito_hifi_recs.fasta
   cd $curdir; cd $bootstrap_dir

   echo $cmpstr > bootstrap.done

   cd $curdir
   touch -h $wdir/mito_hifi_recs.fasta  # put this below the bootstrap_dir for lt
}


################################################################################
#                                 starts here                                  #
################################################################################

MIN_COV=50

curdir=$(pwd)

make_bootstrap_dir

run_if_no_file run_megahit          $bootstrap_dir/megahit_out/megahit_best.fa
run_if_no_file orient_megahit_rec   $bootstrap_dir/megahit_out/bootstrap_rec.fa
run_if_no_file blast_bootstrap_rec  $bootstrap_dir/bootstrap_blast_hits_to_input.tsv
run_if_no_file extract_recs         $bootstrap_dir/bootstrap_matching_recs.fasta

run_if_no_file make_mito_hifi_recs_fasta $bootstrap_dir/bootstrap.done

cd $curdir  # should be unnecessary but just in case
