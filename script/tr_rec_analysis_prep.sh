#!/bin/bash

# from cr_analysis dir:
# use the trf_output/trf_stats_report.tsv to and the *_CR_*.fasta file
# to create directories where the CR recs with tandem repeats of a specific copynum are copied
# for exsmple those CR recs with 2 copies of a tandem repeat (also will segregate the no repeat recs in a dir too)

: ' # these files would be created for copynum 2 by this script:
   tr_copynum_crs/2_copy_tr/cr_recs.info  # has the lines from trf_output/trf_stats_report.tsv for this copynum
   tr_copynum_crs/2_copy_tr/cr_recs.fa    # has the records in the above pulled from the *_CR_*.fasta file
'

: ' # this will also put an entry into tr_copynum_info.tsv file which looks like this:

#Copies	HiFi reads	cr_recs_info	cr_rec_fasta	cr_consensus	cr_anno
2	19	tr_copynum_crs/2_copy_tr/cr_recs.info	tr_copynum_crs/2_copy_tr/cr_recs.fa	tr_copynum_crs/2_copy_tr/cr_consensus.fa	tr_copynum_crs/2_copy_tr/cr_consensus.anno
0       1155    tr_copynum_crs/no_copy_tr/cr_recs.info	tr_copynum_crs/2_copy_tr/cr_recs.fa	tr_copynum_crs/no_copy_tr/cr_consensus.fa       tr_copynum_crs/no_copy_tr/cr_consensus.anno
'

# after the prep we will get a mafft msa consensus and annotate the consensus cr with OH, gh and the tandem repeats with the flanking tRNAs

: ' # trf_output/trf_stats_report.tsv has lines like
m64049_191213_040746/151258491/ccs_RC_Pro_CR_Phe 1898..2919  	1022	676--718 pos	44 len	61 score	2.0 copies	Consensus pattern 22 bp: AATGCTTGATAGACATTAAGTT
'

function msg { # write a msg to stderr
   >&2 echo -e "$@"
}

function copynum_categories {  # validate numbers in roght place in $stats file and returns number of different copynums rounded to ints, excludes no_copy
   awk 'BEGIN{FS="\t"; OFS="\t"; delete categ}
      /^#/{next}
      {cnum = int($6) }
      cnum > 1{ categ[ cnum ]++ }
      END{ for(c in categ) printf("%s ", c); printf("\n")}
   ' $stats
}
function num_categories {
   copynum_categories | awk '{print NF}'
}

function make_category_dirs {
   for cpnum in $(copynum_categories); do
      mkdir -p tr_copynum_crs/${cpnum}_copy_tr
   done
}
function need_no_copy_recs { # if there are tr cr_rec subset files and the no_copy_tr dir does not exist, then we need to create the no_copy_tr files
   ls tr_copynum_crs/*/cr_recs.fa >/dev/null 2>/dev/null && ! ls tr_copynum_crs/no_copy_tr >/dev/null 2>/dev/null
}
function add_info_files_to_copynum_dirs {
   awk '
      BEGIN{FS="\t"; OFS="\t"}
      /^#/{next}
      {
         cpnum = int($6); cpdir = "tr_copynum_crs/"cpnum"_copy_tr/"
         info_file = cpdir "cr_recs.info"
         print $0 > info_file
      }
   ' $stats
}
function pull_recs {
   local info_file=$1

   bawk '
      FNR==NR {
         pull[$name]++
         next
      }

      $name in pull {
         print ">"$name" "$comment
         print $seq
   } ' <(awk '{print ">"$0}' $info_file) $CR_recs
}
function pull_no_copy_recs {
   bawk '
      FNR==NR {  # if it in the stats file we do not want it, but take everything else
         no_pull[$name]++
         next
      }

      ! ($name in no_pull) {
         print ">"$name" "$comment
         print $seq
      } ' <(awk '/^#/{next}{print ">"$0}' $stats) $CR_recs
}
function make_no_copy {
   ! need_no_copy_recs && return

   local dir=tr_copynum_crs/no_copy_tr
   local recs=$dir/cr_recs.fa
   local info=$dir/cr_recs.info
   local prefix=$dir/cr_consensus_no_copy_tr

   mkdir -p $dir
   [ -d $dir ] && pull_no_copy_recs > $dir/cr_recs.fa

   local nrecs=$(numrecs $dir/cr_recs.fa)
   echo -e "0\t$nrecs\t$info\t$recs\t${prefix}.fasta\t${prefix}.anno\t${prefix}.compare" >>tr_copynum_info.tsv
}

function check_prereqs { # also sets CR_recs env var
   # presumes we are in the cr_analysis directory
   [ ! -d trf_output ] && return 1
   [ ! -s trf_output/trf_stats_report.tsv ] && return 2
   ls *_CR_*.fasta >/dev/null 2>/dev/null; [ $? != 0 ] && msg "No CR sequence file found" && return 3

   CR_recs=$(ls ???_CR_???.fasta | head -n 1)
   [ ! -s $CR_recs ] &&  msg "No CR sequence file found" && return 4

   return 0  # prereqs check out
}

function make_copynum_rec_subsets {
   stats=trf_output/trf_stats_report.tsv  # stats var needed in num_categories function
   nm_categs=$(num_categories); [ $nm_categs -lt 1 ] && return 5 # nothing to do

   make_category_dirs
   add_info_files_to_copynum_dirs

   # create the rec fasta file for each copynum category
   #for d in tr_copynum_crs/*_copy_tr; do
   for d in $(ls -dv tr_copynum_crs/*_copy_tr); do  # sort in numeric order 29Dec2022
      info=${d}/cr_recs.info
      recs=${d}/cr_recs.fa
      prefix=$d

      [ ! -s tr_copynum_info.tsv ] && echo -e "#Copies HiFi reads      cr_recs_info    cr_rec_fasta    cr_consensus_fa    cr_consensus_anno     cr_consensus_compare" >tr_copynum_info.tsv

      pull_recs $info >$recs
      nrecs=$(numrecs $recs)

      cpnum=$(echo $(basename $d) | awk '{print int($1); exit}')
      if [ ! "$(grep -s -c "$cpnum\s" tr_copynum_info.tsv)" = "1" ]; then # add line to tr_copynum_info.tsv
         echo -e "$cpnum\t$nrecs\t$info\t$recs\t${prefix}.fasta\t${prefix}.anno\t${prefix}.compare" >>tr_copynum_info.tsv
      fi
   done

   make_no_copy # if need be, that is if there are tr dirs already created
}

# presumes we are in the cr_analysis dir when called
check_prereqs && make_copynum_rec_subsets
