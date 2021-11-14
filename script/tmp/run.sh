#!/bin/bash

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, script_dir, msglog, msglog_module, is_number functions and usage among other things

: '
1 init.sh  -- calls gettaxid.sh and is add_files.sh called with a hyphen for first arg
2 blast_to_mito.sh
  2b  make_mito_rec_cand_tsv.sh
3 pull_fofn_cand_recs.sh
4 select_mito_features.sh
5 blast_features.sh
6 split_linear_cr_recs.sh # to do
7 call megahit on linear recs to create mito cand w/o cr
8 selects subset of CR recs (within 1 SD?) for kalign and consensus calling
  8b call trf soemwehre around here
9 report on mito and CR heteroplasmy / repeats etc
10 annotate trn, lsu, ssu, genes
   trn with mitfi.sh
11 append consensus CR onto mito cand
'

: ' files created, etc
1 init.sh creates:
   files_to_search.fofn
   toptaxid
   taxidlist
   todo: other settings such as threads to use

2 blast_to_mito.sh (uses blastmax5.sh)
   creates dir: hifi_mito_matches
   blasts to create .tsv for the files in files_to_search.fofn
   calls make_mito_rec_cand_tsv.sh with the tsv files as args

2a make_mito_rec_cand_tsv.sh (called by blast_to_mito.sh)
   consolidates all tsvs into file: mito_rec_candidates.tsv
   creates mito coverage stat file: top_mito_species_match_coverage.bins (has ids for ones we want to download features)

3 pull_fofn_cand_recs.sh
  creates mito_rec_candidates.fasta using mito_rec_candidates.tsv to pull records from files in files_to_search.fofn
  creates dir: rec_db and builds blast db from mito_rec_candidates.fasta named mito_rec_candidates

4 select_mito_features.sh
  uses top_mito_species_match_coverage.bins to get list of IDs for which features are desired
  calls mito_analyze.py to download each into mito_feature_sequences.fasta

5 blast_features.sh (calls blastn with -task blastn)
  blasts mito_feature_sequences.fasta to db mito_rec_candidates in dir rec_db/
  output into mito_feature_match_to_cand_recs.tsv
  also create_dist_file trnF_trnP_distances.tsv
  and creates ControlRegion_length.stats from the dist file

6 split_recs_into_2_sets.sh: from mito_rec_candidates.fasta using mito_feature_match_to_cand_recs.tsv info, creates:
     recs_linear_split_nocr.fa
     recs_cr_w_trna_flanks.fa

hfmt.log used for logging information -- work in progress asis it all
'
