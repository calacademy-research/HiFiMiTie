#!/bin/bash

# separate file so we can show extensive usage / help info if we want to, without getting in way of main script

cmd=$1

if [ -z $cmd ]; then
    echo -e '
    HiFiMiTie version 0.001

    Usage: hifimitie <command> [options]

    Commands:
      -- Setup
         init              create working directory and optionally set taxid and HiFi file(s)
           gettaxid        search all mitogenomes or only those associated with a certain taxid (optionally view taxonomy ids based on term)
           addfiles        add some or all of the HiFi files with which to begin the search
           clearfiles      remove the files already in the fofn file list

      -- Execution
         run               run the default pipeline of the commands below, specifying number of threads to use with -t <num_threads>, defaults to 32

         blast_to_mito     blast the mito assemblies associated with the taxids you have chosen against the PacBio HiFi files

         pull_fofn_cand_recs  create file with candidate mito recs based on the blast tsv file pulled from the HiFi files

         select_mito_features  for the most closely matching mito assemblies, pull tRNA, rRNA, gene, CR features from GenBank into records
         blast_features        blast the features against the candidate records using -task blastn instead of default megablast

         rna_search            rrna_cmsearch of SSU, LSU rrna. multi_mitfi for trna, goose hairpin info added into this

         split_recs_into_sets  one set has the CR regions found with its flanking tRNA, an other has the rest of the mito oriented from Phe forward
         CR_analysis           look at the control region (or regions) and assess complexity

         assemble_mito          assembles the mito, will exclude the control region in some cases. uses megahit and an alignment method
            assemble_w_megahit  run the megahit program and pull out the highest scoring match
            assemble_w_mafft    run mafft and a consensus building script to build the mito assembly
            asm_compare         compare megahit and kalign assemblies

         assemble_CR       CR processed separately for hetroplasmy and tandem repeat analysis if earlier analysis shows a need
            trf_analysis   perform tandem repeat analysis on control region excerpts using the trf program

         complete          put a representative version of the control region assemblies together with the rest of the mito assembly

      -- Auxiliary
         rrna_cmsearch
         multi_mitfi
         add_goosehairpin_to_mitfi
         mitoCodeFromAccnOrTaxid
         make_mito_rec_cand_tsv  culls the blasted records based on query coveage percentage -- to exclude NUMTs etc.
'
elif [ $cmd == "init" ]; then
    echo -e '
    Usage: hifimitie init [<taxid number> | <taxonomy term> | - ] [HiFi file]...

    creates the hfmt working directory if necessary
    places fastq or fasta HiFi filename paths into files_to_search.fofn in the hfmt working directory

    taxid number is used to create a taxidlist for blast searches of mitogenomes (use 0 to search them all)
    using a term will show the taxonomy list info. when you find the correct number, enter it after Quitting the view

    use a dash - as the taxid argument to keep the current taxid and set the files with rest of the arguments
'
elif [ $cmd == "gettaxid" ]; then
    echo -e '
    Usage: hifimitie gettaxid [<taxid number | taxonomy term>]

    taxid number is used to create a taxidlist for blast searches of mitogenomes (use 0 to search them all)
    using a term will show the taxonomy list info. when you find the correct number, enter it after Quitting the view

    you also have an opportunity to look at the names of mitochondrial genomes associated with the given taxid
'
elif [ $cmd == "addfiles" ]; then
    echo -e '
    Usage: hifimitie addfiles [HiFi file]...

    creates the hfmt working directory if necessary
    places fastq or fasta HiFi filename paths into files_to_search.fofn in the hfmt working directory
'
else
   echo -e "(usage) command \"$cmd\" not found"
fi

: '
01_init.sh
02_blast_to_mito.sh
02.B_make_mito_rec_cand_tsv.sh
03_pull_fofn_cand_recs.sh
04_select_mito_features.sh
05_blast_features.sh
06_split_recs_into_2_sets.sh
07_assemble_w_megahit
08_align_call_consensus
09_process_control_region_recs
10_meld_mito_w_cr
11_annotate
'
