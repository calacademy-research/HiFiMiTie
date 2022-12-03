#!/bin/bash

# separate file so we can show extensive usage / help info if we want to, without getting in way of main script

cmd=$1
hfmt_version=0.05
hfmt_version_date=01-Dec-2022
hfmt_title="HiFiMiTie version $hfmt_version -- Find & Analyze Metazoan Mitochondria from HiFi reads [$hfmt_version_date]"

# set color vars that we use or might use here
Blue="\033[0;34m"; BBlue="\033[1;34m"; UBlue="\033[4;34m"; NC="\033[0m" # No Color (reset)

if [ -z $cmd ]; then
    echo -e "
    $hfmt_title

    Usage: hifimitie <command> [options]

    Commands:
      -- (1) Setup

         ${BBlue}init${NC}         create working directory and set taxid and HiFi file(s)

           gettaxid      set to search all mitogenomes or only those associated with a certain taxid (optionally view taxonomy ids based on term)
           addfiles      add some or all of the PacBio HiFi file names with which to begin the search
           clearfiles    remove the entries already in the fofn filename list (just removes the names, files themselves are untouched)

      -- (2)...(12) Execution

         ${BBlue}run${NC}          run the default pipeline of the commands below, specifying number of threads to use with -t <num_threads>, defaults to 32

         -- Pipeline

         (2) blast_to_mito         blast the mito assemblies associated with the taxids you have chosen against the PacBio HiFi files (uses local copy of GenBank mito db)

         (3) pull_fofn_cand_recs   create file with candidate mito recs based on the blast tsv file pulled from the HiFi files

         (4) select_mito_features  for the most closely matching mito assemblies, pull tRNA, rRNA, gene, CR features from GenBank into records
         (5) blast_features        blast the features against the candidate records using -task blastn instead of default megablast

         (6) rna_search            rrna_cmsearch of SSU, LSU rrna. multi_mitfi for trna, goose hairpin info added into this

         (7) CR_analysis           look at the control region (or regions) and assess complexity
         (8) split_recs_into_sets  one or more sets has the CR region(s) found with flanking tRNA, another has the rest of the mito oriented from starting tRNA forward (usually Phe)

         (9) assemble_mito         assembles the mitogenome, will exclude the control region in some cases. uses megahit and an alignment method
               assemble_w_megahit    run the megahit program and pull out the highest scoring match
               assemble_by_aligner   run mafft or another aligner and a consensus building script to build the mito assembly

         (10) compare_assemblies   compare megahit and msa consensus created assemblies

         (11) assemble_CR          CR processed separately for hetroplasmy and tandem repeat analysis if earlier analysis shows a need
                trf_analysis         perform tandem repeat analysis on control region excerpts using the trf program

         (12) complete             put a representative version of the control region assemblies together with the rest of the mito assembly

      -- Auxiliary
         ${BBlue}settings${NC}      ${BBlue}check${NC}       mitoCodeFromAccnOrTaxid       mitodb_update.sh
"
elif [ $cmd == "init" ]; then
    echo -e '
    Usage: hifimitie init [<taxid number> | <taxonomy term> | - ] [HiFi file]...

    creates the hfmt working directory if necessary
    places fastq or fasta HiFi filename paths into files_to_search.fofn in the hfmt working directory

    taxid number is used to create a taxidlist for blast searches of mitogenomes (can use 0 to search them all, but better to be more specific)
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
    Usage: hifimitie addfiles [HiFi file(s)]...

    creates the hfmt working directory if necessary
    places fastq or fasta HiFi filename paths into files_to_search.fofn in the hfmt working directory
'
elif [ $cmd == "settings" ]; then
    echo -e "
    Usage: hifimitie settings [ all | <keyword> | <keyword> <value> ]

    view or modify/create settings created at various steps of the pipeline

    the most likely use for creating a setting before running the pipeline
    is if the trna that will be chosen to start the beginning of the assembly
    is not the standard order for your species.
    this setting is ${BBlue}first_trna${NC} and takes a one letter AA code.
"
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
07_split_recs_into_2_sets.sh
08_assemble_w_megahit
09_align_call_consensus
09_process_control_region_recs
10_meld_mito_w_cr
11_annotate
'
