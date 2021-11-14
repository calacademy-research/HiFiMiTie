# HiFiMiTie
HiFiMiTie -- Find &amp; Analyze Metazoan Mitochondria from HiFi reads

This is a work in progress. And is currently made to run on CAS ccg setups. Several programs and scripts have not yet been uploaded.

HiFiMiTie takes advantage of the fact that the PacBio HiFi read cell output can contain a hundred to a thousand mitochondrial reads each representing most to all to sometimes several copies of the subject's mitochondria.

It relies on a local copy of the NCBI mito database being downloaded and several supporting files being constructed from that. Scripts for this will added in the supporting directory soon. The primary tools used are blastn, mitfi, cmsearch, megahit, mafft and a consensus script.

```
    HiFiMiTie version 0.01 -- Find & Analyze Metazoan Mitochondria from HiFi reads

    Usage: hifimitie <command> [options]

    Commands:
      -- (1) Setup

         init         create working directory and set taxid and HiFi file(s)

           gettaxid      set to search all mitogenomes or only those associated with a certain taxid (optionally view taxonomy ids based on term)
           addfiles      add some or all of the PacBio HiFi file names with which to begin the search
           clearfiles    remove the entries already in the fofn filename list (just removes the names, files themselves are untouched)

      -- (2)...(12) Execution

         run          run the default pipeline of the commands below, specifying number of threads to use with -t <num_threads>, defaults to 32

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
         settings      check       mitoCodeFromAccnOrTaxid       mitodb_update.sh
```
