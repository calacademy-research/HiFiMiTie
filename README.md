# HiFiMiTie
HiFiMiTie -- Find &amp; Analyze Metazoan Mitochondria from HiFi reads

**This is a work in progress.** And is currently made to run on CAS ccg setups. Several programs and scripts have not yet been uploaded. But the broad brush of the pipeline is represented in the scripts in the scripts directory. These are primarily bash and awk scripts calling out to other programs.

HiFiMiTie takes advantage of the fact that the PacBio HiFi read cell output can contain a hundred to a thousand mitochondrial reads each representing most to all to sometimes several copies of the subject's mitochondria.

It relies on a local copy of the NCBI mito database being downloaded and several supporting files being constructed from that. Scripts for this will be added in the support_files directory soon. The primary tools used are blastn, mitfi, cmsearch, megahit, trf, mafft and a consensus script.

Each read is run through rrna and trna annotation to determine order. Location of one or more Control Regions (CRs) is also discerned, along with CR tandem repeat analysis, and heteroplasmy analysis. HiFiMiTie was built around scripts that found non-standard tRNA (WACNY instead of WANCY) and substantial tandem repeat heteroplasmy in the CR of a fish mitochondrial genome.

You can see the basic pipeline in the help screen below. However, to execute the pipeline it is only necessary to use hifimitie (or its abbreviation hfmt) with init to set the taxonomy numeric identifier and the PacBio HiFi file(s) -- fastq or fasta accepted; then use run indicating the number of threads to use.

Example forthcoming.

```
    HiFiMiTie version 0.04 -- Find & Analyze Metazoan Mitochondria from HiFi reads

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

### Installation

Not now, but eventually this should require cloning this repo to say ~/bin/HiFiMiTie then creating a softlink `ln -s ~/bin/HiFiMiTie/hifimitie.sh ~/bin/hifimitie` and another one `ln -s ~/bin/hifimitie ~/bin/hfmt` and downloading the dependencies. The command ``hifimitie check`` checks for the existence of the major dependencies.

### Introduction

HiFiMiTie works on low-error long PacBio HiFi reads to identify and annotate metazoan mitochondrial reads and create a consensus mitochondrion from these reads.
This consensus is compared with a result from the program Megahit and based on the annotations of each the best one is chosen as the representative and an annotation
of it is created. Megahit creates a kmer based result very quickly and the sequence consensus based version takes quite a bit longer
-- mostly due to the time it takes to do covariance matrix analysis of the tRNA and rRNA elements of each of the subject mitochondrial reads.
This analysis and the blasts of the local mito database downloaded from GenBank is where execution threads are most useful if you can allow the program to use them.

It might not seem worthwhile for all the extra effort for the additional read analysis and consensus mitochondrion. However, although there is typically a
great deal of agreement between the Megahit result and the sequence consensus result, in almost every instance seen so far the Megahit version
incorrectly assembles some portion or portions of the mito genome, at a minimum the Control Region.

### Pipeline gloss

#### Init

To begin HiFiMiTie (aka hfmt) a taxonomic numeric is provided and one or more fasta or fastq files containing the PacBio HiFi reads (gzipped or not)
are provided in one or more `hfmt init` commands. The program guides the user through the choice of the taxonomic number if an alphabetic string is used.

Also the first time the hfmt init command is executed it creates a directory with the format hfmt_\<date\>, e.g. hfmt_110322.
Subsequent runs of the hfmt command in the directory will create files in this hfmt_\<date\> directory.

The taxonomic identifier has two purposes. One is to limit the blast search of the HiFi reads to likely relevant mitochondrial records in the local mito database.
Although the information about the number of relevant mitochondria is of interest, the entire mito db could easily be searched in about the same time as the subset
(at least at present when there are so few total mitochrondrial genomes in the GenBank mito database).

The other purpose is to determine the first tRNA from which to start the mitochondrial sequence.
It is also assumed this tRNA will often be the ending flank of a Control Region (though not true for Cnidaria for example). Though a second Control Region can be
identified in other positions of the sequence. The mito database has been processed so the distribution of first tRNAs for each taxonomic level can be interrogated.
So the taxonomic level is used as one way to determine tRNA to start the sequence. Another is the first tRNA of the most closely related hit from the mito db search.
And it is also possible using `hfmt settings first_trna <trna sym>` for the user to set this explicitly.

The taxonomic level is usually what is used unless the user manually sets this. One reason for this is that some groups rotate the mito sequences to start at F in all instances.
Though this is true for `Vertebrata: F 5811, L1 99, V 82, S2 45, T 38, I 29, P 10, L2 9, E 4, G 4, W 4, C 2, H 2` it is not so e.g. for
`Lepidoptera: M 523, L2 14, I 11, G 4, V 4, W 4, F 1, L1 1, Q 1` though some groups are adding lepidoptera mitogenomes that start with F to GenBank.

#### Mito DB blast search and Candidate record pull

The HiFi read file(s) set in the init command(s) are used as the query into the mito database with the taxidlist argument used
to limit the search to those matching the taxonomic numeric specified. The blast tsv output uses the qcovus field so that 
the total percentage of a query read that matches a mitochondrion is recorded.

The blast tsv list is read and then the HiFI reads files are read.
HiFi reads in the blast tsv are kept as a candidate set for those with a 50% or greater qcovus match.

#### Create set of Mito records with bootstrapping

The candidate set of HiFi reads is is used as input to Megahit to create a preliminary assembly and 


### Notes

First versions relied upon close mitochondrial references from other taxa. Though this is often the case it is not always a good idea to rely solely upon this.
So the first pass of blasting to the mitodb is going to be supplemented with an additional pass using these matched reads
if less than 500 reads are found in the blast output of the mito db.
The reads matched are used for a bootstrap version of the mitochondrion from megahit. This record is blasted against mito db and revcomped appropriately.
Then the single bootstrap mitochondrion record is made into a blast db and the input reads are blasted against it. The records found here are added to those found in 
the first pass -- though in all likelihood they are a superset of those.

Also, we found that high Control Region heteroplasmy could cause records to be missed in that first pass too. The bootstrap second pass will likely find more and perhaps all.
However, we are still doing yet another pass through the input reads. This time in the split records step where we search the CR flanks in the input reads and
where the begin flank is followed by the end flank in the read we extract those and the putative CR sequence in between. In our fish genome with high CR heteroplasmy
the first pass got only 22 CR sequences but this second pass using the flanks of those 22 yielded 133 sequences, 111 additional sequences.
