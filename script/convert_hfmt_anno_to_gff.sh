#!/bin/bash

# convert our anno style into gff3

: '
Chr1_mm	cmsearch	gene	3172239	3172348	.	+	.	ID=gene-Gm26206;Dbxref=GeneID:115487594,MGI:MGI:5455983;Name=Gm26206;gbkey=Gene;gene=Gm26206;gene_biotype=snRNA

msa_consensus_mitochondrion     14870   16299   .       .       .       cr      Control Region  +       1430    1
msa_consensus_mitochondrion     14926   14940   .       .       CCCCCCCTCCCCCCC gh      goose_hairpin   +       15      56
msa_consensus_mitochondrion     16974   17151   .       .       .       cr      Control Region  +       178     1
msa_consensus_mitochondrion	16981	18376	3489	.	.	repeat	Repeat Region	+	34.0 copies 60 bp consensus: CCCCCCGTTCGGGCTTTGCTTAAGTCCATGCTAATATATTTCCTTTTTTTTTCGTCCGCA

# gff conversion
msa_consensus_mitochondrion	mitfi	tRNA	1	68	1.935E-12	+	.	symbol=F;anticodon=GAA;model=Metazoa_F.cm;distance_from_last=1
msa_consensus_mitochondrion	cmsearch	rRNA	68	1043	2.2e-183	+	.	symbol=rrnS;type=complete;model=12S_rRNA.cm;distance_from_last=0
msa_consensus_mitochondrion	blastn	gene	2793    3760	0.0	.	.	symbol=ND1;model=ND1_NC_007975;distance_from_last=25
msa_consensus_mitochondrion	hifimitie	mitochondrial_control_region	14870   16299   .	+	.	symbol=cr;type=Control Region;distance_from_last=1
msa_consensus_mitochondrion     hifimitie       stem_loop	14926   14940   .	+	.	symbol=gh;type=goose hairpin
msa_consensus_mitochondrion     hifimitie       mitochondrial_control_region	16974   17151   .	+	.	symbol=cr;type=remnant;distance_from_last=1

#header	start	stop	score	evalue	AC	AA	model	strand	length	dist_from_last
msa_consensus_mitochondrion	1	68	62.69	1.935E-12	GAA	F	Metazoa_F.cm	+	68	0
msa_consensus_mitochondrion	68	1043	708.3	2.2e-183	complete	rrnS	12S_rRNA.cm	+	976	0
msa_consensus_mitochondrion	1043	1113	62.49	7.348E-12	UAC	V	Metazoa_V.cm	+	71	0
msa_consensus_mitochondrion	1119	2695	1099.5	1e-281	complete	rrnL	16S_rRNA.cm	+	1577	6
msa_consensus_mitochondrion     2793    3760    1066    0.0     .       ND1     ND1_NC_007975   .       968     25
'

anno=$1

awk '
   BEGIN {
      FS="\t"; OFS="\t"
      print "##gff-version 3"
   }
   /^#/{next}

   {
      line++; start=$2
      dist_from_last = (line==1) ? 1 : start - end # this ones start last ones end
      name=$1; end=$3; score=$5; model=$8; strand=$9; phase="."
   }

   model ~ "^Metazoa" {  # mitfi trna
      id = sprintf("symbol=%s;anticodon=%s;model=%s;distance_from_last=%d", $7, $6, model, dist_from_last)
      print name, "mitfi", "tRNA", start, end, score, strand, phase, id
      next
   }
   model ~ "_rRNA.cm" {  # 12S or 16S rRNA
      id = sprintf("symbol=%s;type=%s;model=%s;distance_from_last=%d", $7, $6, model, dist_from_last)
      print name, "cmsearch", "rRNA", start, end, score, strand, phase, id
      next
   }
   model ~ "Control Region" {
      type = "Control Region" # maybe do more with this later, but depends info being in the anno file
      id = sprintf("symbol=cr;type=%s;distance_from_last=%d", type, dist_from_last)
      print name, "hifimitie", "mitochondrial_control_region", start, end, score, strand, phase, id
      cr_start = start; cr_end = end
      next
   }
   model ~ "goose_hairpin" {
      dist_from_last = start - cr_start # distance from control region start
      id = sprintf("symbol=gh;type=goose hairpin;sequence=%s;dist_from_CR_start=%d", $6, dist_from_last)
      print name, "hifimitie", "stem_loop", start, end, score, strand, phase, id

      end = cr_end  # use control region end to calculate distance from the next one
      next
   }
   model ~ "OL.cm" {  # OL
      id = sprintf("symbol=%s;type=%s;model=%s;distance_from_last=%d", $7, $6, model, dist_from_last)
      print name, "cmsearch", "rep_origin", start, end, score, strand, phase, id
      next
   }
   model ~ "OH.cm" {  # OH
      id = sprintf("symbol=%s;type=%s;model=%s;distance_from_last=%d", $7, $6, model, dist_from_last)
      print name, "blastn", "rep_origin", start, end, score, strand, phase, id
      next
   }
   model ~ "Repeat Region" {  # tandem repeat region
      split($NF, rptinf, " ")
      id = sprintf("symbol=repeat;type=Repeat Region;consensus=%s;bp=%d;copies=%s", rptinf[6], rptinf[3], rptinf[1])
      print name, "trf", "repeat_region", start, end, score, strand, phase, id
      next
   }
   {  # should be a gene
      id = sprintf("symbol=%s;model=%s;distance_from_last=%d", $7, model, dist_from_last)
      print name, "blastn", "gene", start, end, score, strand, phase, id
   }
' $anno
