#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things

set_wdir # set the working dir variable wdir
[ -z $wdir ] && msg "[select_mito_features.sh]: The hfmt_<num> working directory not found" && exit 2

# look in top_mito_species_match_coverage.bins to figure out which IDs to get the features from
: '
cat top_mito_species_match_coverage.bins
NC_033356	 ## 1 chosen for downloading of the mito feature set
100	  45	  45:NC_033356 Neotoma fuscipes
90s	1266	1265:NC_033356 Neotoma fuscipes	   1:NC_035594 Neotoma mexicana
80s	  13	  13:NC_033356 Neotoma fuscipes
70s	   7	   4:NC_033356 Neotoma fuscipes	   2:NC_035594 Neotoma mexicana	   1:NC_035598 Peromyscus pectoralis
60s	   9	   6:NC_035594 Neotoma mexicana	   2:NC_035598 Peromyscus pectoralis	   1:NC_033356 Neotoma fuscipes
'

binfile=${wdir}/top_mito_species_match_coverage.bins

function select_feature_IDs {
   unset IDs
   [ ! -f $binfile ] && msglog "[select_mito_features.sh]: ${wdir}/top_mito_species_match_coverage.bins not found" && exit 2

#   IDs=$(head -n1 $binfile)
   IDs=$(awk '{sub("  *##.*","")}{print; exit}' $binfile)
}

function retrieve_features {

   msg "features being retrieved from $@"

   # for each feature we do like this, but ue -rec instead of -pro
   # input from mito_analyze.py -q KU745736 -pro
   # TCAAGAAGGAAGGACAAACCCCCCACCATCAGCACCCAAAGCTGACATTCTCAATAAACTACTTCCTG	Pro 15356 15423 KU745736 Neotoma fuscipes mitochondrion, complete genome

   # output:
   # >Pro_KU745736 15356 15423 Neotoma fuscipes mitochondrion, complete genome
   # TCAAGAAGGAAGGACAAACCCCCCACCATCAGCACCCAAAGCTGACATTCTCAATAAACTACTTCCTG

   for query in $@; do
      msg "retrieving mito features for $query"
      mito_analyze.py -rec -q "$query" -nh | sort -k3,3n | bioawk_cas '{print ">" $2"_"$5, $3,$4,fldcat(6,NF); print $1}' | fold -w 120
   done
}

if [ ! -s "${wdir}/mito_feature_sequences.fasta" ]; then
   select_feature_IDs
   retrieve_features $IDs > ${wdir}/mito_feature_sequences.fasta
else
   msg "mito_feature_sequences.fasta file already created"
fi
