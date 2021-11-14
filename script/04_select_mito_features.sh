#!/bin/bash

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

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

### define the names

binfile=${wdir}/top_mito_species_match_coverage.bins
feature_output_file=${wdir}/top_match_feature_sequences.fasta

### define the functions

function select_feature_IDs {
   unset IDs
   [ ! -f $binfile ] && msglog "[select_mito_features.sh]: ${wdir}/top_mito_species_match_coverage.bins not found" && exit 2

#   IDs=$(head -n1 $binfile)
   IDs=$(awk '{sub("  *##.*","")}{print; exit}' $binfile)
   update_setting top_mito_matches "$IDs"
}

function retrieve_features {

   msglog_module "features being retrieved from $@"

   # for each feature we do like this, but use -rec instead of -pro
   # input from mito_analyze.py -q KU745736 -pro
   # TCAAGAAGGAAGGACAAACCCCCCACCATCAGCACCCAAAGCTGACATTCTCAATAAACTACTTCCTG	Pro 15356 15423 KU745736 Neotoma fuscipes mitochondrion, complete genome

   # output:
   # >Pro_KU745736 15356 15423 Neotoma fuscipes mitochondrion, complete genome
   # TCAAGAAGGAAGGACAAACCCCCCACCATCAGCACCCAAAGCTGACATTCTCAATAAACTACTTCCTG

   for query in $@; do
      msglog_module "retrieving mito features for $query"
      msglog_module "mito_analyze.py -rec -q $query -nh | sort -k3,3n |" bioawk_cas '{print ">" $2"_"$5, $3, $4, fldcat(6,NF); print $1}' "| fold -w 120"
      mito_analyze.py -rec -q "$query" -nh | sort -k3,3n | bioawk_cas '{print ">" $2"_"$5, $3,$4,fldcat(6,NF); print $1}' | fold -w 120
   done
}

###############################################################################
#                    retrieve the features from GenBank                       #
###############################################################################

if [ ! -s $feature_output_file ]; then

   select_feature_IDs
   retrieve_features $IDs > $feature_output_file
   [ -s $feature_output_file ] && msglog_module "$(basename $feature_output_file) file created"

else

   msg "$(basename $feature_output_file) file already created"

fi
