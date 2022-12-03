#!/bin/bash

: '
>Phe_NC_007975 1 69 Cnemotriccus fuscatus mitochondrion, complete genome
>12S_NC_007975 70 1050 Cnemotriccus fuscatus mitochondrion, complete genome
>Val_NC_007975 1051 1120 Cnemotriccus fuscatus mitochondrion, complete genome
>16S_NC_007975 1121 2700 Cnemotriccus fuscatus mitochondrion, complete genome
>Leu2_NC_007975 2701 2774 Cnemotriccus fuscatus mitochondrion, complete genome
'

# should include, ie source, this file in a file after the shared.sh functions have been loaded
[ -z "$shared_loaded" ] && echo -e "Could not load the shared functions needed from the shared.sh file." >&2

# depending on the organism, there are different canonical formats for the mito layout
# these all (or, most) seem to come the placement of the first trna in the mito sequence
# for example, vertebrates start with Phe (F), moths, ant, butterflies with Met (M)
# beetles look like they start with Ile (I).
# we have pulled down the features from mitos most closely resembling the current focus
# of analysis, so we will find the first trna in the feature set and presume that
# represents the canonical format. if we get know, then we default to F later on

function get_feature_names {
   feature_file=$1
   bawk '
      {
         sub("_.*", "", $name)
         printf(" %s", $name)
      }
      END {
         printf("\n")
   }' $feature_file
}

function get_first_trna_from_features {
   AAsyms.sh $(get_feature_names $1) | grep -v Unknown | head -1 | cut -c1
}

function get_first_trna_from_feat_or_taxname {
   feature_file=$1
   [ -z $feature_file ] && feature_file=$wdir/top_match_feature_sequences.fasta
   feat_start_trna=$(get_first_trna_from_features $feature_file)
   update_setting_if_changed "top_mito_match_first_trna" "$feat_start_trna"

   taxname=$(get_setting "taxname")
   if [ ! -z "$taxname" ]; then
      tax_trna_start_distribution=$(mito_trna_start_dist.sh $taxname | head -n 1)
      update_setting_if_changed "taxname_trna_starts" "$tax_trna_start_distribution"
      taxname_start_trna=$(echo $tax_trna_start_distribution | awk '{print $2; exit}')
      update_setting_if_changed "taxname_first_trna" "$taxname_start_trna"
   fi

   # now handle cases of either feat_start_trna or taxname_start_trna missing or existing and different
   unset trna_start_chosen
   if [ -z $taxname_start_trna ]; then
      if [ -z $feat_start_trna ]; then
         trna_start_chosen="F"
         set_msg="no feature or taxname $taxname starting trna found, setting default trna start of F, i.e. Phe"
      else
         trna_start_chosen=$feat_start_trna
         set_msg="found starting trna for the feature set but none for the taxname $taxname, setting first trna to $feat_start_trna"
      fi
   elif [ -z $feat_start_trna ]; then
      trna_start_chosen=$taxname_start_trna
      set_msg="found starting trna for the taxname $taxname but none for the feature set, setting first trna to $taxname_start_trna"
   elif [ $taxname_start_trna = $feat_start_trna ]; then
      trna_start_chosen=$feat_start_trna
      set_msg="both feature set and taxname $taxname have starting trna of $feat_start_trna, setting first trna to $feat_start_trna"
   else
      trna_start_chosen=$feat_start_trna
      set_msg="feature set has starting trna of $feat_start_trna while taxname $taxname has $taxname_start_trna, setting first trna to feature set value $feat_start_trna"
   fi
}

# calling this should make sure the first_trna variable is set and recorded in the settings file (if already there it is not changed)
function get_first_trna_setting {
   mm_fasta=$wdir_path/top_match_feature_sequences.fasta
   top_mm_first_trna=$(bawk '{print $name; exit}' $mm_fasta)
   [ ! -z $top_mm_first_trna ] && update_setting_if_changed "top_mito_match_first_recname" "$top_mm_first_trna"

   # this will set the one letter code for the first_trna in trna_start_chosen and a message about why it was chosen in set_msg
   get_first_trna_from_feat_or_taxname $mm_fasta

   current_first_trna=$(get_setting "first_trna")
   if [ -z $current_first_trna ]; then
      msglog_module $set_msg
      first_trna=$trna_start_chosen
   else
      first_trna=$current_first_trna
   fi

   # if there is an existing setting for first_trna
   [ ! -z $first_trna ] && update_setting_if_changed "first_trna" $first_trna
}

# call the following from the script you source this file from,
# you can pass in the feature file location or it will use the typical one if you do not
: get_first_trna_setting
