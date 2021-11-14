#!/bin/bash

: '
>Phe_NC_007975 1 69 Cnemotriccus fuscatus mitochondrion, complete genome
>12S_NC_007975 70 1050 Cnemotriccus fuscatus mitochondrion, complete genome
>Val_NC_007975 1051 1120 Cnemotriccus fuscatus mitochondrion, complete genome
>16S_NC_007975 1121 2700 Cnemotriccus fuscatus mitochondrion, complete genome
>Leu2_NC_007975 2701 2774 Cnemotriccus fuscatus mitochondrion, complete genome
'

# depending on the organism, there are different canonical formats for the mito layout
# these all (or, most) seem to come the placement of the first trna in the mito sequence
# for example, vertebrates start with Phe (F), moths, ant, butterflies with Met (M)
# beetles look like they start with Ile (I).
# we have pulled down the features from mitos most closely resembling the current focus
# of analysis, so we will find the first trna in the feature set and presume that
# represents the canonical format. if we get know, then we default to F later on

feature_file=$1
[ -z $feature_file ] && feature_file=mito_feature_sequences.fasta

function get_feature_names {
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
   AAsyms.sh $(get_feature_names) | grep -v Unknown | head -1 | cut -c1
}

get_first_trna_from_features
