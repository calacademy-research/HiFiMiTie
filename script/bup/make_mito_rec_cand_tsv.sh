#!/bin/bash

src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msg "The hfmt_<num> working directory not found" && exit 2

# get tophit for each matching read, which is the best as ranked by blastn
#  and write those where the query coverage was >=60% to mito_rec_candidates.tsv
#  sorting the output by the query coverage

# for query coverage bins, 100 and 60s, 70s, 80s,  90s, show the hit distributions

cutoff=60
top_cutoff=80  # all qcovs scores 80 upto 100 pct of query coverage

tsv_recs="$@"
[ -z "$tsv_recs" ] && msg "need to enter the name(s) of the tsv file(s)" && exit 1

out_dir=${wdir}/hifi_mito_matches

[[ ! "$tsv_recs" == *"hifi_mito_matches/"* ]] && tsv_recs=${out_dir}/*.tsv

sorttab -k1,1 -k12,12nr $tsv_recs | tophit | awk '$NF>9' | sorttab -k17,17nr > ${wdir}/mito_rec_candidates.tsv

awk -v qpct_cutoff=$cutoff -v top_qpct_cutoff=$top_cutoff 'BEGIN{FS="\t"; OFS="\t"}
   function min(a,b) {return (a<=b) ? a : b}; function max(a,b) {return (a>=b) ? a : b}
   function show_bucket(b) {
      s = ((10*b)==100) ? "" : "s"
      printf("%d%s\t%4d", b*10,s, dist[b])

      PROCINFO["sorted_in"] = "@val_num_desc"
      for(sp in species[b])
         printf("\t%4d:%s", species[b][sp], sp)

      printf("\n")
   }
   function show_IDs_for_feat_list_choice() { # we want this output to be the first line in our bin file
      PROCINFO["sorted_in"] = "@val_num_desc"
      first_ones_hits = 0
      for(id in top_bucket) {
         if(first_ones_hits==0) {
            first_ones_hits = top_bucket[id]
            threshold_for_feat_retrieval = max(5, int(first_ones_hits/2))
         }

         if(top_bucket[id] >= threshold_for_feat_retrieval) {
            to_retrieve++
            printf("%s\t", id)
         }
      }
      if(to_retrieve==0)
         printf("## not enough hi quality matches found for retrival of features\n")
      else
         printf(" ## %d chosen for downloading of the mito feature set\n", to_retrieve)
   }

   {qcov = $NF}

   qcov >= qpct_cutoff {
      bucket = int(qcov/10)
      if( ! (bucket in dist)) num_buckets++
      dist[bucket]++
      sub(" *[Mm]itochondrion, complete genome", "", $14)
      sub("\\.[1-9][0-9]*","", $2)
      sp_info = $2" "$14
      species[bucket][sp_info]++
   }
   qcov >= top_qpct_cutoff {  # we will choose from which of the mitogenomes to pull feature lists based on reprsentation in the top bucket
      top_bucket[$2]++
   }
   END {
      show_IDs_for_feat_list_choice()  # this line has what we will read to use in the mito_analyze.py query to download the feature set(s)
      show_bucket(10)
      for(b = 9; b>0; b--) {
         if(b in dist) {
            show_bucket(b)
         }
      }
   }
' ${wdir}/mito_rec_candidates.tsv > ${wdir}/top_mito_species_match_coverage.bins

# top_mito_species_match_coverage.bins has the best fit accession number as first token in first line. use it to find the appropriate mito genetic code
codeFromAccn=$(get_script "mitoCodeFromAccnOrTaxid")
top_accn=$(awk '{print $1; exit}' ${wdir}/top_mito_species_match_coverage.bins)
code_rslt=$($codeFromAccn $top_accn)
mito_code=$(echo $code_rslt | awk '{print int($1)}')

[ ! -z "$code_rslt" ] && msglog_module "top accn info: $code_rslt"
[ -z "$mito_code" ] && mito_code=5
msglog_module "mito code chosen: " $mito_code

# now remember the code in the ${wdir_path}/settings file
update_setting "code" $mito_code
