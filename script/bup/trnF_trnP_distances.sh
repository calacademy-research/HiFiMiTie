#!/bin/bash

# puts the Pro and Phe liens together on a single line and removes the movename from the record names
# using bioawk_cas for the convenience of fldcat()

source $(dirname $(realpath $0))/shared.sh  # sets src_dir, wdir, script_dir, msglog, msglog_module, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

function create_dist_file {
   pos_tsv=${wdir}/mito_feature_match_to_cand_recs.tsv
   [ ! -f $pos_tsv ] && msglog_module "$pos_tsv could not be found" && exit 2

   grep -e trnP -e trnF -e "^Phe" -e "^Pro" $pos_tsv | sed  "s|m[0-9]*_[0-9]*_[0-9]*/||" | \
   bioawk_cas -t '
        function lrgpos() { return ($9<$10) ? $10:$9}
        lstrd==$2{ print lst, fldcat(1,10), "dist: " lrgpos()-lstpos+1 }

        {
           lst = fldcat(1,10)
           lstrd = $2
           lstpos = lrgpos()
        }
   ' > ${wdir}/trnF_trnP_distances.tsv
}

create_dist_file
