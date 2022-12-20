#!/bin/bash

fasta=$1
mitfi=$2  # could be mitfi or cmanno, but is in the mitfi format style

function msg {
   echo -e "$@" >&2
}
function usage {
   msg "\n    usage: mito_pcg_anno.sh <fasta> <mitfi or cm_anno file>\n"
   exit
}

[ -z $fasta ] && usage

db_prefix=/ccg/bin/HiFiMiTie_dir/PCG_AA

pcgs=( nad1 nad2 nad3 nad4 nad4l nad5 nad6 atp6 atp8 cox1 cox2 cox3 cob )
names=( ND1 ND2 ND3 ND4 ND4L ND5 ND6 ATP6 ATP8 COX1 COX2 COX3 CYTB )

: '
# what each line of output looks like, normal outfmt 6 but with PCG name prefixed
ND1 nd1.fas	NODE_1_length_16359_cov_985.237610_RC	NC_007936:nad1-1-2732-3686	82.079	279	50	0	2865	3701	23	301	2.90e-122	387	N/A	NC_007936:nad1-1-2732-3686 Cricetulus griseus	16359
ND2 nd2.fas	NODE_1_length_16359_cov_985.237610_RC	NC_013276:nad2-1-3901-4933	71.859	199	56	0	3963	4559	1	199	5.04e-67	229	N/A	NC_013276:nad2-1-3901-4933 Mesocricetus auratu16359

# mitfi format example
NODE_1_length_16359_cov_985.237610_RC	1	71	60.03	6.417E-12	GAA	F	Metazoa_F.cm	+
'

function best_pcg_hits {
   ix=0  # bash 0 indexes arrays
   for gene in "${pcgs[@]}"; do
      gene_name=${names[$ix]}
      ((ix++))

      db=${db_prefix}/${gene}.fas

      top_hit=$(blastx -db $db -query $fasta -outfmt "6 std staxid stitle qlen" | head -1)
      if [ ! -z "$top_hit" ]; then
         printf "%s\t%s\t%s\n" $gene_name ${gene}.fas "$top_hit"
      else
         msg $gene_name not found.
      fi
   done

   # add the OH element into the mix
   top_hit=$(OH_hits | head -1)
   if [ ! -z "$top_hit" ]; then
      printf "%s\t%s\t%s\n" "OH" "OH.fas" "$top_hit"
   else
      msg OH not found.
   fi
}

function OH_hits {
   local qry=$fasta
   local OH_db=/ccg/bin/HiFiMiTie_dir/data/OrgRepl/OH.fas

   blastn -db $OH_db -query $qry -outfmt  "6 std staxid stitle qlen" -subject_besthit -evalue 9e-70
}

function format_pcg_hits {
   if [ -z $mitfi ]; then  # for debugging if no mitfi file on cmd line
      cat | sort -k9,9n
   else # sort and put into mitfi style format fields
      sort -k9,9n | \
      cawk -t ' NF > 10 {
            name = $3; strand = "+"
            strt = $9; end = $10
            if(end<strt){t=strt; strt=end; end=t; strand="-"}
            evalue = $13; score = $14; pcg = $1; db = $2

            print name, strt, end, score, evalue, ".",  pcg, db, strand
         }
      '
   fi
}

function add_to_mitfi {
   local pcg_formatted_hits=$1
   local mitfi=$2

   if [ -z $mitfi ]; then
      echo mitfi file not found
      cat
   else
      cawk -t '  # basically sorting on start field, field 2
         function next_pcg() { pcg_ix++; strt = ar_strts[pcg_ix] }

         FNR==NR {
            ar_strts[ ++pcgs ] = $2
            ar_row[ pcgs ] = $0
            next
         }

         FNR==1 { next_pcg() }
         /^#/ { print; next }

         $2 > strt {
            for (; $2 > strt && pcg_ix <= pcgs; next_pcg() )
               print ar_row[ pcg_ix ]
         }

         { print } # always print a line from the mitfi file

         END {
            for (; pcg_ix <= pcgs; next_pcg() )
               print ar_row[ pcg_ix ]
      }' $pcg_formatted_hits $mitfi
   fi
}


#################################################################
#                        do the work                            #
#################################################################

best_pcg_hits $fasta | format_pcg_hits | add_to_mitfi - $mitfi
