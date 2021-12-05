#!/bin/bash

# 13811	CCCCCCCTCCCCCCC 	m64044_210611_022728/15534806/ccs_RC
# -	not found       	m64044_210818_004812/165676943/ccs
# 8666	CCCCCCCTCCCCCCC 	m64044_210818_004812/82839461/ccs_RC

# m64044_210818_004812/82839461/ccs_RC    8541    8608    62.64   1.053E-11       UGU     T       Metazoa_T.cm    +
# m64044_210818_004812/82839461/ccs_RC    10038   10106   65.29   3.176E-12       UGG     P       Metazoa_P.cm    -

# add any found goosehairpins into proper place as a comment in mitfi file

mitfi=$1
[ ! -s "$mitfi" ] && mitfi=mito_hifi_recs.mitfi
[ ! -s "$mitfi" ] && msg.sh "add_goosehairpin_to_mitfi: can not find $mitfi" && exit 1

fasta=$2
[ -z $fasta ] && fasta=mito_hifi_recs.fasta
[ ! -s $fasta ] && fasta=../mito_hifi_recs.fasta
[ ! -s $fasta ] && msg.sh "add_goosehairpin_to_mitfi: can not find $fasta" && exit 2

awk '

   function look_for_gh(recid, trna_start_pos) {
      for(g = 1; g <= num_gh[recid]; g++) {

         gh_seq       = gh[recid][g]
         gh_start_pos = gh_pos[recid][g]

         if (gh_start_pos > 0 && gh_start_pos < trna_start_pos) { # insert the goosehairpin info before this trna info
            pos = gh_start_pos; end_pos = pos + length(gh_seq)-1
            printf("%s\t%d\t%d\t.\t.\t%s\tgh\tgoose_hairpin\t+\n", recid, pos, end_pos, gh_seq)
            gh_pos[recid][g] = 0  # so we just add it once
         }
      }
   }

   FNR==NR && int($1)>0 {
      recid = $3; seq = $2; pos = $1
      num_gh[recid]++; n = num_gh[recid]

      gh[recid][n] = seq
      gh_pos[recid][n] = pos
   }
   FNR==NR{ next }

   /^#/ {
      look_for_gh(lst_recid, 10000000) # if goose hairpin at after last trna this will output it
      print
      next
   }

   $1 in gh {
      recid = $1; trna_start_pos = $2
      look_for_gh(recid, trna_start_pos)
   }

   {
      lst_recid = $1
      print
   }

   END { # output if needed for the last record mifi info
      look_for_gh(lst_recid, 10000000) # if goose hairpin at after last trna this will output it

 } ' <(goosehairpin.sh $fasta) $mitfi
