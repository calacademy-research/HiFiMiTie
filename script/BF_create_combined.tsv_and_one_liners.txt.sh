#!/bin/bash

# combine the two and output in our anno format
: '
#header start   stop    score   evalue  AC      AA      model   strand
m64049_230121_044339/36307315/ccs       15204   15640   233     4.49e-60        .       OH      blastn  +       .
'

# inputs
: '
Thr_NC_026102   m64044_201011_075919/459824/ccs 100.000 45      0       0       1       45      299     343     6.77e-17        82.4    N/A     QPct:100        45      100     100

m64044_201011_075919/590660/ccs_RC      NC_026102:OH-1-4347-4909        96.631  564     18      1       8709    9272    1       563     0.0     935     N/A     Demodex folliculorum    14165   4       4
'

# also create one_liners.txt

function add_separators {
   cawk -t '
      (NR==1) { print "# header start   stop    score   evalue  AC      AA      model   strand"; name = $1  }
      $1 != name { print "#"; name = $1 }
      { print }
   '
}

function flip_OH_query_source {
   cawk -t '
      $11 > 1.0e-46 {next}  # this is value mitos uses for OH discimination
      {sub(":.*", "", $2); $2 = "OH_" $2}
      {print $2, $1, fldcat(3,6), $9, $10, $7, $8, fldcat(11,NF)}
   ' OH_blast_to_cand_recs.tsv
}

function fmt_feature_and_OH_tsv {
   cawk -t '
      FILENUM==1 || FILENUM==2 {
         ele = $1; sub("_.*", "", ele)
         score = $12; e_val = $11
         if (e_val == "0.0") e_val = "0.0     "
         if ($9 < $10) {
            start = $9; end = $10; strand = "+"
         } else {
            start = $10; end = $9; strand = "-"
         }

         print $2 " ", start, end, score, e_val, ".", ele, "blastn", strand
      }
   ' mito_feature_match_to_cand_recs.tsv <(flip_OH_query_source) |
   sort -k1,1V -k2,2n -k4,4nr | add_separators
}

function one_liners_from_combined_tsv {
   awk '
      /^#/ { next }
      $1 != lst { if(lst != "") print ""; printf("%s   \t", $1); lst = $1 }

      { first = getsym($7); printf("%s ", first) }
      END { print "" }

      function getsym(sym) {
         if (sym in seen)
            return seen[sym]

         cmd = "AAsyms.sh " sym
         cmd | getline symvals; close(cmd)
         split(symvals, ar); seen[sym] = ar[1]
         return ar[1]
      }
   ' combined_pcg_rna_OH.anno
}


###############################################################################
#           create combined_pcg_rna_OH.anno then combined_pcg_rna_OH.tsv       #
###############################################################################

if [ ! -s combined_pcg_rna_OH.anno ]; then
   fmt_feature_and_OH_tsv > combined_pcg_rna_OH.anno
fi

if [[ -s combined_pcg_rna_OH.anno && ! -s one_liners.txt ]]; then
   one_liners_from_combined_tsv > one_liners.txt
fi
