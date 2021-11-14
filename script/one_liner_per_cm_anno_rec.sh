#!/bin/bash

right_neighbors=cm_anno_right_neighbor.matrix
cm_anno_recs=mito_hifi_recs.cm_anno

cand=$(ls -td hfmt_[0-9]*/ 2>/dev/null | head -1)

function msg { # write a msg to stderr
   >&2 echo -e "$@"
}

[ "$1" = "sort" ] && do_sort="sort"

# find the directory
dir=""
[ -s $right_neighbors ] && [ -s $cm_anno_recs ] && dir="."
[ -z $dir ] && [ -s cm_results/$right_neighbors ] && [ -s cm_results/$cm_anno_recs ] && dir=cm_results
[ -z $dir ] && [ -d $cand ] && [ -s $cand/cm_results/$right_neighbors ] && [ -s $cand/cm_results/$cm_anno_recs ] && dir=$cand/cm_results

[ -z $dir ] && msg "Could not find $right_neighbors and/or $cm_anno_recs" && exit 1

# $7 for trna is the codon, $8 for trna is the model that best fits.
# sometimes these differ and we have been using the codon, better to use model
# but in either case we will add a tilde after the symbol if $7 and $8 conflict

# m64044_210818_004812/82445376/ccs_RC    1221    1294    55.22   2.305E-10       UGA     S2      Metazoa_S2.cm   -       .
# m64044_210818_004812/82445376/ccs_RC    1297    1365    51.83   3.86E-10        GUC     D       Metazoa_D.cm    +       3

# this will just cat them out unless $do_sort is defined and then we will sort based on first field and then remove that field from the output
function output_results {
   if [ -z $do_sort ]; then
      cat
   else
      sort -k1,1n | bioawk_cas -t '{ print fldcat(2,NF) }'
   fi
}

awk -v do_sort=$do_sort '
   BEGIN{rec_examp="m64049_191221_223534/20775426/ccs"; gsub(".", " ", rec_examp); spaces=rec_examp  }

   function set_spaces(nm) {
      spaces = ""
      for(s=1; s<=nm; s++)
      spaces = spaces " "
   }

   /^#/{next}
   NR==1 {
      $1=$1 # to force space as field sep on output
      sub("12S", "rrnS")
      sub("16S", "rrnL")
      pre=""; if (do_sort) pre="-1\t"
      printf("%s%s\t %s %s", pre, spaces, $0, $0)

      pos_str = sprintf("\n %s", $0)
      n = split(pos_str, syms, " ")
      for (i=1; i<=n; i++) {
         sym = syms[i]
         p = index(pos_str, sym)
         pos_ar[sym] = p
      }
      next
   }
   FNR==NR{next}

   $5 > .001 { next }

   { aa_sym=$7; cm_sym=$7 }

   $8 ~ "^Metazoa_" {  # reduce it to aa_sym format
      cm_sym = $8
      sub("Metazoa_", "", cm_sym)
      sub(".cm$" ,"", cm_sym)
   }

   { sym = aa_sym }  # this defines which one of the two we will use
   { addtl= (aa_sym==cm_sym) ? "" : "~" }

   $1!=lst {
      p = pos_ar[sym]-3; set_spaces(p)
      pre=""; if(do_sort) pre = p "\t"
      printf("\n%s%s\t%s %s%s", pre, $1, spaces, sym, addtl);
      lst=$1
      next
   }

   { printf(" %s%s", sym, addtl) }

   END{printf("\n")}
' ${dir/}/$right_neighbors ${dir}/$cm_anno_recs | output_results
