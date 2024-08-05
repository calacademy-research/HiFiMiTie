#!/bin/bash

# 10Aug2022 -- add categorization of annotation for suspect one with comment at end of line
# 06Sep2022 -- remove gh, OL and OH items for template and record to compare with. I.e., ignore those for template order comparison
# 09May2023 -- change too few items to be 3 from 7

right_neighbors=cm_anno_right_neighbor.matrix
cm_anno_recs=mito_hifi_recs.cm_anno

TOO_FEW_DEFAULT=3  # was 7

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

function categorize_cm_annos {
   # compare each read's annotation to the template at the top of the file
   # also remove tildes ~, these are put in when called tRNA does not use the CM model we expect (L1 using l2.cm for example)

   # annotation lines that fail a test have a comment appended that describes the reason
   # to remove these from the list do grep -v "#"
   anno_file=$1

   awk  -v too_few=$TOO_FEW_DEFAULT ' # too_few refers to number of elements annotated
      BEGIN {OFS = "\t"}
      function prepare_anno() {
         too_many = 0 # counts number of first annotation elements in rest of annotation
         anno = $0

         # remove read name and spaces around annotation portion of line, also remove any tildes
         sub("^[^ \t]*[ \t]*", "", anno)
         sub("[ \t]*$", "", anno)
         gsub("~", "", anno)
         len1 = length(anno)

         # remove gh, OL, OH from consideration
         gsub("gh", "", anno)
         gsub("OL", "", anno)
         gsub("OH", "", anno)
         gsub("  *", " ", anno)
         sub("^ *", "", anno)  # 02Apr2023 fixes problem if OH or gh or OL is the very first item
         len2 = length(anno)

         num_annotated =  split(anno, ar, " ")
         if (num_annotated < 1)
            return

         # see if first item is repeated more than once, that means read wraps around
         itm = ar[1]
         for (i = 2; i <= num_annotated; i++) {
            if(ar[i] == itm)
               too_many++
         }
      }

      # first line is the template
      NR==1 {  # template line
         prepare_anno()
         template = anno
         print $0
         next
      }

      # compare lines after the first with template and perform a few other checks
      {
         prepare_anno()
         found_at = index(template, anno)
         if (found_at < 1)
            print $0, "# T template mismatch"
         else if (num_annotated <= too_few)
            print $0, "# N too few items"
         else if (too_many)
            print $0, "# W wraps around"
         else
            print $0
      }
   ' $anno_file
}

# this will just cat them out unless $do_sort is defined and then we will sort based on first field and then remove that field from the output
function output_results {
   if [ -z $do_sort ]; then
      cat
   else # we have put a number at line start to tell us how to sort, sort then remove that number when outputting result, then put categorization comment at end for suspect annos
      sort -k1,1n | bioawk_cas -t '{ print fldcat(2, NF) }' | categorize_cm_annos
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
      p = (sym in pos_ar) ? pos_ar[sym]-3 : length(pos_str)
      set_spaces(p)
      pre=""; if(do_sort) pre = p "\t"
      printf("\n%s%s\t%s %s%s", pre, $1, spaces, sym, addtl);
      lst=$1
      next
   }

   { printf(" %s%s", sym, addtl) }

   END{printf("\n")}
' ${dir/}/$right_neighbors ${dir}/$cm_anno_recs | output_results
