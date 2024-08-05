#!/bin/bash

function sort_oneliners_under_template {

   function get_template { ${dir}/BF_right_neighbor_list.sh $starting_element | head -n 1; }

   template=$(get_template)
   awk -v template="$template" 'BEGIN{print "\t\t\t\t\t " template " " template}'

   cawk -v template="$template" '
      NR == 1 {
         spc_prefix_bld()
         match_template = template " " template  # string 2 together for comparing with each reads elements
         prepare_anno("dummy\t"match_template)
         match_template = anno  # removes OH gh and the like
         next  # first line is the template and we have already printed it above
      }
      {  # this is where we move the 2nd field of the line over with spaces
         set_msg($0)
         print length(spc_prefix[$2]), $1 "  \t" spc_prefix[$2], fldcat(2,NF) msg
      }

      function spc_prefix_bld(   nf, f) {
         nf = split(template, flds)
         for (f = 1; f <= nf; f++) {
            spc_prefix[$f] = sprintf("%-" spaces "s", "")
            spaces += length($f) + 1
         }
      }

      function set_msg(line) {
         prepare_anno(line)
         found_at = index(match_template, anno)
         if (found_at < 1)
            msg = "  # T template mismatch"
         else if (num_annotated <= too_few)
            msg = "  # N too few items"
         else if (too_many)
            msg = "  # W wraps around"
         else
            msg = ""
      }

      function prepare_anno(line) {
         too_few = 7
         too_many = 0 # counts number of first annotation elements in rest of annotation
         anno = line

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
         sub("^ *", "", anno)  # fixes problem if OH or gh or OL is the very first item
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

   ' <(get_template) one_liners.txt | sort -k1,1n | sed -e "s/^[0-9]*//" -e "s/^ //"
}

dir=$(dirname $(realpath $0))
starting_element=$1

[ ! -s right_neighbor_list ] && ${dir}/BF_right_neighbor_list.sh $starting_element > right_neighbor_list.txt

sort_oneliners_under_template
