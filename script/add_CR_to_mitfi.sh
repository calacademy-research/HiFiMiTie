#!/bin/bash

# give settings files and mitfi file at the CR info

: '
num_CRs	1
CR1_flanks	P F
CR1_mean	828
CR1_type	Control Region
Primary_CR	CR1
CR1_stddev	1
CR1_recs_w_CR	1174

m64044_210611_022728/29623447/ccs       1762    1829    62.69   1.098E-12       GAA     F       Metazoa_F.cm    +
'

settings=$1
mitfi=$2
asmlen=$3

if [ -z $3 ]; then

   msg.sh "usage: add_CR_to_mitfi.sh <settings_file> <mitfi_file> <assembly_length>"
   [ -s $mitfi ] && cat $mitfi # pass the original on in the pipeline

else

 awk -v asmlen=$asmlen '
   BEGIN{OFS="\t"; dt="."}
   function write_CR(start, end) {
      print asm_name, start, end, dt,dt,dt, "cr", "Control Region", "+"

      if (gh_str!="") {
         print gh_str; gh_str = ""
      }

      CR_starts_at = 0
   }

   FNR==NR && /^CR[1-9]_flanks/ {
      prev_trna = substr($2, 1,1); succ_trna = substr($3, 1,1)   # BUG: expecting 1 letter which will not always be the case
      id = substr($1, 1,3)
      cr_id[prev_trna] = id
      cr_prev[id] = prev_trna
   }
   FNR==NR && /^CR[1-9]_type/ {
      id = substr(1,3, $1)
      cr_type[id] = $2
   }
   FNR==NR{ next }

   $7 in cr_id {  # this record comes right before a CR
      CR_starts_at = $3+1 # start one gt end of this tRNA
      asm_name = $1
      print
      next
   }

   $7=="gh" {if(gh_str==""){ gh_str=$0 } else{ gh_str=gh_str"\n"$0 };next}
   CR_starts_at > 0 && $7 != "gh" && $7 != "OH" && $2 > CR_starts_at {  # ends one pos before beginning of this feature
       write_CR(CR_starts_at, $2-1)
   }

   { print }

   END {
      if (CR_starts_at > 0)
         write_CR(CR_starts_at, asmlen)
 } ' $settings $mitfi | sort -k2,2n -k3,3nr

fi
