#!/bin/bash

source $(dirname $(realpath $0))/shared.sh # sets src_dir, wdir, msglog, is_number functions and usage among other things
[ -z "$1" ] && $usage gettaxid && exit 1

# gettaxid's ultimate job is to set the taxid to one of these:
# 1) 0 meaning don't restrict mitogenome search, search the complete mito database
# 2) a single species taxid to search just one mitogenome (not particularly recommended
# 3) a taxidlist file containing all the taxids to be included in the mito db search

# the second above really looks like the third one with just a single taxid entry in the taxidlist file

# that is what happens when a number is given as the argument, taxidlist file is made and info about it is reported

taxtok=$1

function view_mitogenome_names {
   tid=$1
   (echo List of $num_mitogenomes mitogenomes associated with taxid $tid. Press Q to exit. && (blastdbcmd -db mito -taxidlist <(make_taxidlist.sh $tid) | grep "^>")) | less -N
}

while true; do

   if is_number $taxtok; then
      num_ids=$(make_taxidlist.sh $taxtok |wc -l)
      lineage=$(grep -m1 "^${taxtok}\s" $fullnamelineage)
      num_mitogenomes=$(blastdbcmd -db mito -taxidlist <(make_taxidlist.sh $taxtok) | grep -c "^>")

      msg "\ntaxid $lineage\n"
      msg "taxid $taxtok has $num_mitogenomes mitogenomes matching one of the $num_ids taxids associated with it. These $num_mitogenomes will be used to limit the mitogenome search.\n"

      while true; do   # reask if V used to view names of the mitogenomes or an invalid option is used
         read -p "OK [Y|N|V]? " resp; firstCharUC_resp
         ( [ "$resp" == "Y" ] || [ "$resp" == "N" ] ) && break
         [ "$resp" == "V" ] && view_mitogenome_names "$taxtok"
      done

      [ "$resp" == "Y" ] && break
      echo ""

      # No, keep looking
      taxtok=""
   else
      [ ! -z $taxtok ] && make_taxidlist.sh "$taxtok" -C1
      read -p "Enter genus or species taxid to see full taxonomy taxids or another term to search more (Q to quit): " species
      resp=$species; firstCharUC_resp; [ $resp == "Q" ] && exit 0
      echo ""

      ! is_number $species && taxtok=$species && continue

      taxid_info.sh $species
      echo ""

      read -p "Taxid to identify mitogenome to search (or another term to search more, Q to quit): " taxtok
      resp=$taxtok; firstCharUC_resp; [ $resp == "Q" ] && exit 0
   fi

done
echo ""

# if we get here we have a taxtok that is a number
set_wdir
[ -z $wdir ] && msg "The hfmt_<num> directory into which to save taxidlist was not found" && exit 2
! is_number $taxtok && msg "Expected $taxtok to be a number" && exit 3

# write 2 files, one with number that got us our list and the other the list of taxid numbers we will use for limiting the mito search
topid_name=$(awk -F "\t" -v taxid=$taxtok '$1==taxid{print $3; exit}' $fullnamelineage)
echo $taxtok $topid_name > ${wdir}/toptaxid

make_taxidlist.sh $taxtok >${wdir}/taxidlist
