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
set_mitodb  # defined in shared.sh, sets mitodb_path
confirmation_timeout=60  #  presume Y if we wait this long

# 22Jul2022 new way to get code using taxid and nodes.dmp file where mito code is in the 9th field of the record where taxid is field 1
function set_mito_code_from_taxid {
   local taxid=$1
   local mito_code_script=$(get_script "mito_genetic_code.sh") # script in script dir searcs nodes.dmp for taxid and returns mito code number

   local mito_code=$($mito_code_script $taxid $taxnodes)
   msglog_module "mito code from taxid $taxid: $mito_code"

   # now remember the code in the ${wdir_path}/settings file
   update_setting_if_changed "code" $mito_code
}

function view_mitogenome_names {
   tid=$1
   [ ! -s ${wdir}/taxidlist ] && make_taxidlist.sh $tid >${wdir}/taxidlist  # 05Dec2023 use existing file or create it
   (echo List of $num_mitogenomes mitogenomes associated with taxid $tid. Press Q to exit. && (blastdbcmd -db $mitodb_path -taxidlist ${wdir}/taxidlist | grep "^>")) | less -N
}

while true; do

   if is_number $taxtok; then
      msg "creating taxidlist for $taxtok..."
      make_taxidlist.sh $taxtok >${wdir}/taxidlist   # 05Dec2023 create the file so we do not have call make_taxidlist.sh multiply
      num_ids=$(wc -l ${wdir}/taxidlist)
      lineage=$(grep -m1 "^${taxtok}\s" $fullnamelineage)
      num_mitogenomes=$(blastdbcmd -db $mitodb_path -taxidlist ${wdir}/taxidlist | grep -c "^>")

      msg "\ntaxid $lineage\n"
      msg "taxid $taxtok has $num_mitogenomes mitogenomes in $mitodb_path matching one of the $num_ids taxids associated with it. These $num_mitogenomes will be used to limit the mitogenome search.\n"

      while true; do   # reask if V used to view names of the mitogenomes or an invalid option is used
         read -t $confirmation_timeout -p "OK [Y|N|V]? " resp; ret_code=$?; firstCharUC_resp
         [ $ret_code -gt 128 ] && resp="Y" && echo $resp  # it is Y if we time out
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
      resp=$species; firstCharUC_resp; [ $resp == "Q" ] && exit 4
      echo ""

      ! is_number $species && taxtok=$species && continue

      taxid_info.sh $species
      echo ""

      read -p "Taxid to identify mitogenome to search (or another term to search more, Q to quit): " taxtok
      resp=$taxtok; firstCharUC_resp; [ $resp == "Q" ] && exit 4
   fi

done
echo ""

# if we get here we have a taxtok that is a number
set_wdir
[ -z $wdir ] && msg "The hfmt_<num> directory into which to save taxidlist was not found" && exit 2
! is_number $taxtok && msg "Expected $taxtok to be a number" && exit 3

# write 2 files, one with number that got us our list and the other the list of taxid numbers we will use for limiting the mito search
topid_name=$(awk -F "\t" -v taxid=$taxtok '$1==taxid{print $3; exit}' $fullnamelineage)

[ ! -s ${wdir}/taxidlist ] && make_taxidlist.sh $taxtok >${wdir}/taxidlist

# log the info and also store the taxid info in the settings file
msglog_module "$topid_name mitogenomes, taxid $taxtok, chosen as search targets. $num_mitogenomes $topid_name mitogenomes in the $mitodb_path db."

update_setting "taxid"   "$taxtok"
update_setting "taxname" "$topid_name"
set_mito_code_from_taxid "$taxtok"  # gets it from nodes.dmp based on taxtok and puts it in settings

update_setting "taxlineage" "$fullnamelineage"
update_setting "taxnodes" "$taxnodes"

update_setting "mitogenomes" $num_mitogenomes
