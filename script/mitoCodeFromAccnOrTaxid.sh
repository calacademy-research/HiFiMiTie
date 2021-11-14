#!/bin/bash

# given a GenBank Accession ID or a numeric taxid, get its taxid from nucl_gb.accession2taxid and then its lineage from fullnamelineage.dmp

: '
2	Vertebrate
4	Mold
5	Invertebrate	Arthropoda
9	Echinoderm	starfish, sea urchins, sand dollars, and sea cucumbers
13	Ascidian	tunicates and sea squirts
14	Alternative-flatworm
'
# Protostomes are divided into the Ecdysozoa (e.g. arthropods, nematodes) and the Spiralia (e.g. molluscs, annelids, platyhelminths, and rotifers)

function msg { # write a msg to stderr
   >&2 echo -e $@
}

function is_number { # checks if argument is a positive integer, return value is true if a number false otherwise
   re='^[0-9]+$'
   [[ "$1" =~ $re ]]
}

dir=/ccg/db_sets/taxdump

accn2taxid=${dir}/nucl_gb.accession2taxid
lineage=${dir}/fullnamelineage.dmp

accn=$1

[ -z "$accn" ] && msg "\n    usage: codeFromAccn.sh <genbank accn id>\n" && exit 2
[ ! -f $accn2taxid ] && msg "\n    could not find: $accn2taxid\n"        && exit 3
[ ! -f $lineage ]    && msg "\n    could not find: $lineage\n"           && exit 4

if is_number $accn; then
   taxid=$accn
else
   line=$(grep "^$accn\b" $accn2taxid -m1)
   [ -z "$line" ] && msg "\n    $accn not found in $accn2taxid\n" && exit 5
   taxid=$(echo $line | awk '{print $3}')
fi

fullname=$(grep "^$taxid\b" $lineage -m1)

code=$(echo $fullname | awk 'BEGIN{code=5} # default is invertebrate
   index($0, "Tunicata;")                                 { code = 13; exit }
   index($0, "Vertebrata;") || index($0, "Chordata;")     { code=2; exit }
   index($0, "Arthropoda;") || index($0, "Panarthopoda;") { code=5; exit }
   index($0, "Echinodermata")                             { code=9; exit }
   END{print code}
')

echo -e "$code\tmito code of $accn $fullname"
