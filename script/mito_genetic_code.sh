#!/bin/bash

# given a taxonomy id show the mitochondrial genetic code

codefile=nodes.dmp
codedir=/ccg/db_sets/taxdump
codepath=$codedir/$codefile

mito_field=9
default_code=5

function msg {
   >&2 echo -e "$@"
   exit 1
}

function usage {
   msg "
       mito_genetic_code.sh <taxid_number>  [<optional dir for $codefile or $codefile fullpath>]
"
}

[ -z $1 ] && usage
taxid=$1

if [ -d "$2" ]; then
   codepath=$2/$codefile
elif [ -s "$2" ]; then
   codepath=$2
fi

# node_line should look something like this: 6073|6072|phylum||1|1|1|1|4|0|0|0||||0|0|1|
node_line=$(grep "^$taxid\s" $codepath | sed "s/\s*//g")

# look in nodes.dmp file for taxid and return code, or if not there return default_code
awk -v default_code=$default_code -v mito_field=$mito_field ' {

   flds = split($0, ar, "|")
   if (flds >= mito_field)
      code = int(ar[9])
   else
      code = int(default_code)

   if (code < 1 || code > 35)
      code = default_code

   print code

   exit # only process the first line we are given (which is probably all we are given but just in case)
} ' <(echo $node_line)
