#!/bin/bash

cd /ccg/blastdbs

rm -f mito.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz
tar xzvf mito.tar.gz

blastdbcmd -db mito -tax_info >mito_taxids.tsv
blastdbcmd -db mito -entry all | grep "^>" >mito_db_entries.txt

echo
blastdbcmd -db mito -info
