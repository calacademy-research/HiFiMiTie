#!/bin/bash

# separate file so we can show extensive usage / help info if we want to, without getting in way of main script

cmd=$1

if [ -z $cmd ]; then
    echo -e '
    HiFiMiTie version 0.001

    Usage: hifimitie <command> [options]

    Commands:
      -- Setup
         init           create working directory and optionally set taxid and HiFi file(s)
         gettaxid       search all mitogenomes or only those associated with a certain taxid (optionally view taxonomy ids based on term)
         addfiles       add some or all of the HiFi files with which to begin the search
         clearfiles     remove the files already in the fofn file list

      -- Execution
         run            run the default pipeline of ...,...,..
         trf_analysis   perform tandem repeat analysis on control region excerpts
'
elif [ $cmd == "init" ]; then
    echo -e '
    Usage: hifimitie init [<taxid number | taxonomy term>] [HiFi file]...

    creates the hfmt working directory if necessary
    places fastq or fasta HiFi files into files_to_search.fofn in the hfmt directory

    taxid number is used to create a taxidlist for blast searches of mitogenomes (use 0 to search them all)
    using a term will show the taxonomy list info. when you find the correct number, enter it after quitting the view
'
else
   echo -e "(usage) command \"$cmd\" not found"
fi
