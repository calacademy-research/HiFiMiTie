#!/bin/bash

# function to set blast major and mino version number variables
# we can use this iwth version 2.12 or greater to use the new -mt_mode 1 option
# this will use the threads among the queries instead of partitioning the database search
# with a -taxidlist in use this is probably faster

function set_blast_major_minor {
   blast_version=$(blastn -version | awk '/[1-9]/{sub("[^\\.0-9]*",""); print; exit}')
   blast_major=$(echo $blast_version | awk '{split($0,ar,"\\."); major=int(ar[1]); print major; exit}')
   blast_minor=$(echo $blast_version | awk '{split($0,ar,"\\."); minor=int(ar[2]); print minor; exit}')
}

function set_blast_extra_args {
   unset blast_extra_args
   exargs="-mt_mode 1"

   set_blast_major_minor
   [ $threads -lt 4 ] && return
   [ $blast_major -gt 2 ] && blast_extra_args=$exargs
   [ $blast_major -eq 2 ] && [ $blast_minor -ge 12 ] && blast_extra_args=$exargs

   unset exargs
}

threads=32

set_blast_extra_args
echo extra args "|"$blast_extra_args"|"


set_blast_major_minor

tst=1

if [ "$tst" -eq 1 ]; then
   echo $blast_version
   echo $blast_major
   echo $blast_minor
fi
