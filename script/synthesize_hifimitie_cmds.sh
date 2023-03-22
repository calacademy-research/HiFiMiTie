#!/bin/bash

# using the info in the hfmt $wdir dir create what the init and run commands looks like

function msg { echo -e "$@" >/dev/stderr; }

function synthesize_init_cmd {
   [ ! -s $wdir/settings.tsv ]         && msg "No $wdir/settings.tsv file found"    && exit 1
   [ ! -s $wdir/files_to_search.fofn ] && msg "No $wdir/files_to_search.fofn found" && exit 2

   taxid=$(grep "^taxid\s" $wdir/settings.tsv | cut -f2)
   [ -z "$taxid" ] && msg "Could not find taxid in settings.tsv file" && exit 3

   taxname=$(grep "^taxname\s" $wdir/settings.tsv | cut -f2)
   [ ! -z "taxname" ] && echo -e "# $taxname taxid $taxid used for mitogenome search and CR ending flank info\n"

   prefix="hifimitie init $taxid "
   len=$(echo $prefix | awk '{print length($0);exit}')

   files=$(awk -v len=$len '
      BEGIN{cm="%-" len "s"; spaces=sprintf(cm, " ")}

      /^$/ || /^#/ { next }
      { f++ }

      NR==1{ files=$1 }
      NR>1 { files = files " \\\n " spaces $1 }
      END { print files }
   ' $wdir/files_to_search.fofn)

   echo -e $prefix "$files"
}

function synthesize_run_cmd {
   threads=$(grep "^threads\s" $wdir/settings.tsv | cut -f2)
   [ -z "threads" ] && msg "Could not find threads in settings.tsv file" && exit 4

   echo hifimitie run -t $threads
}

function set_wdir {
   wdir=$1
   [ ! -z "$wdir" ] && return  # accept what was passed in

   wdir=$(ls -d * | grep -m 1 "^hfmt_[0-9]*")

   [ ! -d $wdir ] && msg "\n    Could not find hfmt_<num> directory\n" && exit 5
}


set_wdir $1

synthesize_init_cmd
echo ""
synthesize_run_cmd
