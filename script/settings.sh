#!/bin/bash

# set color vars that we use or might use here
Blue="\033[0;34m"; BBlue="\033[1;34m"; UBlue="\033[4;34m"; NC="\033[0m" # No Color (reset)

# let user see all the settings, a specific setting, or change a specific setting

this_dir=$(dirname $(realpath $0)) && source ${this_dir}/shared.sh && this=$(basename $0) # sets msg, is_number functions and usage among other things
[ -z $wdir ] && msglog_module "The hfmt_<num> working directory not found" && exit 2

[ -z $1 ] && $usage settings && exit

[ $1 = "all" ] && cat $wdir/settings.tsv && exit

[ ! -z $1 ] && [ -z $2 ] && grep $1 $wdir/settings.tsv && exit

while true; do
    read -p "Are you sure you want to set $1 to $2 [Y/N]? " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Y or N";;
    esac
done

# they said yes change it

update_setting_if_changed $1 $2
grep "^$1" $wdir/settings.tsv
