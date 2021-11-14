#!/bin/bash
# convenient setting of ANSI escape codes into vars
# use source clr_vars.sh in your script to include these
# or copy the contents into your script
Black="\033[0;30m";  BBlack="\033[1;30m";  UBlack="\033[4;30m"
Red="\033[0;31m";    BRed="\033[1;31m";    URed="\033[4;31m"
Green="\033[0;32m";  BGreen="\033[1;32m";  UGreen="\033[4;32m"
Yellow="\033[0;33m"; BYellow="\033[1;33m"; UYellow="\033[4;33m"
Blue="\033[0;34m";   BBlue="\033[1;34m";   UBlue="\033[4;34m"
Purple="\033[0;35m"; BPurple="\033[1;35m"; UPurple="\033[4;35m"
Cyan="\033[0;36m";   BCyan="\033[1;36m";   UCyan="\033[4;36m"
White="\033[0;37m";  BWhite="\033[1;37m";  UWhite="\033[4;37m"
Undl="\033[4m"; Bold="\033[1m"
NC="\033[0m" # No Color (reset)
