#!/bin/bash

source $(dirname $(realpath $0))/shared.sh

#src_dir=$(dirname $(realpath $0)) && source ${src_dir}/shared.sh

echo $usage

function firstCharUC_resp {
   resp=${1:0:1}; resp=$(echo $resp |tr '[:lower:]' '[:upper:]')
}

function test {
   read -p "test: " usertest
   firstCharUC_resp $usertest
   echo $resp
}

# test
