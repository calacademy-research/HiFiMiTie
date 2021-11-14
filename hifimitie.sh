#!/bin/bash

dir=$(dirname $(realpath $0))

# load these preamble settings before sourcing the shared functions and var settings and before executing the script
[ -s ${dir}/script_preamble.sh ] && source ${dir}/script_preamble.sh

# load the shared functions and var settings
[ -f ${dir}/shared.sh ] && source ${dir}/shared.sh
[ -z "$shared_loaded" ] && source ${dir}/script/shared.sh
[ -z "$shared_loaded" ] && echo -e "Could not load the shared functions needed from the shared.sh file."

function is_command {
   awk -v chk=$cmd 'BEGIN{
      cmds="init\trun\tcheck\tsettings\tgettaxid\taddfiles\tclearfiles\tblast_to_mito\tpull_fofn_cand_recs\tselect_mito_features\tblast_features\t"
      cmds = cmds "rrna_cmsearch\tsplit_recs_into_sets\tassemble_mito\tassemble_w_megahit\tassemble_w_mafft_and_consensus\tcompare_assemblies\tassemble_CR\tcomplete\t"
      cmds = cmds "make_mito_rec_cand_tsv\trna_search\ttest\t"

      cmdchk = chk "\t"
      print index(cmds, cmdchk) ? chk : ""
   }'
}

cmd=$1

[ -z $cmd ] && [ ! -f $usage ] && echo "incomplete installation: usage.sh not found" && exit 1
[ -z $cmd ] && $usage && exit 0

[ -z $(is_command) ] && echo -e "command \"$cmd\" not found" && exit 1

export hifimitie_pgmdir=$dir

# look in src_dir first, then if not there look in script_dir, we'll add a .sh if needed and also be ok if name is prefixed with a 2 digit number followed by underscore
script=$(get_script $cmd)
shift

$script $@
