
# the contents of this file are sourced everytime a script is run with hifimitie.sh

# be careful not to include anything that stops scripts from running
# one use is to set or add to the PATH or PYTHONPATH for file or module visibility

#-------------------------------------------------------------------------------------------------#

# we are using this since mito_analyze.py might not have the Bio package accessible
# longer term should convert to python3 and encourage installation of the Bio package

Bio_module_loc="/home/jhenderson/.local/lib/python2.7/site-packages"
[[ "$PYTHONPATH" == *"$Bio_module_loc"* ]] || export PYTHONPATH=$PYTHONPATH:$Bio_module_loc

#-------------------------------------------------------------------------------------------------#
