#!/bin/bash
#

POPNAME=$1
DESTINATION_DIRECTORY='../data/accepted_clones'

# clean up any previous files
rm $DESTINATION_DIRECTORY'/'${POPNAME}'_clone_list.tsv'

#compile all lineage for which first flag is "A" or "B"
for i in {1..10}; do
	file=../data/flags/${POPNAME}-BC${i}_flags.tsv
	echo ${file}
	awk '((index(substr($2,1,1),"A") != 0) || (index(substr($2,1,1),"B") != 0) ) {print $1}' ${file} >> $DESTINATION_DIRECTORY'/'${POPNAME}'_clone_list.tsv'
done

