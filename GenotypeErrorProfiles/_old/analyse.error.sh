#!/bin/bash

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

ROOT=`pwd`

SCRIPT=`echo "${ROOT}/analyse.error.R"`

for FILE in `find -iname 'profile.*.RData'`
do
	cd `dirname $FILE`

	EXEC=`basename $FILE`

	echo "###"
	echo $EXEC
	echo "###"

	Rscript $SCRIPT $EXEC

	cd $ROOT
done

IFS=$SAVEIFS

