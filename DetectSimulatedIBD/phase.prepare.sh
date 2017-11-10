#!/bin/bash

PREFIX=$1
SUFFIX="RData"

ROOT=`pwd`

for FILE in `find ./ -iname "${PREFIX}*${SUFFIX}"`;
do

	cd `dirname $FILE`

	BASE=`basename $FILE`
	BASE=${BASE%.RData}

	echo "###"
	echo $BASE
	echo "###"

	Rscript ${ROOT}/phase.prepare.R $BASE &

	echo ""

	cd $ROOT

done

wait
