#!/bin/bash
#$ -N ooa_naive
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _ooa_naive_err.txt
#$ -o _ooa_naive_log.txt
#$ -cwd


TH="1"

cd ./results_naive


for METHOD in `echo "fgt dgt fgt" | tr " " "\n" | shuf`
do
	DATA="truH errH"

	if [ "$METHOD" == "fgt" ]
	then
		DATA="truH errH truP errP"
	fi

	for MODE in `echo ${DATA} | tr " " "\n" | shuf`
	do

		HMM=""
		if [ "$METHOD" == "hmm" ] && [ "$MODE" == "errH" ]
		then
			HMM="--hmm ../HMM.initial.prob.txt ../HMM.emission.prob.txt"
		fi

		for MUT in `echo 0 1 | tr " " "\n" | shuf`
		do
			for REC in `echo 0 1 | tr " " "\n" | shuf`
			do
				if [ "$MUT" == "0" ] && [ "$REC" == "0" ]; then continue; fi

				STR="${METHOD}_${MODE}_${MUT}_${REC}"

				echo $STR

				for FILE in `ls ../packs_naive/expos.* | shuf`
				do
					PACK=`basename ${FILE}`
					PACK=`echo $PACK | cut -d '.' -f 3`

					FLAG_RUNS="RUNS__${STR}_${PACK}"
					FLAG_DONE="DONE__${STR}_${PACK}"

					if [ -e "../logs_naive/${FLAG_RUNS}" ]; then continue; fi
					if [ -e "../logs_naive/${FLAG_DONE}" ]; then continue; fi

					> ../logs_naive/${FLAG_RUNS}

					echo $PACK

					../rvage age -i ../${MODE}/${MODE}.bin -o $STR.${PACK} -m ${METHOD} --positions ../packs_naive/expos.10000.${PACK}.txt -t ${TH} --useMutClock ${MUT} --useRecClock ${REC} --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08 ${HMM}

					mv ../logs_naive/${FLAG_RUNS} ../logs_naive/${FLAG_DONE}

				done

			done
		done
	done
done
