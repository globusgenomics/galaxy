#!/bin/bash
# Author: Edward Hills and James Boocock
# Date: 8/12/2011
#
# Script takes and runs the variant effect predictor
# from ensemble from the arguments given.
#
# Inputs
# $1 = input file
# $2 = Sift
# $3 = Polyphen
# $4 = Condel
# $5 = Regulatory
# $6 = Protein
# $7 = HGNC
# $8 = CCDS
# $9 = Most Severe
# $10 = Summary
# $11 = Per Gene
# $12 = Coding Only
# $13 = Check Existing
# $14 = Check Alleles

ENSEMBL_RUN_SCRIPT=""
NUM_SAMPLES=$#

if [ $NUM_SAMPLES != 1 ]
then
	for (( i=2; i <= $NUM_SAMPLES; i++ ))
	do
		NUM=${i}
		eval INPUT=\${${NUM}}
		if [ "${i}" -lt "5" ]; then
			if [ "${i}" != "none" ]; then
			    ENSEMBL_RUN_SCRIPT="${ENSEMBL_RUN_SCRIPT} --${INPUT} b"
			fi
		else
		    ENSEMBL_RUN_SCRIPT="${ENSEMBL_RUN_SCRIPT} --${INPUT}"
		fi
	done
# call actual script
        vep -i $1 --cache --dir /mnt/galaxyIndices/genomes/vep/ --force_overwrite --merged --port 3337 --everything -o ./ensemble-TMP.tmp

#	vep -i $1 -o ./ensemble-TMP.tmp $ENSEMBL_RUN_SCRIPT --merged --cache --port 3337 --dir "/mnt/galaxyIndices/genomes/vep/" --hgvs --force_overwrite --buffer 50000 --fork 2

	cat ./ensemble-TMP.tmp
	rm -f ./ensemble-TMP.tmp

else # call defaults 
        vep -i $1 --cache --dir /mnt/galaxyIndices/genomes/vep/ --force_overwrite --merged --port 3337 --everything -o ./ensemble-TMP.tmp
#	vep -i $1 -o ./ensemble-TMP.tmp --check_existing --gene \
#                        --cache --merged --port 3337 --dir "/mnt/galaxyIndices/genomes/vep/" \
#                       --poly b --sift b --hgvs --force_overwrite  \
#                       --buffer 50000 --fork 2
fi

exit 0
