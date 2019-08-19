#!/bin/bash


while getopts "i:o:" optionName; do
case "$optionName" in

i) INPUT="$OPTARG";;
o) OUTPUT="$OPTARG";;

esac
done

rm $OUTPUT

ln -s $INPUT $OUTPUT
