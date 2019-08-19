#!/bin/bash

inputfile=$1
col=$2
addremove=$3
outputfile=$4

echo "args: $@"
echo "inputfile: $inputfile"
echo "column: $column"
echo "addremove: $addremove"
echo "outputfile: $outputfile"

#get column number
column=`expr match "$col" '\([0-9]*\)'`  
echo "colnumber: $column"

if [ $addremove == "add" ]
then
	echo "adding prefix to column $column"
	awk 'BEGIN{
		FS="\t"
		OFS="\t"
		c="'"$column"'"
	}{
		if (index($0,"#")!=1){
			$c="chr"$c			
		}
		print $0
		
	}END{}' $inputfile > $outputfile

else	#remove prefix
	echo "removing prefix from column $column"
	awk 'BEGIN{
		FS="\t"
		OFS="\t"
		c="'"$column"'"
	}{
		if (FNR>1 && index($0,"#")!=1){
			$c=substr($c,4)			
		}
		print $0
	
	}END{}' $inputfile > $outputfile
fi

echo "inputfile: "
head -5 $inputfile

echo "outputfile: "
head -5 $outputfile
