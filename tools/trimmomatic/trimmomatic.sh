#!/bin/sh
#
# Shell wrapper to run Trimmomatic jar file as a Galaxy tool
#echo Arguments:
count=0
for i in $@ ; do
    echo "*" $i
    ((count+=1))
done
echo $count

# figure out if inputs are gz formatted
format=`file $8`
if [ $4 == "PE" ]; then
  if [[ $format == *gzip* ]]; then
      echo $format
      ln -s $8 input1.fastq.gz
      ln -s $9 input2.fastq.gz
      in1="input1.fastq.gz"
      in2="input2.fastq.gz"
      set -- "${@:1:7}" $in1 $in2 "${@:10:count-1}"
  else
      echo "nonzipped PE format"
  fi
else
  if [[ $format == *gzip* ]]; then
      echo $format
      ln -s $8 input1.fastq.gz
      in1="input1.fastq.gz"
      set -- "${@:1:7}" $in1 "${@:9:count-1}"
  else
      echo "nonzipped PE format"
  fi
  echo "SE"
fi
echo java $@

java $@ 2>&1 | tee trimmomatic.log
status=$?
echo "Exit status: $status"
# Check for successful completion
if [ -z "$(tail -1 trimmomatic.log | grep "Completed successfully")" ] ; then
    echo "Trimmomatic did not finish successfully" >&2
    exit 1
fi
exit $status
##
#
