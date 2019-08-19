#!/bin/bash

# crash: Report a problem and exit
crash()
{
    MSG=$1
    echo ${MSG}  >&2
    exit 1
}

# Verify an argument is not null
verify_not_null()
{
   argname=$1; shift
   if [ _$1 != _ ]; then
      return 0;
   else
      echo $0: value for $argname can not be null
      exit 1
   fi
}

# Process what we know and pass the rest to Swift
SWIFTARGS=""
while [ $# -gt 0 ]
do
   case "$1" in
      -level) export LEVEL=$2; verify_not_null level $LEVEL; shift ;;
      -template) export TEMPLATE=$2; verify_not_null template $TEMPLATE; shift;;
      -work) export WORK=$2; verify_not_null work $WORK; shift;;
      -project) export PROJECT=$2; verify_not_null project $PROJECT; shift;;
      -queue) export QUEUE=$2; verify_not_null queue $QUEUE; shift;;
      -o=*) export OUTPUT=`echo $1| cut -d'=' -f2`; verify_not_null output $OUTPUT; SWIFTARGS="$SWIFTARGS $1" ;;
      -jobthrottle) export JOBTHROTTLE=$2; verify_not_null jobthrottle $JOBTHROTTLE; shift;;
      -log) export LOG=$2; verify_not_null log $LOG; shift;;
       *) SWIFTARGS="$SWIFTARGS $1";;
   esac
   shift
done

# Verify level
if [ -z "$LEVEL" ]; then
   crash "Level not specified. Use -level <value>"
fi

case "$LEVEL" in
   0) ;;
   1) ;;
   2) ;;
   3) ;;
   *) crash "Unknown level $LEVEL";;
esac

# Verify a sites template exists
if [ -z "$TEMPLATE" ]; then
   crash "Template not specified. Use -template <name>"
fi

# Create Swift configuration files
pushd /nfs/software/galaxy-globus/swift_scripts/qc > /dev/null 2>&1
export PATH=$PATH:/nfs/software/swift-0.92.1/bin
export SWIFT_HOME=/nfs/software/swift-0.92.1
gensites -L . $TEMPLATE -p cf > sites.xml
export SWIFT_OUTPUT_DIR=`mktemp -d -p /glusterfs/galaxy-data/tmp`

# Run Swift
rm output > /dev/null 2>&1
rm output2 > /dev/null 2>&1
ln -s $SWIFT_OUTPUT_DIR output > /dev/null 2>&1
ln -s $SWIFT_OUTPUT_DIR output2 > /dev/null 2>&1
#echo swift -sites.file sites.xml -tc.file tc.data -config cf -cdm.file fs.alldirect cga.swift -level=$LEVEL $SWIFTARGS >> $LOG 2>&1
/usr/bin/time swift -sites.file sites.xml -tc.file tc.data -config cf -cdm.file fs.alldirect qc.swift -level=$LEVEL $SWIFTARGS >> $LOG 2>&1
#echo ./create_cga_dataset.py -o $OUTPUT -p `dirname $OUTPUT` -i output >> $LOG 2>&1
./create_cga_dataset.py -o $OUTPUT -p `dirname $OUTPUT` -i `ls -ld output|awk {'print $10'}` >> $LOG 2>&1
#popd > /dev/null 2>&1
                       
