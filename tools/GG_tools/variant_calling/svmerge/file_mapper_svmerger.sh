#!/bin/bash
#take care of the mapper args
while [ $# -gt 0 ]; do
  case $1 in
    -location)          location=$2;;
    -grep_file)         search_grep=$2;;
    *) echo "$0: bad mapper args" 1>&2
       exit 1;;
  esac
  shift 2
done

find ${location} -maxdepth 1 -type f -iname "${search_grep}.*" -print | sort | awk '{printf "[%d] %s\n", NR-1, $1}'

