#!/bin/bash

find $2 -maxdepth 1 -type f -iname "*.mzML" -print | sort | awk '{printf "[%d] %s\n", NR-1, $1}' 
