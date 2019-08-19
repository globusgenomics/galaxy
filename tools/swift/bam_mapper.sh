#!/bin/bash

find -L $2 -maxdepth 1 -type f -iname "*.bam" -print | sort | awk '{printf "[%d] %s\n", NR-1, $1}' 
