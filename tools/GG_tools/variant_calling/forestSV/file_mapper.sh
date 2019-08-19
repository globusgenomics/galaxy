#!/bin/bash

find $2 -maxdepth 1 -type f -iname "*.txt" -print | sort | awk '{printf "[%d] %s\n", NR-1, $1}' 

