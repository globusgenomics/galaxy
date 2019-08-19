#!/bin/bash
file=$1
sed -i 's/value"/value="/g' $file
