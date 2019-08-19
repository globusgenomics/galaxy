#!/bin/bash

PICARD_HOME="/nfs/software/picard-tools-1.17"

java -jar "$PICARD_HOME/SamToFastq.jar" "$@"
