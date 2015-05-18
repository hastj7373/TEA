#!/bin/bash

. $(dirname $BASH_SOURCE)/config.sh

RID_PATH="$TEA_BASE/R/rid.r"

if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` /path/to/run/dir"
	exit 1
fi

echo "Running an R analysis in the background!  Just because this script"
echo "has finished DOES NOT mean the analysis has finished running."
echo "Monitor the R job with 'ps' or some other tool to determine when it"
echo "has finished."

RUN_DIR=$1
CHRL='c("chr21", "chr22")'
#CHRL='hchrl'   # For all chromosomes

R --no-save --no-restore --no-environ << MYEOF >& $RUN_DIR/step3.log &
source("$RID_PATH")
call.germline.batch(dir="$RUN_DIR", chrl=$CHRL, verbose=TRUE)
call.somatic.batch(dir="$RUN_DIR", chrl=$CHRL, verbose=TRUE, min.arr=0.5)
MYEOF
