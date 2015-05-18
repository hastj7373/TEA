#!/bin/bash
# Runs Alice's rid.pl script and produces a run.sh bash script.

# User must set set this to the root directory of the Tea installation
. $(dirname $BASH_SOURCE)/config.sh

if [ $# -ne 2 ]; then
	echo "Usage: `basename $0` /path/to/run/dir sample_name"
	exit 1;
fi

RUN_DIR=$1
SAMPLE=$2

# Args:
#  -S do not look for/generate .cl. files.  Have to use this now because
#     the BAM merge command fails on .cl. headers.
#  -b source directory.  contains files produced by preprocess.
#  -g the original BAM file that was previously fed to preprocess.
#  -P additional directories to add to $PATH
if [ ! -d $RUN_DIR/$SAMPLE/preprocess ]; then
	echo "ERROR: for sample $SAMPLE"
	echo "Could not find preprocess output directory $RUN_DIR/$SAMPLE/preprocess"
	echo "aborting"
	exit 1
fi

S_OPT=
if [ $(echo $SAMPLE | grep -c 35bp) != "0" ]; then
        echo "Short read file detected."
        S_OPT="-S"  # Do not look for softclipped reads
fi

$TEA_BASE/scripts/rid.pl $S_OPT -d $RUN_DIR -s $RUN_DIR/$SAMPLE/scratch/ -b $RUN_DIR/$SAMPLE/preprocess -g $RUN_DIR/$SAMPLE/$SAMPLE.sorted.bam $ORIG_BAM -P $TEA_BASE/scripts $SAMPLE > $RUN_DIR/$SAMPLE/run.sh
chmod +x $RUN_DIR/$SAMPLE/run.sh
