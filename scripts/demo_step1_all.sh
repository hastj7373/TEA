#!/bin/bash

. $(dirname $BASH_SOURCE)/config.sh

PP="$TEA_BASE/scripts/preprocess.sh"

if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` /path/to/demo/directory [bsub_command]"
	exit 1
fi

### BSUB = bsub -q whatever -J jobid
RUN_DIR=$1
#XXX: SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850? | sed -e "s|$RUN_DIR/||g"`
SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850* | sed -e "s|$RUN_DIR/||g"`

echo "About to preprocess all of the following samples in parallel:"
for i in $SAMPLES; do
	echo "    $i"
done
echo "These jobs will be run in the background!  Just because this script"
echo "has finished submitting them DOES NOT mean they are finished running."
echo "Monitor your jobs with 'ps' or some other tool to determine when they"
echo "have finished."
echo -n "Continue? (y/N) "
read X
if [ "x$X" != "xy" ]; then
	echo "Aborting"
	exit 1
fi
BSUB=""

for i in $SAMPLES; do
	cd $RUN_DIR/$i
	echo "Launching $PP on sample $i"
	$BSUB $PP $i.sorted.bam >& step1.log &
done
