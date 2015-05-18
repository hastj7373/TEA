#!/bin/bash

. $(dirname $BASH_SOURCE)/config.sh

RR="$TEA_BASE/scripts/run_rid.sh"

if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` /path/to/demo/directory"
	exit 1
fi

RUN_DIR=$1
#XXX: SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850? | sed -e "s|$RUN_DIR/||g"`
SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850* | sed -e "s|$RUN_DIR/||g"`

echo "About to call RAMs for all of the following samples in parallel:"
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
for i in $SAMPLES; do
	cd $RUN_DIR/$i
	echo "Generating run.sh via run_rid.sh on sample $i.."
	$RR $RUN_DIR $i >& step2_rr.log
	echo "Running run.sh"
	./run.sh >& step2_run.log &
done
