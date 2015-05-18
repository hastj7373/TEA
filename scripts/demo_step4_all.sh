#!/bin/bash

. $(dirname $BASH_SOURCE)/config.sh

CAP3="/data/home/jluquette/alice/CAP3/cap3"

if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` /path/to/demo/directory"
	exit 1
fi

RA="$TEA_BASE/scripts/ra.pl"
RUN_DIR=$1
#XXX: SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850? | sed -e "s|$RUN_DIR/||g"`
SAMPLES=`ls -1d $RUN_DIR/*_cancer $RUN_DIR/*_normal $RUN_DIR/na1850* | sed -e "s|$RUN_DIR/||g"`

echo "About to generate contigs for all of the following samples in parallel:"
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
	echo "Launching $RA on sample $i"
	if [ ! -d contig2 ]; then
		mkdir contig2
	fi

	if [ -d call ]; then
		for f in `ls call`; do
			$RA contig2 -a $CAP3 -d $RUN_DIR -o `pwd`/contig2 -t `pwd`/scratch call/$f >& step4.log &
		done
	fi
	if [ -d call3 ]; then
		for f in `ls call3`; do
			$RA contig2 -a $CAP3 -d $RUN_DIR -o `pwd`/contig2 -t `pwd`/scratch call3/$f >& step4.log &
		done
	fi
done
