#!/bin/bash

# User must change this value to point to the root of their installation
. $(dirname $BASH_SOURCE)/config.sh

# Do not change anything below this line
################################################################################

if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` sorted_input_bam_file"
	exit 1
fi
INPUT_FILE=$1

OUT_DIR="preprocess"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TEA_BASE/preprocess/lib"
BIN_DIR="$TEA_BASE/preprocess/bin"
DB_DIR="$TEA_BASE/preprocess/db"

echo "Creating output directory $OUT_DIR"
mkdir -p $OUT_DIR
echo "Linking to input"
ln -s ../$INPUT_FILE $OUT_DIR/$INPUT_FILE

S_OPT=
if [ $(echo $INPUT_FILE | grep -c 35bp) != "0" ]; then
	echo "Short read file detected."
	S_OPT="-l 0"  # Do not generate and realign clipped reads
fi

echo "=======================1======================="
perl $BIN_DIR/pp1.pl $S_OPT -s 20 -k 1500 -q 5 -t 2 -P is -b $OUT_DIR/$INPUT_FILE -R $OUT_DIR/`basename $INPUT_FILE .bam`.black.rg -A $DB_DIR/hg18/hg18.fasta.fai

echo "=======================2======================="
perl $BIN_DIR/pp1.pl $S_OPT -s 20 -k 1500 -q 5 -t 2 -P cl -b $OUT_DIR/$INPUT_FILE  -R $OUT_DIR/`basename $INPUT_FILE .bam`.black.rg -I $DB_DIR/hg18/hg18_bwa_idx/hg18.fasta -A $DB_DIR/hg18/hg18.fasta.fai

echo "=======================3======================="
perl $BIN_DIR/pp2.pl $S_OPT -s 20 -t 10 -p 3 -o 1 -P dc -b $OUT_DIR/$INPUT_FILE -R $OUT_DIR/`basename $INPUT_FILE .bam`.black.rg
