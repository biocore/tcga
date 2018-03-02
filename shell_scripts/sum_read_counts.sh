#!/bin/bash

TCGA_DATA_DIR=$1
OUTPUT_FP=$2
READ_COUNTS_PY="/home/evko1434/tcga/scripts/sum_read_counts.py"

#for DIR in $TCGA_DATA_DIR/*; do
#    if [ -d "$DIR" ]; then
#        BN=$(basename "${DIR}")
#        OUTPUT_FPP="${OUTPUT_FP}_$BN"
#        echo $BN >> ${OUTPUT_FPP}
#        echo "python ${READ_COUNTS_PY} --data-dp $DIR/data/ >> ${OUTPUT_FPP}" | qsub -l nodes=1:ppn=1 -l walltime=24:00:00 -N tcga_$BN
#    fi
#done

BN=$(basename "${TCGA_DATA_DIR}")
OUTPUT_FPP="${OUTPUT_FP}_$BN"
echo $BN >> ${OUTPUT_FPP}
echo "python ${READ_COUNTS_PY} --data-dp $TCGA_DATA_DIR >> ${OUTPUT_FPP}" | qsub -l nodes=1:ppn=1 -l walltime=48:00:00 -N 01022018_tcga_read_count_$BN
