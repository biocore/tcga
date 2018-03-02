#!/bin/bash

# e.g., /home/evko1434/tcga
TCGA_DATA_DIR=$1
ARIA_INSTALL_DIR="/home/evko1434/install_dir/aria2-1.28.0-install/bin/aria2c"

for FILENAME in $TCGA_DATA_DIR/*/download_links.txt; do
	DIR=$(dirname "${FILENAME}")
	OUTPUT_DIR="$DIR"/data
	mkdir -p $OUTPUT_DIR
	echo "${ARIA_INSTALL_DIR} -i ${FILENAME} -d ${OUTPUT_DIR}" | qsub -l nodes=1:ppn=1 -l walltime=48:00:00 -l pmem=4gb -l mem=4gb -k oe
done