#!/bin/bash

# Cancer type dir
CANCER_TYPE_DIR=$1
STATS_DIR=${CANCER_TYPE_DIR}/stats
SCRIPTS_DIR=/home/evko1434/tcga/scripts
ARIA2C_FP=/home/evko1434/install_dir/aria2-1.28.0-install/bin/aria2c
DATE="21062017"

# Generate download links
echo "Generate download links"
python ${SCRIPTS_DIR}/cgc_get_download_links.py --yaml-fp ${SCRIPTS_DIR}/yaml_test_config.yaml --download-links-fp ${CANCER_TYPE_DIR}/dl_stats_$DATE.txt
# Download files
echo "Download files"
mkdir -p ${STATS_DIR}
${ARIA2C_FP} -i ${CANCER_TYPE_DIR}/dl_stats_$DATE.txt -d ${STATS_DIR}
# Compute stats
echo "Compute stats"
python ${SCRIPTS_DIR}/sum_stats_counts.py --stats-dp ${STATS_DIR}
# Double check
echo "Verify stats"
python ${SCRIPTS_DIR}/solve_Stomach_Adenocarcinoma.py --yaml-fp ${SCRIPTS_DIR}/yaml_test_config.yaml --data-dp ${STATS_DIR}