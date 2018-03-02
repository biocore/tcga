#!/bin/bash

#-----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

root_dir=/projects/tcga-data/clean_Kraken_MPA
output_dir=/projects/tcga-data/filtered_biom_tables
disease_type="Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma"
scripts=/home/evko1434/tcga/scripts

for file in ${root_dir}/Viral_${disease_type}_Classified_Kraken_mpa_report_*_clean.txt
do
	bn=${file##*/}
	bn=${bn%.*}
	echo "unset PYTHONPATH; \
		  source activate parse_kraken_to_biom; \
		  python $scripts/parse_kraken_to_biom.py \
				 --kraken-translate-report-fp $file \
				 --taxonomic-rank genus \
				 --biom-output-fp ${output_dir}/${bn}.biom" | qsub -l nodes=1:ppn=1 -l pmem=5Gb -l mem=5Gb -l walltime=48:00:00 -N build_biom_${bn} -m abe -M jenya.kopylov@gmail.com
done