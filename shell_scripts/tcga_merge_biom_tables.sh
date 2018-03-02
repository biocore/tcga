#!/bin/bash

root_dir=/projects/tcga-data

tables2merge=""
for file in ${root_dir}/*/data/Bacterial_*.biom
do
	echo $file
	tables2merge="$file,$tables2merge"
done
# remove last ','
tables2merge=${tables2merge%?}

echo $tables2merge
echo "source activate qiime1; unset PYTHONPATH; merge_otu_tables.py -i $tables2merge -o ${root_dir}/Bacterial_all_cancer_types.biom" | qsub -l nodes=1:ppn=1 -l walltime=24:00:00 -N tcga_BIOM_bacteria

