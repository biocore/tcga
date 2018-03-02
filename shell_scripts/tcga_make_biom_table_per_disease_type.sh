#!/bin/bash

#-----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

root_dir=/projects/tcga-data
output_dir=/projects/tcga-data/biom_tables_kraken_unfiltered
declare -a diseases=("Acute_Myeloid_Leukemia"\
					 "Adrenocortical_Carcinoma"\
					 "Bladder_Urothelial_Carcinoma"\
					 "Brain_Lower_Grade_Glioma"\
					 "Breast_Invasive_Carcinoma"\
					 "Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma"\
					 "Cholangiocarcinoma"\
					 "Colon_Adenocarcinoma"\
					 "Esophageal_Carcinoma"\
					 "Glioblastoma_Multiforme"\
					 "Head_and_Neck_Squamous_Cell_Carcinoma"\
					 "Kidney_Chromophobe"\
					 "Kidney_Renal_Clear_Cell_Carcinoma"\
					 "Kidney_Renal_Papillary_Cell_Carcinoma" \
					 "Liver_Hepatocellular_Carcinoma"
					 "Lung_Adenocarcinoma"\
					 "Lung_Squamous_Cell_Carcinoma"\
					 "Lymphoid_Neoplasm_Diffuse_Large_B-cell_Lymphoma"\
					 "Mesothelioma"\
					 "Ovarian_Serous_Cystadenocarcinoma"\
					 "Pancreatic_Adenocarcinoma"\
					 "Pheochromocytoma_and_Paraganglioma"\
					 "Prostate_Adenocarcinoma"\
					 "Rectum_Adenocarcinoma"\
					 "Sarcoma"\
					 "Skin_Cutaneous_Melanoma"\
					 "Stomach_Adenocarcinoma"\
					 "Testicular_Germ_Cell_Tumors"\
					 "Thymoma"\
					 "Thyroid_Carcinoma"\
					 "Uterine_Carcinosarcoma"\
					 "Uterine_Corpus_Endometrial_Carcinoma"\
					 "Uveal_Melanoma")

all_files=""
for dt in "${diseases[@]}"
do
	commands_fp=${output_dir}/${dt}/bacteria/commands.txt
	echo "Disease: $dt" >> ${commands_fp}
	for file in ${root_dir}/$dt/data/Bacterial_*.biom
	do
		if [ ! -f $file ]
		then
			echo "WARNING: BIOM does not exist for $dt" >> ${commands_fp}
		else
			echo $file >> ${commands_fp}
		fi
		all_files=$all_files,$file
	done
	all_files="${all_files:1}"
	rm -r ${output_dir}/$dt/bacteria
	mkdir -p ${output_dir}/$dt/bacteria
	echo "command: merge_otu_tables.py -i ${all_files} -o ${output_dir}/${dt}/bacteria/Bacterial_${dt}_all.biom" >> ${commands_fp}
	echo "unset PYTHONPATH; \
		  source activate qiime1; \
		  merge_otu_tables.py -i ${all_files} -o ${output_dir}/${dt}/bacteria/Bacterial_${dt}_all.biom; \
		  biom summarize-table -i ${output_dir}/${dt}/bacteria/Bacterial_${dt}_all.biom -o \
		  ${output_dir}/${dt}/bacteria/Bacterial_${dt}_all_table_summary.txt" | qsub -l nodes=1:ppn=1 -l walltime=1:00:00 -N build_biom_unfiltered_bacteria_${dt} -m abe -M jenya.kopylov@gmail.com
	all_files=""
done
