#!/bin/bash

root_dir=/projects/tcga-data/

declare -a diseases=("Acute_Myeloid_Leukemia"\
					 "Liver_Hepatocellular_Carcinoma"\
					 "Adrenocortical_Carcinoma"\
					 "Lung_Adenocarcinoma"\
					 "Lung_Squamous_Cell_Carcinoma"\
					 "Lymphoid_Neoplasm_Diffuse_Large_B-cell_Lymphoma"\
					 "Bladder_Urothelial_Carcinoma"\
					 "Mesothelioma"\
					 "Brain_Lower_Grade_Glioma"\
					 "Ovarian_Serous_Cystadenocarcinoma"\
					 "Breast_Invasive_Carcinoma"\
					 "Pancreatic_Adenocarcinoma"\
					 "Pheochromocytoma_and_Paraganglioma"\
					 "Prostate_Adenocarcinoma"\
					 "Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma"\
					 "Rectum_Adenocarcinoma"\
					 "Cholangiocarcinoma"\
					 "Colon_Adenocarcinoma"\
					 "Sarcoma"\
					 "Skin_Cutaneous_Melanoma"\
					 "Esophageal_Carcinoma"\
					 "Stomach_Adenocarcinoma"\
					 "Testicular_Germ_Cell_Tumors"\
					 "Glioblastoma_Multiforme"\
					 "Thymoma"\
					 "Head_and_Neck_Squamous_Cell_Carcinoma"\
					 "Thyroid_Carcinoma"\
					 "Uterine_Carcinosarcoma"\
					 "Kidney_Chromophobe"\
					 "Uterine_Corpus_Endometrial_Carcinoma"\
					 "Kidney_Renal_Clear_Cell_Carcinoma"\
					 "Uveal_Melanoma"\
					 "Kidney_Renal_Papillary_Cell_Carcinoma")

all_files=""
for dt in "${diseases[@]}"
do
	for file in ${root_dir}/$dt/data/Bacterial_*.biom
	do
		if [ -f $file ] then;
			echo $file
		else
			echo "WARNING: BIOM does not exist for $dt"
		fi
		#all_files=$all_files,$file
	done
done