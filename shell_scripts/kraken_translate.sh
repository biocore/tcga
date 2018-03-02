#!/bin/bash

#disease_type="Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma"
#disease_type="Sarcoma"
#disease_type="Bladder_Urothelial_Carcinoma"
#disease_type="Prostate_Adenocarcinoma"
#disease_type="Uterine_Corpus_Endometrial_Carcinoma"
#disease_type="Breast_Invasive_Carcinoma"
#disease_type="Acute_Myeloid_Leukemia"
#disease_type="Skin_Cutaneous_Melanoma"
#disease_type="Brain_Lower_Grade_Glioma"
#disease_type="Kidney_Renal_Clear_Cell_Carcinoma"
#disease_type="Lung_Squamous_Cell_Carcinoma"
#disease_type="Lung_Adenocarcinoma"
#disease_type="Kidney_Renal_Papillary_Cell_Carcinoma"
#disease_type="Thyroid_Carcinoma"
#disease_type="Colon_Adenocarcinoma"
#disease_type="Glioblastoma_Multiforme"
#disease_type="Stomach_Adenocarcinoma"
#disease_type="Kidney_Chromophobe"
#disease_type="Head_and_Neck_Squamous_Cell_Carcinoma"
#disease_type="Liver_Hepatocellular_Carcinoma"
disease_type="Ovarian_Serous_Cystadenocarcinoma"
root_dir=/projects/tcga-data/${disease_type}/data
db_fp=/databases/repophlan/05_27_2016/kraken-index/viruses-dust/tcgadb_vir
#output_dir=/projects/tcga-data/bwa_alignments/viral_kraken_mpa_format
output_dir=/projects/tcga-data/biom_tables_kraken_unfiltered/${disease_type}/virus

for file in ${root_dir}/Viral_${disease_type}_*.output
do
	# Remove directory
	bn=${file##*/}
	fn=$(echo $bn | sed -e 's/.*_Kraken_\(.*\).output.*/\1/')
	# Remove extension
	bn=${bn%.*}
	# Remove last char (number)
	bn=${bn%?}
	echo "kraken-translate --mpa-format \
		  --db  ${db_fp} \
		  $file > ${output_dir}/${bn%.*}mpa_report_${fn}.txt" | qsub -l nodes=1:ppn=1 -l pmem=10Gb -l mem=10Gb -l walltime=48:00:00 -N kraken_translate_${bn%.*}mpa_report_${fn} -m abe -M jenya.kopylov@gmail.com
done