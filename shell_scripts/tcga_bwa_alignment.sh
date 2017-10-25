#!/bin/bash
#
# script: Align TCGA Kraken-classified reads vs. 60K RepoPhlAn bacteria + virus database.
#

repophlan_db=/databases/repophlan/05_27_2016/bwa-index/bacteria-viruses-dust.fasta
tcga_clean_dir=/projects/tcga-data/human_filter
#disease_of_interest=("Bacterial_Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_Classified" \
#					 "Viral_Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_Classified" \
#					 "Bacterial_Stomach_Adenocarcinoma_Classified" \
#					 "Viral_Stomach_Adenocarcinoma_Classified")
disease_of_interest=("Bacterial_Ovarian_Serous_Cystadenocarcinoma_Classified" \
					 "Viral_Ovarian_Serous_Cystadenocarcinoma_Classified")
output_dp=/projects/tcga-data/bwa_alignments
threads=15

# Loop through Viral and Bacterial results for disease of interest
for disease in "${disease_of_interest[@]}"
do
	# Loop through files for each disease of interest
	for file in ${tcga_clean_dir}/${disease}_*.fasta.bz2
	do
		# 1) Unzip (human DNA filtered) file to working directory, keep original
		# 2) Run BWA
		# 3) Remove unzipped working file (original stored)
		bn=${file##*/};
		echo "bzip2 -dkc $file > ${output_dp}/${bn%.*}; \
			  bwa mem -t ${threads} -o ${output_dp}/${bn%.*}.sam ${repophlan_db} ${output_dp}/${bn%.*} ; \
			  rm ${output_dp}/${bn%.*}" | qsub -m abe -M jenya.kopylov@gmail.com -l nodes=1:ppn=${threads} -l pmem=56gb -l mem=850gb -N ${bn%.*} -l walltime=48:00:00
		sleep 1
	done
done
