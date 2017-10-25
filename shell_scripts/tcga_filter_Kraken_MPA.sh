#!/bin/bash
#
# script: Use BWA alignments to filter false-positive reads from Kraken MPA format.
#

# Project directory
project_dp=/projects/tcga-data

# BWA alignments
bwa_dp=/projects/tcga-data/bwa_alignments

# Kraken MPA format
kraken_dp=/projects/tcga-data/

# Diseases of interest
disease_of_interest=("Stomach_Adenocarcinoma" "Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma")

# Output directory
tcga_clean_dp=/projects/tcga-data/clean_Kraken_MPA

# Loop through each disease of interest
for disease in "${disease_of_interest[@]}"
do
	# Loop through MPA reports for each disease of interest for *Bacterial* results
	for mpa_report in ${project_dp}/${disease}/data/Bacterial_${disease}_Classified_mpa_report_*.txt
	do
		# MPA report file name (no directory path)
		bn=${mpa_report##*/}
		# MPA report file number
		fn=$(echo $bn | sed -e 's/.*_mpa_report_\(.*\).txt.*/\1/')
		# SAM alignment file with same file number
		sam_fp=${bwa_dp}/Bacterial_${disease}_Classified_$fn.fasta.sam
		if [ -f ${sam_fp} ] then;
			# TODO: modify command to run filtering
			echo "run_filter.py --mpa-report ${mpa_report} \
				  --sam ${sam_fp} \
				  --output-fp ${tcga_clean_dp}/${bn%.*}_clean.txt" | qsub -m abe -M jenya.kopylov@gmail.com -l nodes=1:ppn=1 -l pmem=100gb -l mem=100gb -N ${bn%.*} -l walltime=48:00:00
			sleep 1
		fi
	done
done
