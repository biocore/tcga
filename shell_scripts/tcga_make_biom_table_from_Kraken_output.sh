#!/bin/bash

#-----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

root_dir=/projects/tcga-data/biom_tables_kraken_unfiltered/
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
scripts=/home/evko1434/tcga/scripts

for file in ${root_dir}/${disease_type}/virus/Viral_${disease_type}_Classified_Kraken_2mpa_report_*.txt
do
    echo $file
    bn=${file##*/}
    bn=${bn%.*}
    echo $bn
    echo "unset PYTHONPATH; \
          source activate parse_kraken_to_biom; \
          python $scripts/parse_kraken_to_biom.py \
                 --kraken-translate-report-fp $file \
                 --taxonomic-rank genus \
                 --biom-output-fp ${root_dir}/${disease_type}/virus/${bn}.biom" | qsub -l nodes=1:ppn=1 -l pmem=5Gb -l mem=5Gb -l walltime=48:00:00 -N build_biom_${bn}
done