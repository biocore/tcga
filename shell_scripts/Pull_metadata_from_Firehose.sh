#!/bin/sh

# Adapted from https://github.com/cczysz/tcga/blob/master/mirna_dl.sh

# Example FTP dir
# http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/

# Example expression filename
# gdac.broadinstitute.org_BRCA.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0.tar.gz

# Example clinical phenotypes link
# http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz

diseasetypes='ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM'

for i in $(echo $diseasetypes); do
	# explink=http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$i/20160128/gdac.broadinstitute.org_$i.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0.tar.gz

	# outf=$i.tar.gz;
	# curl -L $explink -o ~/tmp/$i.tar.gz;
	# mkdir -p mirna_expression/$i;
	# tar -xvf ~/tmp/$i.tar.gz -C mirna_expression/$i/ --strip-components=1;
	# rm ~/tmp/*.tar.gz;

	phenlink=http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$i/20160128/gdac.broadinstitute.org_$i.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz
	outf=$i.tar.gz;
	curl -L $phenlink -o ./$outf;
	mkdir -p ./$i;
	tar -xvf ./$i.tar.gz -C ./$i --strip-components=1;
	cd ./$i
	for file in *.txt; do mv "$file" "${file/.txt/_$i.txt}"; done
	cd ..
done;

rm *.tar.gz
