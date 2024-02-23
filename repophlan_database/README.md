# TCGA Database Genomes

On June 14, 2016, a total of 71,782 microbial genomes were downloaded from NCBI using the [RePoPhlAn](https://github.com/SegataLab/repophlan) tool. Among these, 54,470 bacterial/archaeal genomes with a normalized quality score of at least 0.8 (as determined by RePoPhlAn's `screen.py` script, following the method outlined by [Land et al. (2014)](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/1944-3277-9-20), and normalization performed using [standarization](https://github.com/tanaes/script_bin/blob/master/filter_repophlan.py) since four scores are output per genome) were identified. Additionally, all 5,493 viral genomes, which included those with a `.fna` file extension and required no further filtering, were used to construct the Kraken 1 database. This database was created following the instructions provided on the [official Kraken website](https://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases). 

The genome IDs, including accession numbers, taxonomy IDs, NCBI FTP links, and more, for the 54,470 bacterial/archaeal genomes can be found in the file "repophlan_microbes_wscores_k__Bacteria_k__Archaea_screen.txt". Similarly, the genome IDs for the 5,493 viral genomes can be accessed in the file "repophlan_microbes_wscores_k__Viruses_k__Viroids_all.txt".


