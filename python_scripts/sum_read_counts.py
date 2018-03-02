#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
import os


def compute_microbial_read_counts(data_dp):
    """ Compute number of microbial classified reads and length interval.
    """
    total_bac_classified = 0
    total_vir_classified = 0
    total_unclassified = 0
    total_bac_files = 0
    total_vir_files = 0
    read_len_min = 100000
    read_len_max = 0
    read_len_avg = 0
    for fn in os.listdir(data_dp):
        if os.path.isfile(os.path.join(data_dp,fn)):
            # corrupt file
            if "Viral_Esophageal_Carcinoma_Classified_Kraken_1.output" in fn:
                continue
            elif "Bacterial_Uterine_Corpus_Endometrial_Carcinoma_Classified_Kraken_0.output" in fn:
                continue
            if "_Classified_Kraken" in fn:
                print("Processing %s" % fn)
                if "Bacterial_" in fn:
                    total_bac_files+=1
                    with open(os.path.join(data_dp,fn)) as f:
                        for line_n, line in enumerate(f):
                            line = line.strip().split()
                            if line[0] == 'C':
                                total_bac_classified+=1
                                read_len = int(line[3])
                                if read_len > read_len_max:
                                    read_len_max = read_len
                                elif read_len < read_len_min:
                                    read_len_min = read_len
                                read_len_avg += read_len
                            elif line[0] == 'U':
                                total_unclassified+=1
                            if line_n % 10000000 == 0:
                                print("%s ..." % line_n)
                elif "Viral_" in fn:
                    total_vir_files+=1
                    with open(os.path.join(data_dp,fn)) as f:
                        for line_n, line in enumerate(f):
                            line = line.strip().split()
                            if line[0] == 'C':
                                total_vir_classified+=1
                                read_len = int(line[3])
                                if read_len > read_len_max:
                                    read_len_max = read_len
                                elif read_len < read_len_min:
                                    read_len_min = read_len
                                read_len_avg += read_len
                            else:
                                print("Unknown letter code: %s" % line[0])
                            if line_n % 10000000 == 0:
                                print("%s ..." % line_n)
    total_unclassified = total_unclassified - total_vir_classified
    total_non_human_reads = total_bac_classified + total_unclassified
    read_len_avg = read_len_avg / (total_bac_classified + total_vir_classified)
    return (total_non_human_reads, total_bac_classified, total_vir_classified,
            total_unclassified, total_bac_files, total_vir_files, read_len_avg,
            read_len_min, read_len_max)


@click.command()
@click.option('--data-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with data files per disease')
def main(data_dp):
    tup = compute_microbial_read_counts(data_dp)
    print("Total_non_human_reads\ttotal_bac\ttotal_vir\ttotal_unclassified\ttotal_bac_files\ttotal_vir_files\tread_len_avg\tread_len_min\tread_len_max")
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % tup)


if __name__ == "__main__":
    main()