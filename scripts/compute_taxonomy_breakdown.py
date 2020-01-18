#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Compute number of unique taxonomies and their counts

import sys


def get_genome_paths(scores_average_fp, scores_repophlan_fp):
    """Return taxonomy strings and their counts.
    """
    genomes = []
    # Get genome IDs
    with open(scores_average_fp) as scores_average_f:
        next(scores_average_f)
        for line in scores_average_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            if genome_id not in genomes:
                genomes.append(genome_id)
            else:
                raise ValueError("Duplicate genome IDs %s" % genome_id)
    # Get taxonomies based on used genome IDs
    taxonomies = {}
    genomes_without_taxonomy = []
    with open(scores_repophlan_fp) as scores_repophlan_f:
        # header
        line = scores_repophlan_f.readline()
        line = line.strip().split('\t')
        tax_idx = line.index('taxonomy')
        for line in scores_repophlan_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            # only want tax_ids for genomes passing quality filter
            if genome_id in genomes:
                taxonomy = line[tax_idx]
                if (('k__Bacteria' not in taxonomy) and
                        ('k__Archaea' not in taxonomy) and
                        ('k__Viruses' not in taxonomy) and
                        ('k__Viroids' not in taxonomy)):
                    genomes_without_taxonomy.append(genome_id)
                if taxonomy not in taxonomies:
                    taxonomies[taxonomy] = 1
                else:
                    taxonomies[taxonomy] += 1
    return ([(k, taxonomies[k]) for k in sorted(taxonomies, key=taxonomies.get, reverse=True)],
            genomes_without_taxonomy)


def main():
    """Parse output of RepoPhlAn's repophlan_get_microbes.txt to obtain
       unique taxonomies and their counts.
    """
    # Output of RepoPhlAn's repophlan_get_microbes.py
    repophlan_scores_fp = sys.argv[1]
    # Output of Jon's score average script
    repophlan_scores_average_fp = sys.argv[2]

    taxonomies, genomes_without_taxonomy = get_genome_paths(
        repophlan_scores_average_fp, repophlan_scores_fp)
    if len(genomes_without_taxonomy) != 0:
        print("Genomes without taxonomy: %s" % len(genomes_without_taxonomy))

    for tup in taxonomies:
        print("%s\t%s" % (tup[0], tup[1]))


if __name__ == '__main__':
    main()
