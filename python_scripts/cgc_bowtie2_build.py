#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova, Jad Kanbar
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Build Bowtie2 database on all reference genomes in Kraken report.
"""

import click
from os import listdir
from os.path import join, dirname


def load_kraken_mpa_report(kraken_mpa_report_fp,
                           read_per_taxa,
                           taxonomic_rank):
    """Absolute abundance of number of reads matching a defined taxa level.

    Parameters
    ----------
    kraken_mpa_report_fp: tuple
        filepath(s) to output of "kraken mpa report"
    taxa_levels: list
        list of two elements that includes the taxonomic rank at which
        to generate summary and rank below to split by
    read_per_taxa: int
        integer of number of minimum number of reads to keep per taxa
    taxonomic_rank: str
        taxonomic rank at which to generate summary

    Returns
    -------
    taxonomies: set
        set of taxonomies from kraken_mpa_report_fp
    """
    taxa_levels = {"domain": "d__", "phylum": "|p__", "class": "|c__",
                   "order": "|o__", "family": "|f__", "genus": "|g__",
                   "species": "|s__", "strain": "|t__"}

    taxa_levels_idx = {"d__": 0, "|p__": 1, "|c__": 2, "|o__": 3, "|f__": 4,
                       "|g__": 5, "|s__": 6, "|t__": 7, "7": "|t__",
                       "6": "|s__", "5": "|g__", "4": "|f__", "3": "|o__",
                       "2": "|c__", "1": "|p__", "0": "d__"}

    taxa_levels_str = taxa_levels[taxonomic_rank]
    taxa_levels_idx_int = taxa_levels_idx[taxa_levels_str]

    if taxa_levels_idx_int < 7:
        split_on_level = taxa_levels_idx[str(taxa_levels_idx_int + 1)]
    else:
        split_on_level = '\t'

    taxonomic_abundances = {}
    for report_fp in kraken_mpa_report_fp:
        with open(report_fp) as report_f:
            for line in report_f:
                label, taxonomy = line.strip().split('\t')
                # read assigned taxonomy at least up to desired level
                if taxa_levels_str in taxonomy:
                    # keep taxonomy string up to specified level
                    taxonomy_parse = taxonomy.split(split_on_level)[0]
                    taxonomy_parse = taxonomy_parse.replace('d__', 'k__')
                    # format taxonomy strings to RepoPhlAn's style
                    if '.' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace('.', '')
                    if '(' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace('(', '_')
                    if ')' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace(')', '')
                    if '-' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace('-', '_')
                    if '[' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace('[', '')
                    if ']' in taxonomy_parse:
                        taxonomy_parse = taxonomy_parse.replace(']', '')
                    if taxonomy_parse not in taxonomic_abundances:
                        taxonomic_abundances[taxonomy_parse] = 1
                    else:
                        taxonomic_abundances[taxonomy_parse] += 1

    taxonomic_set = set([k for k, v in taxonomic_abundances.items()
                         if v > read_per_taxa])

    return taxonomic_set, split_on_level


def parse_repophlan_genome_ids(repophlan_genome_id_taxonomy_fp,
                               split_on_level):

    """Parse genome IDs and corresponding taxonomies based on taxa level.

    Parameters
    ----------
    repophlan_genome_id_taxonomy_fp: str
        filepath to output of repophlan file genome IDs and associated taxa
    split_on_level: str
        string that determines the level to split the taxonomy string

    Returns
    -------
    repophlan_taxa_and_genomes: dict
        Keys are taxonomies and values are lists of genomes associated
        to taxonomy
    """
    repophlan_taxa_and_genomes = {}
    # assign genome IDs to taxonomies truncated at desired level
    with open(repophlan_genome_id_taxonomy_fp) as genome_id_taxonomy_f:
        for line in genome_id_taxonomy_f:
            genome_id, path, taxonomy =\
                line.strip().split('\t')
            run = 0
            taxonomy = taxonomy.split(split_on_level)[0]
            # remove "noname" filler from RepoPhlAn's taxonomy
            if "_noname" in taxonomy:
                tax = taxonomy.split('|')
                edit_tax = [x for x in tax if '_noname' not in x]
                taxonomy = '|'.join(edit_tax)
            if taxonomy not in repophlan_taxa_and_genomes:
                repophlan_taxa_and_genomes[taxonomy] = [genome_id]
            else:
                repophlan_taxa_and_genomes[taxonomy].append(genome_id)

    return repophlan_taxa_and_genomes


def collect_candidate_repophlan_genomes(taxonomic_set,
                                        repophlan_taxa_and_genomes,
                                        genome_dir):
    """Assign genomes based on Kraken taxonomy breakdown.

    Parameters
    ----------
    taxonomic_set: set
        set of taxonomies from kraken_mpa_report_fp
    repophlan_taxa_and_genomes: dict
        Keys are taxonomies and values are lists of genomes associated
        with taxonomy
    genome_dir: str
        genome directory path

    Returns
    -------
    genomes_to_align_to: list
        list of genome IDs
    genome_fps: list
        bt2 input index files
    """
    genomes_to_align_to = []
    for taxa in taxonomic_set:
        if taxa in repophlan_taxa_and_genomes:
            genomes_to_align_to.extend(repophlan_taxa_and_genomes[taxa])
        else:
            print("[WARNING] %s not found in RepoPhlAn taxonomies" % taxa)
    # Generate bt2 command
    genome_fps = []
    candidate_genomes = []
    for filename in listdir(genome_dir):
        if filename.endswith('.fna'):
            genome_id = filename.split('.fna')[0]
            if genome_id in genomes_to_align_to:
                genome_fps.append(join(genome_dir, filename))
                candidate_genomes.append(genome_id)
    if set(candidate_genomes) != set(genomes_to_align_to):
        raise ValueError('Missing genomes: %s' % (
            set(candidate_genomes) - set(genomes_to_align_to)))
    return genomes_to_align_to, genome_fps


def output_candidate_repophlan_genomes(genomes_to_align_to,
                                       genome_fps,
                                       output_fp,
                                       bt2_index_base,
                                       threads):
    """
    Parameters
    ----------
    genomes_to_align_to: list
        list of genome IDs
    genome_fps_str: str
        bt2 input index files str
    output_fp: str
        string with output filename
    bt2_index_base: str
        Bowtie2 index dir/basename
    """
    # write genomes file
    output_dp = dirname(output_fp)
    output_genomes_fp = join(output_dp, "reference_db.fasta")
    with open(output_genomes_fp, 'w') as f:
        for fname in genomes_fps:
            with open(fname) as infile:
                for line in infile:
                    f.write(line)
    genomes_fps_str = ",".join(genome_fps)
    print("Total genomes: %s\n" % len(genomes_to_align_to))
    #with open(output_fp, 'w+') as output_f:
    #    output_f.write("echo '")
    #    output_f.write("bowtie2-build ")
    #    output_f.write(genomes_fps_str)
    #    output_f.write(" %s" % bt2_index_base)
    #    output_f.write(" --threads %s" % threads)
    #    output_f.write(" ' | qsub -l nodes=1:ppn=32 -q highmem -l walltime=72:00:00 -N bowtie2_build_%s_genomes" % len(genomes_to_align_to))


@click.command()
@click.option('--kraken-mpa-report-fp', required=True, multiple=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to Kraken report')
@click.option('--repophlan-genome-id-taxonomy-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to RepoPhlAn genome ID and taxonomy list')
@click.option('--genome-dir', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Genome directory path')
@click.option('--taxonomic-rank', type=click.Choice(['genus', 'species',
                                                     'family', 'order',
                                                     'class', 'phylum',
                                                     'domain']),
              required=False, default=['genus'], show_default=True,
              help="Taxonomic rank at which to generate summary")
@click.option('--read-per-taxa', required=False, type=int,
              default=10, show_default=True,
              help="Minimum number of reads for each taxa")
@click.option('--threads', required=False, type=int, default=1, show_default=True,
              help="Bowtie2-build threads option")
@click.option('--output-fp', required=False,
              default='subset_repophlan_genomeID_taxonomy.good',
              show_default=True,
              help="Output filepath for genomes to use in reduced index")
@click.option('--bt2-index-base', required=False,
              default='rep_genomes.idx',
              help="Write Bowtie2 data to files with this dir/basename")
def main(kraken_mpa_report_fp,
         repophlan_genome_id_taxonomy_fp,
         genome_dir,
         taxonomic_rank,
         read_per_taxa,
         threads,
         output_fp,
         bt2_index_base):

    taxonomic_set, split_on_level =\
        load_kraken_mpa_report(
            kraken_mpa_report_fp=kraken_mpa_report_fp,
            read_per_taxa=read_per_taxa,
            taxonomic_rank=taxonomic_rank)

    repophlan_taxa_and_genomes = parse_repophlan_genome_ids(
        repophlan_genome_id_taxonomy_fp=repophlan_genome_id_taxonomy_fp,
        split_on_level=split_on_level)

    genomes_to_align_to, genome_fps =\
        collect_candidate_repophlan_genomes(
            taxonomic_set=taxonomic_set,
            repophlan_taxa_and_genomes=repophlan_taxa_and_genomes,
            genome_dir=genome_dir)

    output_candidate_repophlan_genomes(
        genomes_to_align_to=genomes_to_align_to,
        genome_fps=genome_fps,
        output_fp=output_fp,
        bt2_index_base=bt2_index_base,
        threads=threads)


if __name__ == "__main__":
    main()
