#!/usr/bin/env python3

"""
Author: Bekir Erguner
"""

import argparse, os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scan the CODEX2 result files per chromosome and extract the detected CNVs for each sample")
    parser.add_argument('--project_folder', '-p',
                        dest='project_folder',
                        help='Project output folder',
                        type=str)
    parser.add_argument('--chromosome_prefix', '-c',
                        dest='cp',
                        help='Set this as chr for hg38',
                        default='',
                        type=str)
    args = parser.parse_args()

    # Open and parse the file containing the list of clusters and the sample names
    sample_clusters = {}
    cluster_file = os.path.join(args.project_folder, "sample_clusters.tsv")
    with open(cluster_file, 'r') as cf:
        for line in cf.readlines():
            line = line.rstrip('\r').rstrip('\n')
            cluster = line.split('\t')[0]
            samples = line.split('\t')[1].split(',')
            for sample in samples:
                sample_clusters[sample] = cluster

    # Create the results folder f it doesn't exist already
    if not os.path.exists(os.path.join(args.project_folder, 'codex_results', 'sample_results')):
        os.mkdir(os.path.join(args.project_folder, 'codex_results', 'sample_results'))

    # Do the magic for all the samples found in the list
    for sample_name in sample_clusters:
        name_length = len(sample_name)
        sample_cluster = sample_clusters[sample_name]
        sample_cnv_file = os.path.join(args.project_folder, 'codex_results', 'sample_results',
                                       '{}_CODEX2.tsv'.format(sample_name))

        # Do not overwrite if there is results from a previous run
        if os.path.exists(sample_cnv_file):
            continue

        # Aggregate the results from the autosomal chromosomes
        sample_cnv_lines = ['sample_name\tchr\tcnv\tst_bp\ted_bp\tlength_kb\tst_exon\ted_exon\traw_cov\tnorm_cov\tcopy_no\tlratio\tmBIC\n']
        for chromosome in range(1, 23):
            cluster_file = os.path.join(args.project_folder, 'codex_results',
                                        'cls{}_CODEX2_chr{}{}.tsv'.format(sample_cluster, args.cp, chromosome))
            try:
                with open(cluster_file, 'r') as my_file:
                    for line in my_file.readlines():
                        if line[:name_length] == sample_name:
                            sample_cnv_lines.append(line)
            except FileNotFoundError:
                print('Could not open {}, skipping'.format(cluster_file))
        with open(sample_cnv_file, 'w') as out_file:
            out_file.write(''.join(sample_cnv_lines))
