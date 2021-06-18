#!/usr/bin/env python3

"""
Author: Bekir Erguner
"""

import argparse
import csv
import datetime
import gzip
import multiprocessing
import os
import subprocess
import sys
import threading
import time
from collections import OrderedDict

import hdbscan
import pandas
import plotly.express as px
from plotly.offline import plot
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler


def distance_to_cluster(index, cluster_indices, df):
    """
    This function calculates the distance of the outlier sample with the given index to the cluster with the given
    cluster_indices
    :param index: Index of the outlier sample
    :param cluster_indices: Indices of the cluster of interest
    :param df: pandas DataFrame containing the t-SNE components
    :return: Average distance between outlier sample and the samples in the cluster of interest
    """
    mysum = 0.0
    for i in cluster_indices:
        mysum += abs(df.loc[[index], ['component 1']].values[0] - df.loc[[i], ['component 1']].values[0])
        mysum += abs(df.loc[[index], ['component 2']].values[0] - df.loc[[i], ['component 2']].values[0])
    mysum /= len(cluster_indices)
    return mysum


def minimum_distance_label(distances):
    """
    This function finds the label of closest cluster
    :param distances: Dictionary of the label: distance information
    :return: key (label) of the closest cluster
    """
    min_value = -1
    min_key = ''
    for key in distances:
        if min_value == -1 or min_value >= distances[key]:
            min_value = distances[key]
            min_key = key
    return min_key


# Print iterations progress
def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clustering exome samples based on their read counts similarities")
    parser.add_argument('--project_folder', '-p',
                        dest='project_folder',
                        help='Project output folder',
                        type=str)
    parser.add_argument('--minimum_mean_read_count', '-m',
                        dest='min_mean_rc',
                        help='Minimum accepted mean read count per sample',
                        type=int)
    parser.add_argument('--min_cluster_size', '-c',
                        dest='min_cluster_size',
                        help='Minimum accepted cluster size',
                        type=int)
    args = parser.parse_args()

    rc_folder = os.path.join(args.project_folder, "read_counts")
    rc_files = [f for f in os.listdir(rc_folder) if f.find('_exome_RC.csv.gz') > 0]
    rc_dict = {}
    failed_samples = {}
    counter = 0
    print('Now retrieving read counts of the first 10k regions/exons for clustering the samples')
    for rc_file in rc_files:
        counter += 1
        # print_progress_bar(counter, len(rc_files))
        print('{}/{} read counts were retrieved'.format(counter, len(rc_files)), end='\r')
        sample_name = rc_file.replace('_exome_RC.csv.gz', '')
        rc_file = os.path.join(args.project_folder, 'read_counts', rc_file)
        rc = []
        try:
            rc_csv = gzip.open(rc_file)
            for line in rc_csv.readlines(500000)[1:10001]:
                rc.append(int(line.decode("utf-8").rstrip('\n').split(',')[1]))
        except:
            continue
        sum_rc = sum(rc)
        mean_rc = sum_rc / len(rc)
        if mean_rc > args.min_mean_rc:
            sum_rc_per_kb = sum_rc / 1000
            tmp_rc = []
            for c in rc:
                tmp_rc.append(c / sum_rc_per_kb)
            rc_dict[sample_name] = tmp_rc
        else:
            failed_samples[sample_name] = mean_rc
    if len(failed_samples) > 0:
        print('\n{} out of {} samples didn\'t have adequate coverage for CNV analysis:\n{}'.format(len(failed_samples),
                                                                                                   len(rc_files),
                                                                                                   '\n'.join(failed_samples.keys())))

    print('Now performing the clustering')
    samples_list = list(rc_dict.keys())
    cohort_df = pandas.DataFrame.from_dict(rc_dict)
    # Remove 0 coverage regions from the data frame
    cohort_df = cohort_df[(cohort_df.T != 0).any()]
    x = cohort_df.loc[:, samples_list].transpose().values
    x = StandardScaler().fit_transform(x)

    # Run tSNE
    time_start = time.time()
    RS = 123
    exome_tsne = TSNE(random_state=RS, n_components=2, n_iter=1000).fit_transform(
        PCA(n_components=50).fit_transform(x))
    print('t-SNE done! Time elapsed: {} seconds'.format(time.time() - time_start))

    # Cluster tSNE results
    clusterer = hdbscan.HDBSCAN(min_cluster_size=args.min_cluster_size)
    clusterer.fit(exome_tsne)

    tsneDf = pandas.DataFrame(data=exome_tsne, columns=['component 1', 'component 2'])

    # Count the number of samples for each cluster
    label_count_dict = {}
    for l in clusterer.labels_:
        if l not in label_count_dict:
            label_count_dict[l] = {'count': 1}
        else:
            label_count_dict[l] = {'count': label_count_dict[l]['count'] + 1}

    # Create a dictionary of samples and clusters they belong
    sample_labels = OrderedDict()
    for idx in range(len(clusterer.labels_)):
        sample_labels[samples_list[idx]] = {'label': clusterer.labels_[idx]}
        sample_labels[samples_list[idx]]['likely'] = clusterer.labels_[idx]
    print('Not clustered samples: {}'.format(label_count_dict[-1]['count']))
    label_df = pandas.DataFrame.from_dict(label_count_dict).transpose()

    # Assign labels/clusters to the non-clustered samples
    for i in range(len(clusterer.labels_)):
        if clusterer.labels_[i] == -1:
            print('Finding the closest cluster for sample {}'.format(samples_list[i]))
            distances = {}
            for label in label_count_dict:
                if label == -1:
                    continue
                indices = []
                for k in range(len(clusterer.labels_)):
                    if clusterer.labels_[k] == label:
                        indices.append(k)
                distances[label] = distance_to_cluster(i, indices, tsneDf)
            clusterer.labels_[i] = minimum_distance_label(distances)
            sample_labels[samples_list[i]]['likely'] = clusterer.labels_[i]

    # Plot the clustering results and save it as an html file
    tsneDf['labels'] = [str(l) for l in clusterer.labels_]
    tsneDf['samples'] = samples_list
    fig = px.scatter(tsneDf,
                     x='component 1',
                     y='component 2',
                     color='labels',
                     hover_data=['samples'])
    fig.update_layout(showlegend=False)
    now = datetime.datetime.now()
    plot_file = os.path.join(args.project_folder, now.strftime("%Y_%m_%d_%H.%M.%S_CNV_t_SNE.html"))
    plot(fig, filename=plot_file, auto_open=False)

    # Get the unique labels of clusters with new samples
    unique_labels = {}
    for l in range(len(clusterer.labels_)):
        new_sample_exists = False
        for sample_name in rc_dict:
            sample_result_file = os.path.join(args.project_folder,
                                              'annotated_results',
                                              '{}_AnnotSV.tsv'.format(sample_name))
            if not os.path.exists(sample_result_file) and sample_name in samples_list:
                if sample_labels[sample_name]['label'] == int(clusterer.labels_[l]) \
                        or sample_labels[sample_name]['likely'] == int(clusterer.labels_[l]):
                    new_sample_exists = True
        if str(clusterer.labels_[l]) not in unique_labels and new_sample_exists:
            unique_labels[str(clusterer.labels_[l])] = True

    output_clusters = {}
    for sample_name in sample_labels:
        sample_label = str(sample_labels[sample_name]['label'])
        if sample_label in unique_labels and unique_labels[sample_label]:
            if str(sample_label) in output_clusters:
                output_clusters[sample_label].append(sample_name)
            else:
                output_clusters[sample_label] = [sample_name]

    output_file = os.path.join(args.project_folder, "sample_clusters.tsv")
    with open(output_file, 'w') as output:
        for label in output_clusters:
            output.write("{}\t{}\n".format(label, ','.join(output_clusters[label])))
