#!/usr/bin/env python3

"""
Author: Bekir Erguner
"""

import os, sys
import pandas
import numpy as np
import subprocess
import multiprocessing
import time
import datetime
import threading
import gzip
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import hdbscan
import argparse
import yaml
import csv
import threading
import plotly.express as px
from plotly.offline import plot

class sacct_watcher(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.sacct_dict = {}
        self._alive = True
        self.update_count = 0

    def run(self):
        print('sacct_watcher is running')
        while self._alive:
            self.update_count += 1
            self.update_sacct()
            time.sleep(10)

    def update_sacct(self):
        cmd = ['sacct', '-p']
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        lines = result.stdout.decode("utf-8").split('\n')
        field_names = lines[0].rstrip('|').split('|')
        for l in range(1, len(lines) - 1):
            line = lines[l].rstrip('|').split('|')
            self.sacct_dict[line[0]] = {}
            for i in range(len(field_names)):
                self.sacct_dict[line[0]][field_names[i]] = line[i]

    def get_job_state(self, job_id):
        if job_id in self.sacct_dict:
            return self.sacct_dict[job_id]['State']
        else:
            return 'JOB NOT FOUND'

    def stop(self):
        print('sacct_watcher is stopped')
        self._alive = False


def sbatch_read_counting(new_samples, sample_annotations, config):
    the_watcher = sacct_watcher()
    the_watcher.start()
    coverage_script = os.path.join(config['pipeline_folder'],
                                   'codex_exome_coverage.R')
    bed_file = config['bed_file']
    project_folder = config['project_folder']
    coverage_jobs = {}
    for sample in new_samples:
        cmd = ['sbatch', f'--job-name={sample}_exome_RC', '--mem=8000', '--cpus=1', '--partition=shortq',
               '--error={}/read_counts/{}_exome_RC.csv.log'.format(project_folder, sample),
               '--output={}/read_counts/{}_exome_RC.csv.log'.format(project_folder, sample),
               coverage_script,
               '--bed', bed_file,
               '--bam', sample_annotations[sample]['bam_file'],
               '--project_folder', project_folder,
               '--sample_name', sample]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        job_id = result.stdout.decode("utf-8").rstrip('\n').split(' ')[3]
        coverage_jobs[job_id] = ''

    complete_count = 0
    fail_count = 0
    while complete_count + fail_count < len(coverage_jobs):
        complete_count = 0
        fail_count = 0
        for j in coverage_jobs:
            coverage_jobs[j] = the_watcher.get_job_state(j)
            if coverage_jobs[j] == 'COMPLETED':
                complete_count += 1
            elif coverage_jobs[j] == 'FAILED' or coverage_jobs[j] == 'CANCELLED':
                fail_count += 1
        #print(str(coverage_jobs))
        time.sleep(60)
    print('{} out of {} read count jobs were completed and {} of them had failed'.format(complete_count,
                                                                                         len(coverage_jobs),
                                                                                         fail_count))
    the_watcher.stop()


def local_read_counting(rc_args):
    new_sample = rc_args['new_sample']
    bam_file = rc_args['bam_file']
    config = rc_args['config']
    coverage_script = os.path.join(config['pipeline_folder'],
                                   'codex_exome_coverage.R')
    bed_file = config['bed_file']
    project_folder = config['project_folder']
    logfile = open('{}/read_counts/{}_exome_RC.csv.log'.format(project_folder, new_sample), 'w')
    cmd = [coverage_script,
           '--bed', bed_file,
           '--bam', bam_file,
           '--project_folder', project_folder,
           '--sample_name', new_sample]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in proc.stdout:
        logfile.write(line.decode('utf-8'))
    proc.wait()
    logfile.close()


def distance_to_cluster(index, cluster_indices, df):
    mysum = 0.0
    for i in cluster_indices:
        mysum += abs(df.loc[[index], ['component 1']].values[0] - df.loc[[i], ['component 1']].values[0])
        mysum += abs(df.loc[[index], ['component 2']].values[0] - df.loc[[i], ['component 2']].values[0])
    mysum /= len(cluster_indices)
    return mysum


def minimum_distance_label(distances):
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
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CNV Calling Pipeline")
    parser.add_argument('--config', '-c',
                        dest='config',
                        help='Project specific .yaml config file',
                        type=str)
    parser.add_argument('--engine', '-e',
                        dest='engine',
                        help='slurm or local',
                        default='slurm',
                        type=str)
    parser.add_argument('--threads', '-t',
                        dest='threads',
                        help='Number of threads to be used in case running locally. ~4GB of ram per thread is needed',
                        default=2,
                        type=int)
    parser.add_argument('--min-cluster-size', '-z',
                        dest='min_cluster_size',
                        help='Minimum number of samples for HDBSCAN clustering',
                        default=12,
                        type=int)
    parser.add_argument('--min-mean-rc', '-m',
                        dest='min_mean_rc',
                        help='Minimum accepted average read count per region',
                        default=20,
                        type=int)
    args = parser.parse_args()

    config = {}
    with open(args.config, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exception:
            sys.stderr.write(str(exception))
            exit(1)

    project_folder = config['project_folder']
    if not os.path.exists(project_folder):
        os.mkdir(project_folder)

    rc_folder = os.path.join(project_folder, 'read_counts')
    if not os.path.exists(rc_folder):
        os.mkdir(rc_folder)

    codex_folder = os.path.join(project_folder, 'codex_results')
    if not os.path.exists(codex_folder):
        os.mkdir(codex_folder)

    exomedepth_folder = os.path.join(project_folder, 'exomedepth_results')
    if not os.path.exists(exomedepth_folder):
        os.mkdir(exomedepth_folder)

    sample_annotations = {}
    with open(config['sample_sheet']) as csv_file:
        try:
            reader = csv.DictReader(csv_file)
            for row in reader:
                sample_annotations[row['sample_name']] = row
        except csv.Error as err:
            sys.stderr.write('Unable to read sample sheet\n')
            sys.stderr.write(str(err))
            exit(1)

    # Check which samples have read counts and which are new
    # and generate read counts for the new samples
    new_samples = []
    for sample_name in sample_annotations:
        bam_file_check = False
        sample_bam_file = sample_annotations[sample_name]['bam_file']
        bai1 = sample_bam_file[:-4] + '.bai'
        bai2 = sample_bam_file[:-4] + '.bam.bai'
        if os.path.exists(sample_bam_file):
            if os.path.exists(bai1) or os.path.exists(bai2):
                bam_file_check = True
        rc = os.path.join(rc_folder, '{}_exome_RC.csv.gz'.format(sample_name))
        if not os.path.exists(rc) and bam_file_check:
            # Check if the bam file is still being written or not
            mod_date = os.path.getmtime(sample_bam_file)
            d = datetime.datetime.now() - datetime.datetime.fromtimestamp(mod_date)
            if d.seconds > 3600:
                new_samples.append(sample_name)

    print('There are {} sample(s) without read counts'.format(len(new_samples)))
    if len(new_samples) > 0:
        print('Now counting the reads for the new sample(s)')
        if args.engine == 'slurm':
            sbatch_read_counting(new_samples, sample_annotations, config)
        elif args.engine == 'local':
            list_of_arg_dicts = []
            for new_sample in new_samples:
                list_of_arg_dicts.append({'new_sample': new_sample,
                                          'bam_file': sample_annotations[new_sample]['bam_file'],
                                          'config': config})
            rc_pool = multiprocessing.Pool(args.threads)
            rc_pool.map(local_read_counting, list_of_arg_dicts)

    rc_files = [f for f in os.listdir(rc_folder) if f.find('_exome_RC.csv.gz') > 0]
    rc_dict = {}
    failed_samples = {}
    counter = 0
    print('Now retrieving read counts of first 10k regions/exons for clustering the samples')
    for rc_file in rc_files:
        counter += 1
        print_progress_bar(counter, len(rc_files))
        sample_name = rc_file.replace('_exome_RC.csv.gz', '')
        rc_file = os.path.join(project_folder, 'read_counts', rc_file)
        rc = []
        with gzip.open(rc_file) as rc_csv:
            for line in rc_csv.readlines(500000)[1:10001]:
                rc.append(int(line.decode("utf-8").rstrip('\n').split(',')[1]))
        sum_rc = sum(rc)
        mean_rc = sum_rc / len(rc)
        if sample_name in sample_annotations:
            sample_annotations[sample_name]['mean_read_count'] = mean_rc
        if mean_rc > args.min_mean_rc:
            sum_rc_per_kb = sum_rc / 1000
            tmp_rc = []
            for c in rc:
                tmp_rc.append(c / sum_rc_per_kb)
            rc_dict[sample_name] = tmp_rc
        else:
            failed_samples[sample_name] = mean_rc
    if len(failed_samples) > 0:
        print('\n{} out of {} samples did not have adequate coverage for CNV analysis:\n{}'.format(len(failed_samples),
                                                                                                   len(rc_files),
                                                                                                   '\n'.join(failed_samples.keys())))
    cohort_df = pandas.DataFrame.from_dict(rc_dict)

    samples_list = list(rc_dict.keys())
    # Remove not covered regions from the data frame
    cohort_df = cohort_df[(cohort_df.T != 0).any()]
    x = cohort_df.loc[:, samples_list].transpose().values
    x = StandardScaler().fit_transform(x)

    # Run tSNE
    time_start = time.time()
    RS = 123
    exome_tsne = TSNE(random_state=RS, n_components=2, n_iter=1000).fit_transform(PCA(n_components=50).fit_transform(x))
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
    sample_labels = {}
    for idx in range(len(clusterer.labels_)):
        sample_labels[samples_list[idx]] = {'label': clusterer.labels_[idx]}
    print('Not clustered samples: {}'.format(label_count_dict[-1]['count']))
    label_df = pandas.DataFrame.from_dict(label_count_dict).transpose()

    # Assign labels/clusters to the non-clustered samples
    likely_labels = {}
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
            likely_labels[i] = minimum_distance_label(distances)
            clusterer.labels_[i] = likely_labels[i]
            sample_labels[samples_list[i]]['likely'] = likely_labels[i]

    # Add clustering labels to sample annotations for saving them later
    for sample_name in sample_annotations:
        if sample_name in sample_labels:
            sample_annotations[sample_name]['cluster_label'] = sample_labels[sample_name]['label']
            if 'likely' in sample_labels[sample_name]:
                sample_annotations[sample_name]['closest_cluster_label'] = sample_labels[sample_name]['likely']
            else:
                sample_annotations[sample_name]['closest_cluster_label'] = 'NA'
        else:
            sample_annotations[sample_name]['cluster_label'] = 'NA'
            sample_annotations[sample_name]['closest_cluster_label'] = 'NA'

    # Plot the clustering results and save it as an html file
    tsneDf['labels'] = clusterer.labels_
    tsneDf['samples'] = samples_list
    fig = px.scatter(tsneDf,
                     x='component 1',
                     y='component 2',
                     color='labels',
                     hover_data=['samples'])
    fig.update_layout(showlegend=False)
    plot_file = os.path.join(project_folder, 't_SNE.html')
    plot(fig, filename=plot_file, auto_open=False)
