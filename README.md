# Germline CNV Detection Pipeline for Whole Exome data

## Introduction
This pipeline is designed to detect CNV in many exome samples prepared with various exome capture kits. The pipeline includes a clustering step
which attempts to separate samples into smaller groups which have similar coverage patterns. Given that there are enough number of samples from
each exome panel/kit, the clustering should reduce the batch effects and increase the accuracy.

## Software Requirements

* Python 3.6
    * hdbscan, scikit-learn, plotly, pandas, sci-py
* R 3.5 or greater
* [AnnotSV](https://lbgi.fr/AnnotSV/)
* Slurm (sbatch & sacct) - if running on a HPC cluster

## How to run
1. Clone this repository.
2. Run `Rscript install_CNV_prerequisites.R` to install required R packages.
2. Collect the aligned and duplicate marked .bam and .bai files in a folder. Symbolic links would also work.
3. Update the `data/config.yaml` file according to your setup.
4. Run `python3 run_cnv_pipeline.py -h` and follow the instructions

## Pipeline steps
1. Count reads on exonic regions.
2. Clustering samples based on coverage (read count) patterns.
3. Running CNV calling on each cluster of samples.
4. Merging results from [CODEX2](https://github.com/yuchaojiang/CODEX2) and [ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html).
5. Annotation of merged results using AnnotSV.