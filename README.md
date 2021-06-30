# Germline CNV Detection Pipeline for Whole Exome data

## Introduction
This pipeline is designed to detect CNV in many exome samples prepared with various exome capture kits.
The pipeline includes a clustering step which attempts to separate samples into smaller groups which have similar coverage patterns.
Given that there are enough number (>12) of samples from each exome panel/kit, the clustering should reduce the batch effects and increase the sensitivity and specificity.
The pipeline currently supports GRCh37 genome version only.

## Software Requirements
If you can't use Docker/Singularity:
* Python >3.6 
    * hdbscan, scikit-learn, plotly, pandas, sci-py
* R >3.5
* [AnnotSV](https://lbgi.fr/AnnotSV/) for annotation

## How to run
* Clone this repository.
* Download a copy of [Cromwell](https://github.com/broadinstitute/cromwell)
* Run `Rscript install_CNV_prerequisites.R` to install required R packages.\
   or\
   pull the Docker container:\
   `docker pull berguner/bsf_cnv_pipeline:latest`\
   `singularity pull docker://berguner/bsf_cnv_pipeline:latest`
* Collect the aligned and duplicate marked .bam and .bai files in a folder. Symbolic links would also work.
* Update the `test/test.inputs.json` file according to your setup.
  * Create a sample annotation sheet (CSV) with a column named `Sample Name` containing the sample names. 
* Create a backend configuration file for Cromwell. You can modify the one of the provided backends in this repository.
* `cd` to the `project_folder` and run the pipeline:
```
java -Xmx4g -Dconfig.file=/path/to/backend.conf \
  -jar /path/to/cromwell-59.jar \
  run /path/to/cnv_pipeline/cnv_pipeline.wdl \
  --inputs /path/to/my.inputs.json
```

## Pipeline steps
1. Count reads on exonic regions.
2. Clustering samples based on coverage (read count) patterns.
3. Running CNV calling on each cluster of samples.
4. Merging results from [CODEX2](https://github.com/yuchaojiang/CODEX2) and [ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html).
5. Annotation of merged results using AnnotSV.
6. Aggregation of deleteions found in the cohort.