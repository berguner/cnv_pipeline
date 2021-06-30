#!/usr/bin/env bash

PROJECT_FOLDER=$1
BAM_FOLDER=$2
BAM_PREFIX=$3
BAM_SUFFIX=$4
COUNT_SCRIPT="~/src/cnv_pipeline/codex_exome_coverage.R"
BED_FILE="/research/lab_bsf/resources/intervals/exome_cnv_calling_regions_b37.bed"
TRHEADS=32

for i in ${BAM_FOLDER}/*${BAM_SUFFIX}; do \
  BN=$(basename $i);
  SN=${BN/${BAM_PREFIX}/};
  SN=${SN/${BAM_SUFFIX}/};
  echo -e "Rscript --vanilla  ${COUNT_SCRIPT} --bed ${BED_FILE} --project_folder ${PROJECT_FOLDER} --bam ${BAM_FOLDER}/${BN} --sample_name ${SN}";
  done | parallel -j ${TRHEADS} :::