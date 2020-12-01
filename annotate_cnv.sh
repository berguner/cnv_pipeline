#!/bin/bash

PROJECT_FOLDER=$1
SAMPLE_NAME=$2
CODEX_TSV="${PROJECT_FOLDER}/codex_results/sample_results/${SAMPLE_NAME}_CODEX2.tsv"
CODEX_BED="${PROJECT_FOLDER}/temp/${SAMPLE_NAME}_CODEX2.bed"
EDP_TSV="${PROJECT_FOLDER}/exomedepth_results/${SAMPLE_NAME}_ExomeDepth.tsv"
EDP_BED="${PROJECT_FOLDER}/temp/${SAMPLE_NAME}_ExomeDepth.bed"
UNION_BED="${PROJECT_FOLDER}/temp/${SAMPLE_NAME}_union.bed"
ANNO_TSV="${PROJECT_FOLDER}/temp/${SAMPLE_NAME}_AnnotSV.tsv"
FINAL_TSV="${PROJECT_FOLDER}/annotated_results/${SAMPLE_NAME}_AnnotSV.tsv"

if [[ ! -d "${PROJECT_FOLDER}/temp" ]]; then
  mkdir "${PROJECT_FOLDER}/temp"
fi

if [[ ! -d "${PROJECT_FOLDER}/annotated_results" ]]; then
  mkdir "${PROJECT_FOLDER}/annotated_results"
fi

if [[ -z "${ANNOTSV}" ]]; then
  echo "\$ANNOTSV environment variable was not set, please setup AnnotSV first."
  exit -1
fi

if [[ -f "$CODEX_TSV" && -f "$EDP_TSV" ]]; then
awk -v bed=${CODEX_BED} '{
    if(NR == 1) {
    print "#chr\tstart\tend\tCODEX2_ID|type|length(Kbp)|raw_cov|norm_cov|copy_no|lratio|mBIC" > bed;
    } else {
    print $2"\t"$4"\t"$5"\t"$2":"$4"-"$5"|"$3"|"($5-$4)/1000"|"$9"|"$10"|"$11"|"$12"|"$13;
    }
}' ${CODEX_TSV} | sort -k1,1 -k2,2n -k3,3n >> ${CODEX_BED} ;
awk -v bed=${EDP_BED} '{
    if(NR == 1){
    print "#chr\tstart\tend\tExomeDepth_ID|type|length(Kbp)|reads_expected|reads_observed|reads_ratio|BF|num_exons|Conrad" > bed;
    } else {
    print $7"\t"$5"\t"$6"\t"$7":"$5"-"$6"|"substr($3,1,3)"|"($6-$5)/1000"|"$10"|"$11"|"$12"|"$9"|"$4"|"$13;
    }
}' ${EDP_TSV} | sort -k1,1 -k2,2n -k3,3n >> ${EDP_BED} ;
bedtools unionbedg -i ${CODEX_BED} -i ${EDP_BED} -filler NA | \
  awk '{if($4 != "NA" || $5 != "NA") print $0}' > ${UNION_BED} ;
${ANNOTSV}/bin/AnnotSV -SVinputFile ${UNION_BED} -outputFile ${ANNO_TSV} > ${FINAL_TSV}.log 2>&1;
awk 'BEGIN {FS = "\t"; OFS = "\t";} {
    if(NR == 1) {
    $6="CODEX2_pos\tCDX_type\tCDX_length(Kbp)\tCDX_raw_coverage\tCDX_normalized_coverage\tCDX_copy_number\tCDX_likelihood";
    $7="ExomeDepth_pos\tEDP_type\tEDP_length(Kbp)\tEDP_reads_expected\tEDP_reads_observed\tEDP_reads_ratio\tEDP_BayesFactor\tEDP_num_exons\tEDP_Common_CNV";
    print $0;
    } else {
        $5 = ($4 - $3) / 1000;
        if($6 == "NA") {
        $6="NA\tNA\tNA\tNA\tNA\tNA\tNA";
        } else {
        split($6, codex, "|");
        $6=codex[1]"\t"codex[2]"\t"codex[3]"\t"codex[4]"\t"codex[5]"\t"codex[6]"\t"codex[7];
        }
        if($7 == "NA") {
        $7="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        } else {
        split($7, exdepth, "|");
        $7=exdepth[1]"\t"exdepth[2]"\t"exdepth[3]"\t"exdepth[4]"\t"exdepth[5]"\t"exdepth[6]"\t"exdepth[7]"\t"exdepth[8]"\t"exdepth[9];
        }
        print $0;
    }
}' ${ANNO_TSV} > ${FINAL_TSV};
# Clean up
rm ${CODEX_BED} ${EDP_BED} ${ANNO_TSV} ${UNION_BED};

elif [[ -f "$CODEX_TSV" ]]; then
awk -v bed=${CODEX_BED} '{
    if(NR == 1) {
    print "#chr\tstart\tend\tCODEX2_ID|type|length(Kbp)|raw_cov|norm_cov|copy_no|lratio|mBIC" > bed;
    } else {
    print $2"\t"$4"\t"$5"\t"$2":"$4"-"$5"|"$3"|"($5-$4)/1000"|"$9"|"$10"|"$11"|"$12"|"$13;
    }
}' ${CODEX_TSV} | sort -k1,1 -k2,2n -k3,3n >> ${CODEX_BED} ;
${ANNOTSV}/bin/AnnotSV -SVinputFile ${CODEX_BED} -outputFile ${ANNO_TSV} > ${FINAL_TSV}.log 2>&1;
awk 'BEGIN {FS = "\t"; OFS = "\t";} {
    if(NR == 1) {
    $6="CODEX2_pos\tCDX_type\tCDX_length(Kbp)\tCDX_raw_coverage\tCDX_normalized_coverage\tCDX_copy_number\tCDX_likelihood";
    print $0;
    } else {
        $5 = ($4 - $3) / 1000;
        if($6 == "NA") {
        $6="NA\tNA\tNA\tNA\tNA\tNA\tNA";
        } else {
        split($6, codex, "|");
        $6=codex[1]"\t"codex[2]"\t"codex[3]"\t"codex[4]"\t"codex[5]"\t"codex[6]"\t"codex[7];
        }
        print $0;
    }
}' ${ANNO_TSV} > ${FINAL_TSV};
# Clean up
rm ${CODEX_BED} ${ANNO_TSV};

elif [[ -f "$EDP_TSV" ]]; then
awk -v bed=${EDP_BED} '{
    if(NR == 1){
    print "#chr\tstart\tend\tExomeDepth_ID|type|length(Kbp)|reads_expected|reads_observed|reads_ratio|BF|num_exons|Conrad" > bed;
    } else {
    print $7"\t"$5"\t"$6"\t"$7":"$5"-"$6"|"substr($3,1,3)"|"($6-$5)/1000"|"$10"|"$11"|"$12"|"$9"|"$4"|"$13;
    }
}' ${EDP_TSV} | sort -k1,1 -k2,2n -k3,3n >> ${EDP_BED} ;
${ANNOTSV}/bin/AnnotSV -SVinputFile ${EDP_BED} -outputFile ${ANNO_TSV} > ${FINAL_TSV}.log 2>&1;
awk 'BEGIN {FS = "\t"; OFS = "\t";} {
    if(NR == 1) {
    $6="ExomeDepth_pos\tEDP_type\tEDP_length(Kbp)\tEDP_reads_expected\tEDP_reads_observed\tEDP_reads_ratio\tEDP_BayesFactor\tEDP_num_exons\tEDP_Common_CNV";
    print $0;
    } else {
        $5 = ($4 - $3) / 1000;
        if($6 == "NA") {
        $6="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        } else {
        split($6, exdepth, "|");
        $6=exdepth[1]"\t"exdepth[2]"\t"exdepth[3]"\t"exdepth[4]"\t"exdepth[5]"\t"exdepth[6]"\t"exdepth[7]"\t"exdepth[8]"\t"exdepth[9];
        }
        print $0;
    }
}' ${ANNO_TSV} > ${FINAL_TSV};
# Clean up
rm ${EDP_BED} ${ANNO_TSV};

fi