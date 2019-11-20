#!/usr/bin/env Rscript

library(CODEX2)
library(argparser)


args <- arg_parser(description='Calculate sample exome coverage for CODEX')
args <- add_argument(args,arg='--bed', default='/scratch/lab_bsf/scratch/exome_panel_comparison/gnomad/exome_cnv_calling_regions.bed',
                     help='The path of the target regions bed file')
args <- add_argument(args,arg='--bam', short='-i',
                     help='Path of the input bam file')
args <- add_argument(args,arg='--project_folder', short='-p',
                     help='Project folder path')
args <- add_argument(args,arg='--sample_name', short='-s',
                     help='Sample name')
p <- parse_args(args)


bedFile <- p$bed
bamFile <- p$bam
sampname <- p$sample_name 

out_folder <- file.path(p$project_folder, 'read_counts')
if(!dir.exists(out_folder)) {
  dir.create(out_folder)
}

exomtarg <- read.table(bedFile, sep = "\t")
ref <- GRanges(seqnames = exomtarg[, 1], ranges = IRanges(start = exomtarg[, 2], end = exomtarg[, 3]))
ref <- sort(ref)
Y <- matrix(NA, nrow = length(ref), ncol = 1)
bamurl <- bamFile
what <- c("rname", "pos", "mapq", "qwidth")
flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                    isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
param <- ScanBamParam(what = what, flag = flag)
bam <- scanBam(bamurl, param = param)[[1]]
readlength.i = round(mean(bam$qwidth))
if (is.nan(readlength.i)) {
  flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                      isNotPassingQualityControls = FALSE)
  param <- ScanBamParam(what = what, flag = flag)
  bam <- scanBam(bamurl, param = param)[[1]]
  readlength.i = round(mean(bam$qwidth))
}

message("Getting coverage for sample ", sampname, 
        ": ", "read length ", readlength.i, ".", 
        sep = "")
bam.ref = GRanges(seqnames = bam$rname, ranges = IRanges(start = bam$pos, width = bam$qwidth))
bam.ref = bam.ref[bam$mapq >= 20]
Y[, 1] <- countOverlaps(ref, bam.ref)
rownames(Y) = paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
colnames(Y) = sampname
outz <- gzfile(file.path(p$project_folder, 'read_counts', paste(sampname, '_exome_RC.csv.gz', sep='')))
write.csv(Y, file = outz, quote = FALSE)

