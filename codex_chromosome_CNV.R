#!/usr/bin/env Rscript

library(CODEX2)
library(argparser)
library(rtracklayer)

args <- arg_parser(description='Calculate sample exome coverage for CODEX')
args <- add_argument(args,arg='--bed', default='Y:/lab_bsf/scratch/exome_panel_comparison/gnomad/exome_cnv_calling_regions.bed',
                     help='The path of the target region bed file')
args <- add_argument(args,arg='--project_folder', short='-p',
                     default='Y:/lab_bsf/scratch/exome_panel_comparison/BSA_0005/',
                     help='Project folder path')
args <- add_argument(args,arg='--sample_names', short='-s',
                     default = '2012_4006_PB,2013_2597_PB,2013_2598_PB,2013_2599_PB,2013_2600_PB,2013_4593_TU,2013_4594_TU,IBD_392,IBD_484,PID_1113,PID_1126,PID_1224,PID_558,PID_587,PID_602,PID_612,PID_621,PID_622,PID_623,PID_681,PID_689,PID_694,PID_695,PID_697,PID_845,PID_889',
                     help='Comma separated sample names in the same order as bam files')
args <- add_argument(args,arg='--cluster', default='0',
                     help='Identifier for the batch of samples')
args <- add_argument(args,arg='--chromosome', short='-c', default='1',
                     help='Chromosome/contig name')
args <- add_argument(args,arg='--force', short='-f', default='y',
                     help='Set to n if you do not want to overwrite the existing results')
args <- add_argument(args,arg='--write_sample_output', default='n',
                     help='Set to y if you want the results per sample')
args <- add_argument(args, arg='--gtf',
                     default = 'Y:/lab_bsf/resources/genomes/GRCh37_e87/Homo_sapiens.GRCh37.87.gtf.gz',
                     help='The path of the Ensembl GTF file')
args <- add_argument(args, arg='--assembly',
                     default = 'b37',
                     help='The reference genome version; b37 or hg38')
p <- parse_args(args)

# Create output folder
if(!dir.exists(file.path(p$project_folder, 'codex_results'))){
  dir.create(file.path(p$project_folder, 'codex_results'))
}

# Generate the list of sample names
sample_names <- as.matrix(strsplit(p$sample_names,split=',')[[1]])
bedFile <- p$bed

# Generate the exome IRanges object
exomtarg <- read.table(bedFile, sep = "\t")
ref <- GRanges(seqnames = exomtarg[, 1], ranges = IRanges(start = exomtarg[, 2], end = exomtarg[, 3]))
ref <- sort(ref)

# Genome is b37/hg19 by default
if(p$assembly == 'b37'){
  genome = BSgenome.Hsapiens.UCSC.hg19
} else if(p$assembly == 'hg38'){
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome = BSgenome.Hsapiens.UCSC.hg38
} else{
  quit(status=1)
}

# Gather the read counts of the exome samples from previously generated csv files
Y <- matrix(NA, nrow = length(ref), ncol = length(sample_names))
rownames(Y) <- paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
colnames(Y) <- sample_names
for(i in 1:length(sample_names)) {
  fp <- file.path(p$project_folder, 'read_counts', paste(sample_names[i], '_exome_RC.csv.gz', sep=''))
  if(file.exists(fp)){
    message(paste('Reading coverage file of the sample:', sample_names[i]))
    tmp <- read.csv(file = gzfile(fp), col.names = cbind('region', sample_names[i]), check.names = FALSE, header = TRUE)
    Y[,i] <- unlist(tmp[,sample_names[i]])
  }
  else {
    message(paste('Coverage file for sample ', sample_names[i], ' was not found!'))
    quit(status=1)
  }
}
rm(tmp)

# Get the GC content and mapabilitiy values for the exonic ranges
gc <- getgc(ref, genome = genome)
mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))

# Perform QC for samples and regions
qcObj <- qc(Y, sample_names, ref, cov_thresh = c(20, 4000),
            length_thresh = c(20, 2000), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))

Y_qc <- qcObj$Y_qc
sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc
qcmat <- qcObj$qcmat
gc_qc <- ref_qc$gc

qc_path <- file.path(p$project_folder, 'codex_results', paste('codex_qc_cls_', p$cluster, '.tsv', sep=''))
if(!file.exists(qc_path)) {
  write.table(qcmat, file = qc_path, sep = '\t', quote = FALSE, row.names = FALSE)
}

# Estimate library size for the samples
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero, 1, function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)

# Grab the indices of the exons of the chromosomes of interest
chr.index <- which(seqnames(ref_qc) == p$chromosome)

# Perform normalization based on the mixture of poissons model
normObj.null <- normalize_null(Y_qc = Y_qc[chr.index,],
                               gc_qc = gc_qc[chr.index],
                               K = 1:5, N = N)
Yhat.null <- normObj.null$Yhat
AIC.null <- normObj.null$AIC
BIC.null <- normObj.null$BIC
RSS.null <- normObj.null$RSS

# Plot AIC, BIC and RSS vs the latent factors
choiceofK(AIC.null, BIC.null, RSS.null, K = 1:5,
          filename = file.path(p$project_folder, 'codex_results', paste('codex_cls_', p$cluster, '_chr', p$chromosome, '_choiceOfK.pdf', sep='')))

# CBS (recommended) segmentation and CNV calling
finalcall.CBS <- segmentCBS(Y_qc[chr.index,],
                            Yhat.null, optK = which.max(BIC.null),
                            K = 1:5,
                            sampname_qc = colnames(Y_qc),
                            ref_qc = ranges(ref_qc)[chr.index],
                            chr = p$chr, lmax = 400, mode = "integer")

if(p$write_sample_output == 'y') {
  exon_ranges <- import(con = p$gtf,
                        genome = "hs37d5",
                        feature.type = "exon")
  # Annotate (by Michael Schuster) and save CNV results
  for(s in sample_names) {
    codex_table <- subset(as.data.frame(finalcall.CBS), sample_name==s)
    if(length(codex_table$sample_name) > 0) {
      codex_gr <- with(codex_table, GRanges(codex_table$chr, ranges = IRanges(st_bp,ed_bp), copy_no=copy_no, cnv=cnv))
      overlap_frame <- mergeByOverlaps(query = codex_gr, subject = exon_ranges)
      overlap_diagnose_ranges <- overlap_frame$codex_gr
      overlap_diagnose_frame <- mcols(x = overlap_diagnose_ranges)
      overlap_diagnose_frame$gene_id <- mcols(x = overlap_frame$exon_ranges)$gene_id
      overlap_diagnose_frame$gene_name <- mcols(x = overlap_frame$exon_ranges)$gene_name
      mcols(x = overlap_diagnose_ranges) <- overlap_diagnose_frame
      
      # This returns a Grouping object (CompressedManyToOneGrouping) of the IRanges package,
      # specifying which groups contain which indices to the original object.
      message("Group annotated diagnose ranges by region.")
      overlap_diagnose_grouping <-
        as(object = overlap_diagnose_ranges, "Grouping")
      
      grouped_ranges <- unlist(x = GRangesList(lapply(
        X = overlap_diagnose_grouping,
        FUN = function(x) {
          sub_ranges <- overlap_diagnose_ranges[x]
          sub_mcols <- mcols(x = sub_ranges)
          selected_range <- sub_ranges[1L]
          
          mcols(x = selected_range) <- DataFrame(
            "copy_type" = paste(unique(x = sort(x = sub_mcols$cnv)), collapse = ","),
            "copy_number" = paste(unique(x = sort(x = sub_mcols$copy_no)), collapse = ","),
            "gene_names" = paste(unique(x = sort(x = sub_mcols$gene_name)), collapse = ","),
            "gene_ids" = paste(unique(x = sort(x = sub_mcols$gene_id)), collapse = ",")
          )
          return(selected_range)
        }
      )))
      out_name <- paste(s, '_CODEX2_chr', p$chromosome, '.tsv', sep='')
      out_path <- file.path(p$project_folder, 'codex_results', out_name)
      if(!file.exists(out_path) || p$force != 'n' ){
        message(paste('Saving results to', out_path))
        write.table(grouped_ranges, file = out_path, sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)
      }
    }
    else {
      message(paste('There were no CNV calls for sample:', s))
    }
  }
}

out_name <- paste('cls', p$cluster, '_CODEX2_chr', p$chromosome, '.tsv', sep='')
out_path <- file.path(p$project_folder, 'codex_results', out_name)
message(paste('Saving results to', out_path))
write.table(finalcall.CBS, file = out_path, sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)

