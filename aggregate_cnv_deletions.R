##################################
## The aim of this script is combine and annotate the CNV calls in a cohort.
## For this purpose, single sample CNV calls are collected, the number of unique samples per GENE counted.
## Finally, the sample counts per GENE are then annotated for each sample CNV call.
## Author: Matthias Haimel
##################################

library(argparser)
options(stringsAsFactors=FALSE)

## get variables from ENV or STOP
get_checked_env <- function(id, dofail=TRUE){
  val = Sys.getenv(id)
  if(val == "" & dofail){
    stop(paste(id,"variable not found"))
  }
  return(val)
}

args <- arg_parser(description='Aggregate CNV results')
args <- add_argument(args, arg='--project_id',
                     default='BSA_0005_KBL_Exome',
                     help='Project name/ID')
args <- add_argument(args, arg='--sample_sheet',
                     default='X:\\home\\mschuster\\src\\bsfconfiguration\\variant_calling\\BSA_0005_KBL_Exome_variant_calling_samples.csv',
                     help='Sample annotation sheet')
args <- add_argument(args, arg='--project_folder',
                     default='Y:\\lab_bsf\\projects\\BSA_0005_KBL_Exome\\b37_cnv_pipeline',
                     help='Path of the folder containing the cnv_pipeline results')
args <- add_argument(args, arg='--gtf',
                     default = 'Y:/lab_bsf/resources/genomes/GRCh37_e87/Homo_sapiens.GRCh37.87.gtf.gz',
                     help='The path of the Ensembl GTF file')
args <- add_argument(args, arg='--all_gene_symbol_file',
                     default='X:\\groups\\lab_bsf\\resources\\cnv_data\\boztug_reference\\2019-11-22_ibd_pid_hm_genes.found.txt',
                     help='genes from gene panel')
args <- add_argument(args, arg='--canonical_file',
                     default='X:\\groups\\lab_bsf\\resources\\cnv_data\\boztug_reference\\canonical_transcripts.csv',
                     help='Canonical transcripts exported through REST api')
args <- add_argument(args, arg='--pli_score_file',
                     default='X:\\groups\\lab_bsf\\resources\\cnv_data\\boztug_reference\\fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                     help='ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt')
args <- add_argument(args, arg='--hi_score_file',
                     default='X:\\groups\\lab_bsf\\resources\\cnv_data\\boztug_reference\\HI_Predictions_Version3.bed',
                     help='https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz')
args <- parse_args(args)


project <- args$project_id #get_checked_env("PROJECT_ID") # e.g. BSA_0005_KBL_Exome
mft.file <- args$sample_sheet #get_checked_env("MANIFEST") # /path/to/BSA_0005_KBL_Exome_variant_calling_samples.csv
## Output directory, where the results will be stored in
output_dir <- file.path(args$project_folder, 'aggregated') #get_checked_env("CNV_OUT_DIR") # /path/to/BSA_.../b37/cnv_analysis/aggregated
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

project_folder <- args$project_folder

## optional parameters for cross cohort annotation
external_annot_file <- get_checked_env("EXT_ANN_FILE", dofail=FALSE)
external_annot_name <- get_checked_env("EXT_ANN_PREFIX", dofail=FALSE)

all_gene_symbol_file <- args$all_gene_symbol_file  ## genes from panel panel
canonical_file <- args$canonical_file ## exported through REST api
pli_score_file <- args$pli_score_file ## ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
hi_score_file <- args$hi_score_file ## https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz

gtf <- args$gtf

### Load libraries
print(paste("Process",project,"..."))
suppressPackageStartupMessages(expr = require(reshape2))
suppressPackageStartupMessages(expr = require(GenomicRanges))
suppressPackageStartupMessages(expr = require(package = "rtracklayer"))

## Read configuration files
external_annot_db <- data.frame(gene_name=c())
external_annot_db_columns <- c()
if (file.exists(external_annot_file)) {
  external_annot_db <- read.delim(external_annot_file)
  external_annot_db_columns <- c(
    paste(external_annot_name,"sample_hom_cnt",sep=""), 
    paste(external_annot_name,"sample_het_cnt",sep="")
  )
  colnames(external_annot_db) <- c("gene_name", external_annot_db_columns)
}

mft <- read.delim(mft.file, sep=",", stringsAsFactors = FALSE)
pli <- read.delim(pli_score_file)
genelist <- read.delim(all_gene_symbol_file, header=F, stringsAsFactors = FALSE)$V1
hi <- as.data.frame(do.call("rbind",strsplit(read.delim(hi_score_file, skip = 1, header=F, stringsAsFactors=FALSE)$V4, "\\|")), stringsAsFactors=FALSE)
colnames(hi) <- c("gene_name", "hi_score", "pc")
hi$hi_score <- as.numeric(hi$hi_score)
hi <- hi[order(hi$hi_score, decreasing=T),]

canonical_transcript <- read.delim(canonical_file, sep=",")
canonical_transcript$is_canonical <- "YES"

#####
## Load resources
## for annotation of gene level
gene_ranges <-
  import(con = gtf,
         genome = "hs37d5",
         feature.type = "gene") 
geneid_to_annot <- mcols(gene_ranges)[c("gene_id", "gene_name", "gene_biotype")]

## Load ensembl gene list to find protein coding genes
exon_ranges <-
  import(con = gtf,
         genome = "hs37d5",
         feature.type = "exon")

# exon_ranges.ensembl <- exon_ranges[mcols(exon_ranges)$source == "ensembl" & mcols(exon_ranges)$gene_biotype == "protein_coding"]
exon_ranges.ensembl <- exon_ranges[mcols(exon_ranges)$source %in% c("ensembl","ensembl_havana", "havana")]
# exon_ranges.ensembl

######
## find CNV files for each pipeline
unique_cols <- c("Sample.Name","Sample.Project.Name")
mft <- mft[!mft$PairedReads.Exclude,]
target <- unique(mft[,unique_cols])
target$codex <- file.path(project_folder, 'codex_results', 'sample_results', paste(target$Sample.Name,"_CODEX2.tsv", sep=""))
target$exome_depth <- file.path(project_folder, 'exomedepth_results', paste(target$Sample.Name,"_ExomeDepth.tsv", sep=""))
# check if they exist
target$codex_exist <- file.exists(target$codex)
target$exome_depth_exist <- file.exists(target$exome_depth)

#######
### READ in data

## helper function to read and annotate data with sample name
combine_data <- function(x, file_param){
  target_file <- x[file_param]
  if(!file.exists(target_file)) {
    return(NA)
  }
  
  d <- read.delim(target_file, quote='"', sep="\t")
  if(nrow(d) < 1){
    return(NA)
  }
  d$name <- x["Sample.Name"]
  return(d)
}

## load & combine data into one dataframe per analysis
data.codex <- do.call("rbind", apply(target, 1, function(x){combine_data(x, "codex")}))
data.exome_depth <- do.call("rbind", apply(target, 1, function(x){combine_data(x, "exome_depth")}))

## only select deletions
data.codex.del <- subset(data.codex, cnv == "del")
data.exome_depth.del <- subset(data.exome_depth, type == "deletion")

######
## Find overlaps

## transfer into GRanges

remove_chr <- function(chr){
  return (sub("chr","",chr))
}

range.exome_depth <- GRanges(
  remove_chr(data.exome_depth.del$chromosome), 
  IRanges(data.exome_depth.del$start, data.exome_depth.del$end),'*',
  sample_name=data.exome_depth.del$name,
  reads.ratio=data.exome_depth.del$reads.ratio,
  BF=data.exome_depth.del$BF
)

range.codex <- GRanges(
  remove_chr(data.codex.del$chr), 
  IRanges(data.codex.del$st_bp, data.codex.del$ed_bp),'*',
  sample_name=data.codex.del$name,
  copy_no=data.codex.del$copy_no,
  lratio=data.codex.del$lratio,
  seqinfo=seqinfo(range.exome_depth)
)

#####################
## find Overlaps between Exons and Deletions

## First CodeX
hits <- findOverlaps(query = exon_ranges.ensembl, subject = range.codex)
hits.ens <- exon_ranges.ensembl[queryHits(hits)]
hits.codex <- range.codex[subjectHits(hits)]
## save as data frame for later merging
codex.df <- data.frame(chr=seqnames(hits.ens), start=start(hits.ens), end=end(hits.ens),
                          codex_del_id=paste(seqnames(hits.codex),":", start(hits.codex),"-", end(hits.codex), sep=""),
                          sample_name = mcols(hits.codex)$sample_name,
                          codex_lratio=mcols(hits.codex)$lratio, codex_copy_number=mcols(hits.codex)$copy_no,
                          mcols(hits.ens)[c("gene_id", "transcript_id", "exon_number", "gene_name", "transcript_biotype", "ccds_id", "gene_biotype")])

codex.gene_calls <- unique(codex.df[, c("sample_name", "gene_id", "codex_del_id", "codex_lratio", "codex_copy_number")])
codex.gene_calls$call_id <- codex.gene_calls$codex_del_id
codex.gene_calls$codex_del_id <- NULL
# codex.gene_calls.del_ids <- aggregate(codex_del_id ~  sample_name + gene_id , codex.gene_calls, paste)
# codex.gene_calls.ratio <- aggregate(codex_lratio ~  sample_name + gene_id , codex.gene_calls, paste)

## Then ExomeDepth
hits <- findOverlaps(query = exon_ranges.ensembl, subject = range.exome_depth)
hits.ens <- exon_ranges.ensembl[queryHits(hits)]
hits.edep <- range.exome_depth[subjectHits(hits)]
edep.df <- data.frame(chr=seqnames(hits.ens), start=start(hits.ens), end=end(hits.ens),
                       exomedepth_del_id=paste(seqnames(hits.edep),":", start(hits.edep),"-", end(hits.edep), sep=""),
                       sample_name = mcols(hits.edep)$sample_name,
                       exomedepth_bf=mcols(hits.edep)$BF, exomedepth_ratio=mcols(hits.edep)$reads.ratio,
                       mcols(hits.ens)[c("gene_id", "transcript_id", "exon_number", "gene_name", "transcript_biotype", "ccds_id", "gene_biotype")])

edep.gene_calls <- unique(edep.df[, c("sample_name", "gene_id", "exomedepth_del_id", "exomedepth_bf", "exomedepth_ratio")])
edep.gene_calls$call_id <- edep.gene_calls$exomedepth_del_id
edep.gene_calls$exomedepth_del_id <- NULL

## get complete list of position sets
gene_call_with_ids <- unique(rbind(
  codex.gene_calls[, c("sample_name", "gene_id", "call_id")],
  edep.gene_calls[, c("sample_name", "gene_id", "call_id")]
))

#######
## Merge gene level
## add corresponding calls
gene_data_annot <- merge(gene_call_with_ids, codex.gene_calls, by=c("sample_name", "gene_id", "call_id"), all=TRUE)
gene_data_annot <- merge(gene_data_annot, edep.gene_calls, by=c("sample_name", "gene_id", "call_id"), all=TRUE)

## add Gene annotation
gene_data_annot <- merge(gene_data_annot, geneid_to_annot, by="gene_id")

## add stats
gene_sample <- as.data.frame(unique(gene_data_annot[,c("gene_id", "gene_name","sample_name")]))
cnt_genes <- dcast(gene_sample, gene_name + gene_id ~ . , length, value.var = "gene_id")
colnames(cnt_genes) <- c("gene_name", "gene_id","sample_cnt")
cnt_genes <- merge(cnt_genes, pli[,c("gene","pLI", "lof_z", "mis_z", "pRec", "pNull", "oe_lof_upper")], by.x="gene_name", by.y="gene", all.x=TRUE)
cnt_genes <- merge(cnt_genes, hi[,c("gene_name","hi_score")], by.x="gene_name", by.y="gene_name", all.x=TRUE)
## add gene panel flag for easier filtering
cnt_genes$in_gene_panel <- cnt_genes$gene_name %in% genelist
cnt_genes$gene_name <- NULL

## add annotation
gene_data_annot <- merge(gene_data_annot, cnt_genes, by.x="gene_id", by.y="gene_id", all.x=TRUE)

cols=c("call_id", "gene_id", "gene_name", "sample_name", "sample_cnt", "in_gene_panel",
       "codex_lratio", "codex_copy_number", "exomedepth_bf", "exomedepth_ratio", 
       "hi_score", "pLI", "pRec", "pNull", "oe_lof_upper",
       "gene_biotype"
)

## comprehensive version - with all information
write.table(
  file=file.path(output_dir, "samples.cnv_gene_all.tsv"),
  sep="\t",
  na = "",
  row.names = FALSE,
  gene_data_annot[,  cols]
)

cols=c("call_id", "gene_id", "gene_name", "sample_name", "sample_cnt", "in_gene_panel",
       "codex_lratio", "codex_copy_number", "exomedepth_bf", "exomedepth_ratio", 
       "hi_score", "pLI", "pRec", "pNull", "oe_lof_upper"
)

## comprehensive version - with all information
write.table(
  file=file.path(output_dir, "samples.cnv_gene_protein_coding.tsv"),
  sep="\t",
  na = "",
  row.names = FALSE,
  subset(gene_data_annot, gene_biotype == "protein_coding",  cols)
)

