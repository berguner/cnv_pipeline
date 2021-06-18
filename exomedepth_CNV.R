#!/usr/bin/env Rscript

library(ExomeDepth)
library(argparser)
library(CODEX2)
library(rtracklayer)
data(Conrad.hg19)

args <- arg_parser(description='Run ExomeDepth CNV analysis for the given samples')
args <- add_argument(args,arg='--bed', default='Y:/lab_bsf/scratch/exome_panel_comparison/gnomad/exome_cnv_calling_regions.bed',
                     help='The path of the target region bed file')
args <- add_argument(args,arg='--project_folder', short='-p',
                     default='Y:/lab_bsf/scratch/exome_panel_comparison/BSA_0005/',
                     help='Project folder path')
args <- add_argument(args,arg='--sample_names', short='-s',
                     default = 'HM_244_S56388,SCID_177_S53796,HM_091_S56377,PID_3141_S53797,HM_246_S56379,HM_020_S56382,HM_209_S53800,PID_3178_S56386,HM_206_S53799,HM_242_S56383,PID_3064_S56380,HM_245_S56385,PID_3187_S56381,HM_212_S53801,PID_3171_S56384,HM_088_S53798,HM_243_S56387',
                     help='Comma separated sample names')
args <- add_argument(args,arg='--cluster', short='-c', default='0',
                     help='Identifier for the batch of samples')
args <- add_argument(args,arg='--force', short='-f', default='y',
                     help='Set this to n if you do not want to overwrite the existing results')
args <- add_argument(args, arg='--assembly',
                     default = 'b37',
                     help='The reference genome version; b37 or hg38')
p <- parse_args(args)

out_folder <- file.path(p$project_folder, 'exomedepth_results')
if(!dir.exists(out_folder)){
  dir.create(out_folder)
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

# Get the GC content and mapabilitiy values for the exonic ranges
gc <- getgc(ref, genome = genome)
#mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc))

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

# Convert the ref dataframe to ExomeDepth's input format
exome_count <- as(ref, 'data.frame')
colnames(exome_count)[colnames(exome_count)=='gc'] <- 'GC'
exome_count$chromosome <- exome_count$seqnames
exome_count$names <- rownames(Y)
for(sample_name in colnames(Y)){
  exome_count$new <- Y[,sample_name]
  colnames(exome_count)[colnames(exome_count) == 'new'] <- sample_name
}

data(exons.hg19)
exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
                                             IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
                                             names = exons.hg19$name)
#case_sample <- 'SCID_177_S53796'

for(case_sample in colnames(Y)){
  output.file <- file.path(out_folder, paste(case_sample, '_ExomeDepth.tsv', sep = ''))
  if(!file.exists(output.file) || p$force != 'n' ){
  message(paste('Running ExomeDepth for sample:', case_sample))
  my.test <- exome_count[[case_sample]]
  my.ref.samples <- colnames(Y)[,1]
  my.reference.set <- as.matrix(exome_count[, my.ref.samples[my.ref.samples != case_sample ]])

  my.choice <- select.reference.set(test.counts = my.test,
                                  reference.counts = my.reference.set,
                                  bin.length = exome_count$width/1000,
                                  n.bins.reduced = 20000)
  message(paste(length(my.choice[[1]]), 'reference samples were picked for sample:', case_sample))
  my.matrix <- as.matrix( exome_count[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  all.exons <- new('ExomeDepth',
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = exome_count$chromosome,
                        start = exome_count$start,
                        end = exome_count$end,
                        name = exome_count$names)
  if(p$assembly == 'b37'){all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.5,
                             column.name = 'Conrad.hg19')
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.hg19.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'exons.hg19')}
  write.table(file = output.file, x = all.exons@CNV.calls, sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)
  }
}
