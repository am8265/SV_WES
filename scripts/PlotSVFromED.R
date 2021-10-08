#!/home/am5153/miniconda3/envs/gatk/bin/Rscript

#Load Libraries
library(optparse)
library(ExomeDepth)


option_list <- list(
  make_option(c("-i","--variantFile"), type="character", help="tsv file of SV and sample name, e.g. 7_102412845_108524599_DEL 
              /nfs/seqscratch12/am5153/SV_WES/ExomeDepth/IDT_CMA/female/Diagseq2450f892.181998.rds"),
  make_option(c("-f","--from"), type="integer", default="1", help="Index of first test sample"),
  make_option(c("-t","--to"), type="integer", default="1", help="Index of last test sample")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

variantFile = opt$variantFile
#l/outputFile = opt$outputFile
from=opt$f
to=opt$t

variantsBatch = read.table(variantFile, header = F, as.is = T)[from:to, ]

for (row in 1:nrow(variantsBatch)){
  variant = variantsBatch[row, 1]
  print(variant)
  casePath = variantsBatch[row, 2]
  sample = strsplit(basename(casePath), split= ".rds")[[1]][1]
  variants = strsplit(variant, split = "_")[[1]]
  chr = variants[1]
  start = as.numeric(variants[2])
  end = as.numeric(variants[3])
  case = readRDS(casePath)
  outFile = paste(sample, variant, "png", sep = ".")
  png(outFile)
  plot(case, sequence = chr, xlim = c(start - 100000,  end + 100000), count.threshold = 10, cex.lab = 0.8, with.gene =T)
  dev.off()
}
  
  
  
