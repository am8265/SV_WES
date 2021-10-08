#!/home/am5153/miniconda3/envs/gatk/bin/Rscript

#Load Libraries
library(optparse)
library(data.table)
library(tidyverse)


option_list <- list(
  make_option(c("-i","--inputFile"), type="character", help="panel.cnMOPS & Exome Depth OuterJoin/Innerjoin .csv file"),
  make_option(c("-f","--from"), type="integer", default="1", help="Index of first test sample"),
  make_option(c("-t","--to"), type="integer", default="1", help="Index of last test sample")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


inputFile = opt$inputFile
#l/outputFile = opt$outputFile
from=opt$f
to=opt$t


inputFilePath = readLines(inputFile)[from:to] # pass to the main functions # run 1:all

## group samples 
samplesPerBatch = 10 # 10 samples/batch ----for. e.g. in 40 samples list there would be 10 samples/batch 
B=samplesPerBatch

noOfSamples = to - from + 1

batch_no = noOfSamples/B # total batches to run
totalBatch = floor(batch_no) 
frac = batch_no - totalBatch 

##How to run a set of samples

if (frac != 0){# i.e 0.7 or 5.7 kinds 
  
  if (totalBatch != 0){
    end=totalBatch
    run_group = lapply(c(1:end), 
                       FUN=function(i=x){j = 10*i - 9; c(j:(i*10))})
    run_group[[end+1]] = c((end*10+1):(end*10 + frac*10))
  } else{
    end=totalBatch
    run_group = list()
    run_group[[end+1]] = c((end*10+1):(end*10 + frac*10))
  }
} else {
  end = totalBatch
  run_group = lapply(c(1:end), 
                     FUN=function(i=x){j = 10*i - 9; c(j:(i*10))}) # there will be totalBatch/end no. of batches each having 10 samples 
}





translateCN <- function(cn) {
  copyNumber = as.integer(str_extract(cn, '[0-9]+'))
  if ((copyNumber == 0)|(copyNumber == 1)) {
    return("deletion")
  } else if (copyNumber >= 3){
    return("duplication")
  } else if (copyNumber == 2){
    return("ref")
  }
}

nonChrXorYAnalysis <- function(sv_exon){
  
  if (!is.na(sv_exon[,"CN.mixed.cnMOPS"])) {
    svtype = translateCN(sv_exon[,"CN.mixed.cnMOPS"])
    
  } else if (!is.na(sv_exon[,"CN.matchG.cnMOPS"])) {
    svtype = translateCN(sv_exon[,"CN.matchG.cnMOPS"])
    
  } else if (!is.na(sv_exon[,"mcols.type.mixed.ED"])) {
    svtype = sv_exon[,"mcols.type.mixed.ED"]
    
  } else if (!is.na(sv_exon[,"mcols.type.matchG.ED"])) {
    svtype = sv_exon[,"mcols.type.matchG.ED"]
  
  }
  
  return(svtype)
}


ChrXorYAnalysis <- function(sv_exon){
  
  if (!is.na(sv_exon[,"CN.matchG.cnMOPS"])) {
    svtype = translateCN(sv_exon[,"CN.matchG.cnMOPS"])
    
  } else if (!is.na(sv_exon[,"mcols.type.matchG.ED"])) {
    svtype = sv_exon[,"mcols.type.matchG.ED"]
    
  } else {
    svtype = "ref"
  }
  
  return(svtype)
}

consensus <- function(intv_mixed, intv_matchG) {
  if (intv_mixed[[3]] == intv_matchG[[3]]) {
    min_start = min(intv_mixed[[1]], intv_matchG[[1]])
    max_stop = max(intv_mixed[[2]], intv_matchG[[2]])
  }
  return(c(min_start,max_stop))
}
    
    


parallelizeAcrossSamples <- function(inputFile) {

  sv_exon = read.csv(inputFile, header = T)

#  header = paste0(c("Chr", "Start", "End", "CN.mixed.cnMOPS", "mcols.type.mixed.ED", "CN.matchG.cnMOPS", "mcols.type.matchG.ED", "svtype"), collapse = "\t")
 
  header = paste0(c(colnames(sv_exon)[c(1:3,10,11)],"svtype"), collapse = "\t")
  
  outputFile = paste0(strsplit(basename(inputFile), split = '.csv')[[1]],".cnMOPsBased.bed")
  
  writeLines(header, outputFile, sep = '\n')

  for (row in 1:nrow(sv_exon)) {
    
    chr = sv_exon[row,"Chr"]
    if ((chr != 'X') & (chr != 'Y')) {
      svtype = nonChrXorYAnalysis(sv_exon[row,]) # return chr start end other fields, svtype: 
      
    } else {
      svtype = ChrXorYAnalysis(sv_exon[row,])
    }
    
    # write output rows here -----
    
    annotSvInput = sv_exon[row,c(1:3,10,11)]
    annotSvInput = cbind(annotSvInput, svtype = svtype)
    write.table(annotSvInput, file = outputFile, append = T, sep = '\t', row.names = F, quote = F, col.names = F)
    
  }
  
}


## parallelizeAcrossSamples2 for ED intervals

parallelizeAcrossSamples2 <- function(inputFile) {
  
  
  
  sv_exon = read.csv(inputFile, header = T)
  
  #  header = paste0(c("Chr", "Start", "End", "CN.mixed.cnMOPS", "mcols.type.mixed.ED", "CN.matchG.cnMOPS", "mcols.type.matchG.ED", "svtype"), collapse = "\t")
  
  #header = paste0(c("chrom", "chromStart", "chromEnd", "svtype", colnames(sv_exon)), collapse = "\t")
  #header = paste0(c("chrom", "chromStart", "chromEnd", "svtype"), collapse = "\t")
  
  
  
  outdir = dirname(inputFile)
   
  sample = strsplit(basename(inputFile), split = '.csv')[[1]]
  outputFile = c(outdir, "/", sample, ".bed")
  outputFile = paste0(outputFile, collapse = "")
  outputFilecnMOPs = c(outdir, "/", sample, "cnMOPs.bed")
  outputFilecnMOPs = paste0(outputFilecnMOPs, collapse = "")
  
  
  
  #writeLines(header, outputFile, sep = '\n')
  
  for (row in 1:nrow(sv_exon)) {
    
    chr = sv_exon[row,"Chr"]
    start_cnMOPs = sv_exon[row, "Start"]
    end_cnMOPs = sv_exon[row, "End"]
    
    if ((chr != 'X') & (chr != 'Y')) {
      
      if (!is.na(sv_exon[row,"names.1.mixed.ED"])) {
        interval = sv_exon[row, "names.1.mixed.ED"]
        
        
        chromStart = as.integer(strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][1]) - 1 # bed adjusted offset by 1
        chromEnd = strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][2]
        
        svtype = sv_exon[row, "mcols.type.mixed.ED"]
        
#        intv_mixed = list(chromStart, chromEnd, svtype)
        
      
    #   } #else { # we assume !is.na(sv_exon[row,"names.1.matchG.ED"])
    #   interval_matchG = sv_exon[row, "names.1.matchG.ED"]
    #   
    #   
    #   chromStart_matchG = as.integer(strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][1]) - 1 # bed adjusted offset by 1
    #   chromEnd_matchG = strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][2]
    #   
    #   svtype_matchG = sv_exon[row, "mcols.type.matchG.ED"]
    #   intv_matchG = list(chromStart_matchG, chromEnd_matchG, svtype_matchG)
    #   
    #   conconsensus(intv_mixed, intv_matchG)
    # }
      
        lowQual_cnMOPS = sv_exon[row, "lowQual.mixed.cnMOPS"]
        CN_cnMOPS = sv_exon[row, "CN.mixed.cnMOPS"]
        
      } else if (!is.na(sv_exon[row,"names.1.matchG.ED"])) {
          
          interval = sv_exon[row, "names.1.matchG.ED"]
        
        
          chromStart = as.integer(strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][1]) - 1 # bed adjusted offset by 1
          chromEnd = strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][2]
        
          svtype = sv_exon[row, "mcols.type.matchG.ED"]
          
          lowQual_cnMOPS = sv_exon[row, "lowQual.mixed.cnMOPS"]
          CN_cnMOPS = sv_exon[row, "CN.mixed.cnMOPS"]
          
      
      } else {
          warning("NA names in mixed.ED and matched.ED")
      }
        
        
      
        
    # For Chr X or chr Y  
    
    } else {
      
      interval = sv_exon[row, "names.1.matchG.ED"]
      
      chromStart = as.integer(strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][1]) - 1
      chromEnd = strsplit(str_extract(interval, '[0-9]+-[0-9]+'), split = '-')[[1]][2]
      
      svtype = sv_exon[row, "mcols.type.matchG.ED"]
      
      lowQual_cnMOPS = sv_exon[row, "lowQual.matchG.cnMOPS"]
      CN_cnMOPS = sv_exon[row, "CN.matchG.cnMOPS"]
    }
    
    
    if (!is.na(interval)) { # We only write nonNA things...a 2nd check/assertion
    
    # write output rows here -----
    
      #annotSvInput = sv_exon[row, ]
      chr = as.character(chr)
      chromStart = as.character(chromStart)
      chromEnd = as.character(chromEnd)
      svtype = as.character(svtype)
      lowQual_cnMOPS = as.character(lowQual_cnMOPS)
      CN_cnMOPS = as.character(CN_cnMOPS)
      start_cnMOPs = as.character(start_cnMOPs)
      end_cnMOPs = as.character(end_cnMOPs)
      
      
      if (row == 1){
        
        annotSvInput = data.frame(chrom = chr, chromStart = chromStart, chromEnd = chromEnd, svtype = svtype, lowQual_cnMOPS = lowQual_cnMOPS, CN_cnMOPS = CN_cnMOPS, start_cnMOPs = start_cnMOPs, end_cnMOPs = end_cnMOPs, stringsAsFactors = F)
        
      } else {
        
        annotSvInput = rbind(annotSvInput, c(chr,chromStart,chromEnd,svtype,lowQual_cnMOPS, CN_cnMOPS, start_cnMOPs, end_cnMOPs), stringsAsFactors = F)
      }
      
    }
  
  } 
  
  annotSvInput$chrom = factor(annotSvInput$chrom, levels = c(c(1:22), "X", "Y"))
  annotSvInput$chromStart = as.numeric(annotSvInput$chromStart)
  annotSvInput$chromEnd = as.numeric(annotSvInput$chromEnd)
  
  # Creating BED based on cnMOPs intervals
  annotSvInput.cnMOPs = annotSvInput[,c("chrom", "start_cnMOPs", "end_cnMOPs", "svtype", "chromStart","chromEnd","lowQual_cnMOPS", "CN_cnMOPS")]
  annotSvInput.cnMOPs$start_cnMOPs = as.numeric(annotSvInput.cnMOPs$start_cnMOPs)
  annotSvInput.cnMOPs$end_cnMOPs = as.numeric(annotSvInput.cnMOPs$end_cnMOPs)
  
  
  annotSvInput = annotSvInput[order(annotSvInput$chrom, annotSvInput$chromStart, annotSvInput$chromEnd), ]
  annotSvInput.cnMOPs = annotSvInput.cnMOPs[order(annotSvInput.cnMOPs$chrom, annotSvInput.cnMOPs$start_cnMOPs, annotSvInput.cnMOPs$end_cnMOPs), ]
  
  
  annotSvInput = annotSvInput[!duplicated(annotSvInput), ]
  annotSvInput.cnMOPs = annotSvInput.cnMOPs[!duplicated(annotSvInput.cnMOPs), ]
    
  write.table(annotSvInput, file = outputFile, append = F, sep = '\t', row.names = F, quote = F, col.names = T)
  write.table(annotSvInput.cnMOPs, file = outputFilecnMOPs, append = F, sep = '\t', row.names = F, quote = F, col.names = T)
  
}




##
    

for (run in run_group) { # i.e. for loop is iterated for each batch, i.e. if we have 40 samples and thus 4 batches, 4 times the for loop is executed and 10 samples are run in parallel/iteration

  run_sampleIndex = run

  runParallelizeAcrossSamples = parallel::mclapply(inputFilePath[run_sampleIndex],
                                                   FUN = parallelizeAcrossSamples2, mc.cores = 30)

}

# for (run in run_group) {
#   
#   run_sampleIndex = run
#   
#   runParallelizeAcrossSamples = parallel::mclapply(inputFilePath[run_sampleIndex],
#                                                    FUN = parallelizeAcrossSamples2, mc.cores = 30)
#   
# }


# for (run in run_group) {
# 
#   run_sampleIndex = run
# 
#   runParallelizeAcrossSamples = lapply(inputFilePath[run_sampleIndex],
#                                                    FUN = parallelizeAcrossSamples2)
# 
# }

  
 