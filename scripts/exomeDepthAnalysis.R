#!/home/am5153/miniconda3/envs/gatk/bin/Rscript

#Load Libraries
library(optparse)
library(ExomeDepth)
library(data.table)
library(tidyverse)

#male.control.parents.bamLike.txt
path2 <- "/nfs/seqscratch12/am5153/SV_WES/ExomeDepth/IDT_CMA/LOGS/"

#Opt parse
option_list <- list(
  make_option(c("-k","--kit"), type="character", default="IDTERPv1", help="exomeKit"),
  make_option(c("-g","--gender"), type="character", default="mixed", help="seqGender"),
  make_option(c("-m","--mapFile"), type="character", help="Proband/Sibling map to parents .tsv file, Proband/Sibling - 1st field, Parents - 2nd field onwards"),
  make_option(c("-a","--case"), type="character", help="file with case sample names, no header"),
  make_option(c("-b","--control"), type="character", help="file with control sample names, no header"),
  make_option(c("-c","--countsAll"), type="character", help="control and case counts RDS file"),
  make_option(c("-f","--from"), type="integer", default="1", help="Index of first test sample"),
  make_option(c("-t","--to"), type="integer", default="1", help="Index of last test sample")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

kit <- opt$kit
gender <- opt$gender
readCountsFile <- opt$countsAll
proMapFile <- opt$mapFile
ctrlParents <- opt$control
caseProSib <- opt$case # specify the M, F or mixed group of cases

#coschiz.count 
readCounts = readRDS(readCountsFile)


ctrl.parents.bamLike = readLines(ctrlParents)

case.proSib.bamLike = readLines(caseProSib) #Case



#proMapFile=readLines("../../mixed.ProbandOrSiblingMapParents.txt") 

proMapFile=readLines(proMapFile) 

# Convert a proMap File to a named list

probandMapToParents <- function(x) {
  x=strsplit(x, split = '\t', perl = T)
  names_of_x = sapply(x,'[[',1)
  names(x) = names_of_x
  return(x)
}

proMapFile = probandMapToParents(proMapFile)

samples = colnames(readCounts)[6:length(colnames(readCounts))] # cases and ctrl

# Determine M, F or mixed and create a case list based on that 

if (gender == "F") {
  sex <- "female"
} else if (gender == "M") {
  sex <- "male"
  
} else {
  sex <- "mixed"
}

cases = case.proSib.bamLike # has ending like .realn.recal.bam.rds supplied by user
#cases = sapply(strsplit(cases, split = "\\.", perl =T), 
#               FUN=function(x){paste0(x[1:2],collapse = ".")}) # remove .realn.recal.bam.rds or similar extensions

cases = cases[opt$from : opt$to] # choose cases to run 

allCases = names(proMapFile) # No .realn.recal.bam.rd ending "samples and allCases are similar"--    not used in the current code

###### run analysis ###########

for(i in cases) {
  case_sample = samples[grepl(pattern = i, samples, perl = T)]
  
  i_trunc = strsplit(i, split ="\\.", perl =T)[[1]][1]
  
  family_members = proMapFile[[i_trunc]] 
  
  # ctrl.parents.bamLike is M , F or mixed as provided by the user. Note for Mixed proMapFile could be used 
  logfile <- paste0(path2, "messages_", kit, "_", gender, "_samples_", opt$from, "_" , opt$to, "_", 
                    Sys.Date(), ".log")
  
  con <- file(logfile, "w")
  sink(con, append = TRUE, type = "message")
  message(timestamp())
  message(c("case_sample:",case_sample))
  
  message(c("family_members:",paste0(family_members),collapse = ','))
  
  
  search_family_members = grepl(pattern = paste(family_members, collapse="|"), ctrl.parents.bamLike, perl = T)
  
  control_agg = ctrl.parents.bamLike[!search_family_members]
  control_agg = samples[grepl(pattern = paste(control_agg, collapse="|"), samples, perl = T)]
  message(c("control_agg:",paste0(control_agg, collapse = ',')))
  
  
  readCounts.mat = as.matrix(readCounts[, c(case_sample,control_agg)])
  
  if (is.null(readCounts.mat)) {
    
    message("error in readCounts.mat")
  }
  
  my.choice <- select.reference.set (test.counts = readCounts.mat[,1], reference.counts = readCounts.mat[,-1],
                                     bin.length = (readCounts$end - readCounts$start)/1000, n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = readCounts.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1,
                                 FUN = sum)
  
  message('Now creating the ExomeDepth object')
  
  all.exons <- new('ExomeDepth',
                   test = readCounts.mat[,1],
                   reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')
  
  
  ##Now call the CNVs
  all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
                        chromosome = readCounts$chromosome, start = readCounts$start,
                        end = readCounts$end,
                        name = readCounts$exon)
  
  saveRDS(all.exons, file=paste0(i,".rds"))
  
  CNV.calls = all.exons@CNV.calls[order(all.exons@CNV.calls$BF, decreasing = TRUE),]
 
#  case_name = strsplit(case_sample, ".realn.recal.bam")[[1]]
  sink(type="message")
  close(con)
  

  output.file <- paste(i, 'csv', sep = '.')
  
  write.csv(file = output.file, x = CNV.calls, row.names = FALSE)
}
 
  
  
