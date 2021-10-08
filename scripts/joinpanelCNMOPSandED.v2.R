#!/home/am5153/miniconda3/envs/gatk/bin/Rscript

#Load Libraries
library(optparse)
library(panelcn.mops)
library(data.table)
library(tidyverse)
library(dplyr)
#coschiz416400.167340
#Opt parse
#/nfs/seqscratch12/am5153/SV_WES/scripts/joinpanelCNMOPSandED.R -g M -m 
#/nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/LOGS_2/resulttable_IDTERPv1_mixed_samples_1_50_0.985_50_quant_quant_del058_rw150.Rds -a 
#/nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/LOGS_2/resulttable_IDTERPv1_M_samples_1_32_0.985_50_quant_quant_del058_rw150.Rds -b ED.mixed.txt -c ED.male.txt 
#-f 1 -t 32

option_list <- list(
  make_option(c("-g","--gender"), type="character", help="seqGender M or F "),
  make_option(c("-m","--cnMOPSmixed"), type="character", help="panel.cnMOPS mixed result Rds file path"),
  make_option(c("-a","--cnMOPSmatch"), type="character", help="panel.cnMOPS gender matched result Rds file path"),
  make_option(c("-b","--EDmixed"), type="character", help="file with ED mixed gender.csv full files path, no header"),
  make_option(c("-c","--EDmatch"), type="character", help="file with ED matched gender .csv full files path, no header"),
  make_option(c("-f","--from"), type="integer", default="1", help="Index of first ED matched gender sample"),
  make_option(c("-t","--to"), type="integer", default="1", help="Index of last ED matched gender sample")
)

#-m /nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/LOGS_2/resulttable_IDTERPv1_mixed_samples_1_50_0.985_50_quant_quant_del058_rw150.Rds 
#-a /nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/LOGS_2/resulttable_IDTERPv1_M_samples_1_32_0.985_50_quant_quant_del058_rw150.Rds

# -b ED.mixed.txt
# -c ED.male.txt

# -f 1 
# -t 32


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

gender <- opt$gender
cnMOPS.mixed.path <- opt$cnMOPSmixed
cnMOPS.match.path <- opt$cnMOPSmatch
ED.mixed.list <- opt$EDmixed
ED.match.list <- opt$EDmatch # use this as the reference file

# read panel.cnMOPS data

cnMOPS.mixed.table = readRDS(cnMOPS.mixed.path)
cnMOPS.mixed.samples = sapply(cnMOPS.mixed.table, FUN = function(x){ unique(as.character(x$Sample)) }) # vector of mixed gender sample names

cnMOPS.match.table=readRDS(cnMOPS.match.path)
cnMOPS.match.samples = sapply(cnMOPS.match.table, FUN = function(x){ unique(as.character(x$Sample)) }) # vector of gender matched samples names


# read ED data

ED.mixed.file  <- readLines(ED.mixed.list)
ED.match.file <- readLines(ED.match.list)

## Mixed data
#cnMOPS.mixed.path='/nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/mixed_2/resulttable_IDTERPv1_mixed_samples_1_50_0.985_50_quant_quant_del058_rw150.Rds'

## female data 

#cnMOPS.female.path='/nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/female_2/resulttable_IDTERPv1_F_samples_1_18_0.985_50_quant_quant_del058_rw150.Rds'

# Read the Exome Seq BED file
bed='/nfs/seqscratch12/am5153/SV_WES/data/IDT_xGen_Exome_panelcn.mops_new_10bp.bed'
countWindows <- panelcn.mops::getWindows(bed)


segmentED <- function(res, countWindows){
  
  res$start = res$start - 1 # adjust the start to 0 based to match with panel.cnMOPS
  resGR <- GRanges(seqnames=res$chromosome, IRanges(start=res$start, end=res$end), names=res$id, 
                 mcols=res[,c("start.p", "end.p", "type", "nexons", "BF", "reads.expected", "reads.observed", "reads.ratio")]) 
  
  countWindowsGR = GRanges(seqnames=countWindows$chromosome, IRanges(start=countWindows$start, end=countWindows$end), 
                           names=countWindows$name)
  
  
  hits=findOverlaps(resGR, countWindowsGR)
  
  countWindowsGRHits=countWindowsGR[subjectHits(hits), ]
  countWindowsGRHits@elementMetadata <- cbind(countWindowsGRHits@elementMetadata, 
                                              resGR[queryHits(hits), ]@elementMetadata)
  
  countWindowsGRHits.df = as.data.frame(countWindowsGRHits)
  
  return(countWindowsGRHits.df)
  
  
}

#ED.mixed.cases = basename(ED.mixed.file)

#ED.match.cases = basename(ED.match.file)


# choose matched cases to run
ED.match.file = ED.match.file[opt$from : opt$to]


# For loop for female cases
for (casePath in ED.match.file){
  
  case = basename(casePath)
  
  case = paste0(strsplit(case, split ="\\.")[[1]][c(1,2)], collapse = '.') # case name without extension
  
  
  indexCaseIncnMOPS.mixed.samples = which(grepl(case, cnMOPS.mixed.samples))
  indexCaseIncnMOPS.match.samples = which(grepl(case, cnMOPS.match.samples))

  # indexFemaleIncnMOPS.female.samples = which(grepl(female.case, cnMOPS.male.samples))
  
  case.cnMOPS.mixed = cnMOPS.mixed.table[[ indexCaseIncnMOPS.mixed.samples ]]
  case.cnMOPS.match = cnMOPS.match.table[[ indexCaseIncnMOPS.match.samples ]]
  
# Read ED data and segment it
  case.ED.match = read.table(casePath, sep=',', header =T)
  
  ED.mixed.path = ED.mixed.file[grepl(pattern = case, x = ED.mixed.file, perl = T)]
  
  case.ED.mixed = read.table(ED.mixed.path, sep=',', header =T)
  
  
  case.ED.match = segmentED(res = case.ED.match, countWindows)
  case.ED.mixed = segmentED(res = case.ED.mixed, countWindows)
  
  
  
#join EDs --outerJoin  #Note we  join on svtype --- 

  case.ED.joined = full_join(x = case.ED.match, y = case.ED.mixed, by = c("seqnames", "start", "end"),
                             suffix = c(".matchG.ED", ".mixed.ED"))


# join panel.cnMOPS --outerJoin---we now use CN to join too 

  case.panelcnMOPS.joined = full_join(x= case.cnMOPS.match, y = case.cnMOPS.mixed, by = c("Chr", "Start", "End"), 
                                      suffix = c(".matchG.cnMOPS", ".mixed.cnMOPS"))


# Remove redundant Columns

  case.panelcnMOPS.joined$Sample.mixed.cnMOPS = NULL

# Remove uninformative Columns and sort case.ED.joined
  
  case.ED.joined[c("mcols.start.p.matchG.ED",  "mcols.start.p.mixed.ED")] = NULL

  case.ED.joined = case.ED.joined[order(case.ED.joined[, "seqnames"], case.ED.joined[, "start"], case.ED.joined[, "end"]), ]

  
# Final Join panel.cnMOPS and ED ---more specifically should be joined on svtype i.e. CN with mcols.type
# Another improvement would be to pick may be the larger of 2 columns for names.1.mixed.ED and names.1.matchG.ED (for chr 1 22), 
# for X,Y it should be names.1.matchG.ED --- this will come up when 

  case.cnMOPS.ED.OutJoined = full_join(x= case.panelcnMOPS.joined, y = case.ED.joined, 
                                         by = c("Chr" = "seqnames", "Start" = "start", "End" = "end"), suffix = c("cnMOPS", "ED"))

  case.cnMOPS.ED.InnerJoined = inner_join(x= case.panelcnMOPS.joined, y = case.ED.joined, 
                                            by = c("Chr" = "seqnames", "Start" = "start", "End" = "end"), suffix = c("cnMOPS", "ED"))
  
  
  
  header = colnames(case.cnMOPS.ED.OutJoined) # "names.consensus.ED")
  
  #print(c("header_OUT",header))
  
  ## Reorder Columns and group them by mixed gender or matched gender
  matchG.cnMOPSinCols = which(grepl("matchG.cnMOPS", header, perl =T))
  
  #print(header[matchG.cnMOPSinCols])
  matchG.EDinCols = which(grepl("matchG.ED", header, perl =T))
  
  #print(header[matchG.EDinCols])
  mixed.cnMOPSinCols = which(grepl("mixed.cnMOPS", header, perl =T))
  #print(header[mixed.cnMOPSinCols])
  mixed.EDSinCols = which(grepl("mixed.ED", header, perl =T))
  #print(header[mixed.EDSinCols])
  
  colsOrder = c(2, 5, 6, mixed.cnMOPSinCols, mixed.EDSinCols, matchG.cnMOPSinCols, matchG.EDinCols)
  
  case.cnMOPS.ED.OutJoined = case.cnMOPS.ED.OutJoined[, colsOrder]
  case.cnMOPS.ED.InnerJoined = case.cnMOPS.ED.InnerJoined[, colsOrder]
  

  
  # order by Chr, Start, End 
  
  case.cnMOPS.ED.OutJoined = case.cnMOPS.ED.OutJoined[order(case.cnMOPS.ED.OutJoined[, "Chr"], case.cnMOPS.ED.OutJoined[, "Start"], 
                                                            case.cnMOPS.ED.OutJoined[, "End"]), ]
  
  case.cnMOPS.ED.InnerJoined = case.cnMOPS.ED.InnerJoined[order(case.cnMOPS.ED.InnerJoined[, "Chr"], case.cnMOPS.ED.InnerJoined[, "Start"], 
                                                            case.cnMOPS.ED.InnerJoined[, "End"]), ]

# Write to disk

  write.csv(case.cnMOPS.ED.OutJoined, file = paste0(case, ".outJoin.csv"), quote = F, row.names = F)
  write.csv(case.cnMOPS.ED.InnerJoined, file = paste0(case, ".innerJoin.csv"), quote = F, row.names = F)
}





