#!/home/am5153/miniconda3/envs/gatk/bin/Rscript

#Load Libraries
library("optparse")
library(panelcn.mops)
library(data.table)
library(tidyverse)


#Opt parse
option_list <- list(
  make_option(c("-k","--kit"), type="character", default="IDTERPv1", help="exomeKit"),
  make_option(c("-g","--gender"), type="character", default="mixed", help="seqGender"),
  make_option(c("-a","--case"), type="character", help="file with case RDA files, no header"),
  make_option(c("-b","--control"), type="character", help="file with control RDA files, no header"),
  make_option(c("-f","--from"), type="integer", default="1", help="Index of first test sample"),
  make_option(c("-t","--to"), type="integer", default="1", help="Index of last test sample")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


path <- "/nfs/seqscratch12/am5153/SV_WES/scripts/"
path2 <- "/nfs/seqscratch12/am5153/SV_WES/panelcnMOPS/female/"

source(paste0(path, "panelcn.mops.R"))
source(paste0(path, "runPanelcnMops.R"))

kit <- opt$kit
gender <- opt$gender

# kit <- "IDTERPv1"
# gender <- "F"
# gender <- "M"
bed <- "/nfs/seqscratch12/am5153/SV_WES/data/IDT_xGen_Exome_panelcn.mops_new_10bp.bed"
countWindows <- getWindows(bed)

############ load RCS ###############

cr <- fread(opt$case, header = FALSE)$V1

for (f in cr) {
  print(f)
  load(f)
  if (!exists("caseRangesAll")) {
    caseRangesAll <- caseRanges
  } else {
    caseRangesAll@elementMetadata <- cbind(caseRangesAll@elementMetadata, caseRanges@elementMetadata)
  }
}

caseRangesAll@elementMetadata <- caseRangesAll@elementMetadata[,unique(colnames(caseRangesAll@elementMetadata))] 

cr <- fread(opt$control, header = FALSE)$V1

for (f in cr) {
  print(f)
  load(f)
  if (!exists("controlRangesAll")) {
    controlRangesAll <- caseRanges
  } else {
    controlRangesAll@elementMetadata <- cbind(controlRangesAll@elementMetadata, caseRanges@elementMetadata)
  }
}

controlRangesAll@elementMetadata <- controlRangesAll@elementMetadata[,unique(colnames(controlRangesAll@elementMetadata))] 

##### set parameters ########
del <- 0.58 # new (was 0.58)
dup <- 1.5 # new (was 1.46)
I = c(0.025, del, 1, dup, 2)
normType = "quant" # was "quant"
sizeFactor = "quant" # was "quant"
maxControls <- 50
minMedianRC = 30 # was 10
corrThresh = 0.985


selectedGenes <- NULL

if (gender == "F") {
  sex <- "female"
} else if (gender == "M") {
  sex <- "male"
} else {
  sex <- "mixed" 
}

###### run analysis ###########

testiv <- opt$from:opt$to 

logfile <- paste0(path2, "messages_", kit, "_", gender, "_samples_", opt$from, "_" , opt$to, "_", 
                  corrThresh, "_", maxControls, "_", normType, "_", sizeFactor, "_", Sys.Date(), "_del058.log")

con <- file(logfile, "w")
sink(con, append = TRUE, type = "message")
message(timestamp())


XandCB <- caseRangesAll[,testiv]
XandCB@elementMetadata <- cbind(XandCB@elementMetadata, controlRangesAll@elementMetadata)


resultlist <- runPanelcnMops(XandCB, testiv = testiv, countWindows = countWindows, selectedGenes = selectedGenes, I = I, 
                             normType = normType, sizeFactor = sizeFactor, maxControls = maxControls, sex = sex, corrThresh = corrThresh)

message(timestamp())

# sessionInfo()

saveRDS(resultlist, file = paste0(path2, "resultlist_", kit, "_", gender, "_samples_", opt$from, "_" , opt$to, "_",
                                  corrThresh, "_", maxControls, "_", normType, "_", sizeFactor, "_del058_rw150.Rds"))

sink(type="message")
close(con)

sampleNames <- colnames(caseRangesAll@elementMetadata)[testiv]
resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, countWindows = countWindows,
                                 selectedGenes = selectedGenes, sampleNames = sampleNames)


saveRDS(resulttable, file = paste0(path2, "resulttable_", kit, "_", gender, "_samples_", opt$from, "_" , opt$to, "_", 
                                   corrThresh, "_", maxControls, "_", normType, "_", sizeFactor, "_del058_rw150.Rds"))


