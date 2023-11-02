# Setting directory paths -------------------------------------------------
rm(list=ls(all=TRUE))
here::i_am("702sieve/rootDir.R")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "702sieve/code/preScreen")
outputDir <- file.path(repoDir, "702sieve/code/preScreen/output")
figureDir <- file.path(repoDir, "702sieve/figures/preScreen")
tableDir <- file.path(repoDir, "702sieve/tables/preScreen")

library(tidyverse)
library(plyr)
source(file.path(codeDir, "screenUtils.R"))
source(file.path(repoDir, "702sieve/code/common.R"))


sieveData <- read.csv(file.path(dataDir, datFile))

alvacMatch <-  colnames(sieveData)[grepl("1mer.96ZM651.match",colnames(sieveData))]
c1086match <- colnames(sieveData)[grepl("1mer.1086.match",colnames(sieveData))]
cTV1match <- colnames(sieveData)[grepl("1mer.TV1.match",colnames(sieveData))]
posIsAA <- colnames(sieveData)[grepl(".is.",colnames(sieveData)) & grepl("hxb2",colnames(sieveData)) &
                               (!grepl("is.sequon.tier1",colnames(sieveData)))]
  
sequon <- colnames(sieveData)[grepl("is.sequon.tier1",colnames(sieveData))]




#screening heatmaps for match/mismatch
  mos1matchList <- strsplit(mos1match, split = ".", fixed = TRUE)
  mos1matchNewName <- plyr::laply(mos1matchList , function(list){paste0(list[2])})
  mos2matchList <- strsplit(mos2match, split = ".", fixed = TRUE)
  mos2matchNewName <- plyr::laply(mos2matchList , function(list){paste0(list[2])})
  c97zaMatchList <- strsplit(c97zaMatch, split = ".", fixed = TRUE)
  c97zaMatchNewName <- plyr::laply(c97zaMatchList , function(list){paste0(list[2])})
  v1v21428matchList <- strsplit(v1v21428match, split = ".", fixed = TRUE)
  v1v21428matchNewName <- plyr::laply(v1v21428matchList , function(list){paste0(list[2])})
  posIsAAlist <- strsplit(posIsAA, split = ".", fixed = TRUE)
  #remove residues that are vaccine insert in vaccine insert region
  posIsPosLetter <- plyr::ldply(posIsAAlist , function(list){c(list[2],list[4])})
  subInd <- NULL
  
  for(i in 1:length(posIsAA)){
    pos = posIsPosLetter$V1[i]
    residue = posIsPosLetter$V2[i]
    insertMap <- filter(vaccineInsertMap, HXB2.Position==pos)
    if(residue %in% insertMap[1, 3:5]){
      subInd[i] = FALSE
    }else{
      subInd[i] = TRUE
    }
  }
  subPosIsAA <- posIsAA[subInd]
  subPosIsAAlist <- strsplit(subPosIsAA, split = ".", fixed = TRUE)
  posIsAANewName <- plyr::laply(subPosIsAAlist , function(list){paste0(list[2], ".",list[4])})
  
  
  
  seqonAbsenceList <- strsplit(seqonAbsence, split = ".", fixed = TRUE)
  seqonAbsenceNewName <- plyr::laply(seqonAbsenceList , function(list){paste0(list[2])})
  
  
  featureScreen.MatchvsMismatch(sieveData, mos1match, mos1matchNewName, "Mos1", tableDir, figureDir, fileTag)
  featureScreen.MatchvsMismatch(sieveData, mos2match, mos2matchNewName, "Mos2S", tableDir, figureDir, fileTag)
  featureScreen.MatchvsMismatch(sieveData, c97zaMatch, c97zaMatchNewName, "C97ZA", tableDir, figureDir, fileTag)
  featureScreen.MatchvsMismatch(sieveData, v1v21428match, v1v21428matchNewName, "v1v21428", tableDir, figureDir, fileTag)
  featureScreen.is.aa(sieveData, subPosIsAA, posIsAANewName, "posIsAA",tableDir, figureDir, fileTag)
  featureScreen.presenceVSabsence(sieveData, seqonAbsence, seqonAbsenceNewName, tableDir, figureDir, fileTag)
  
  Mos1posScreenedIn <- read.csv(file.path(tableDir, paste0("Mos1posScreenedIn", fileTag,".csv")))
  Mos2posScreenedIn <- read.csv(file.path(tableDir, paste0("Mos2SposScreenedIn", fileTag,".csv")))
  C97ZAposScreenedIn <- read.csv(file.path(tableDir, paste0("C97ZAposScreenedIn", fileTag,".csv")))
  
  posScreenedIn <- rbind(C97ZAposScreenedIn, Mos1posScreenedIn, Mos2posScreenedIn)
  posScreenedIn$Insert <- c("C97ZA Clade C", "Mos1", "Mos2S" )
  write.csv(posScreenedIn, file.path(tableDir, paste0("posScreenedIn", fileTag,".csv")),row.names=FALSE)
  




