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
source(file.path(codeDir, "get_minvar.R"))
source(file.path(repoDir, "702sieve/code/common.R"))

dat <- read_csv(file.path(dataDir, datFile))
dat$subjid <- dat$SUBJID
dat$armdesc <- ifelse(dat$TRT01P == "T1", "Vaccine", "Placebo")
dat$hiv1event <- dat$HIV24

tte <- read_csv(file.path(dataDir, "tte_702.csv"))

dat_c <- filter(dat, HIV24==1&!is.na(hxb2.1.1mer.96ZM651.match.tier2))

alvacMatch_tier1 <-  colnames(dat)[grepl("1mer.96ZM651.match.tier1",colnames(dat))]
alvacMatchNewName_tier1 <- plyr::laply(strsplit(alvacMatch_tier1, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
posIsAA_tier1 <- colnames(dat)[grepl(".is.",colnames(dat)) & grepl("hxb2",colnames(dat)) & grepl("tier1",colnames(dat)) &(!grepl("is.sequon.tier1",colnames(dat)))]
posIsAAnewName_tier1 <- plyr::laply(strsplit(posIsAA_tier1, split = ".", fixed = TRUE) , function(list){paste0(list[2],list[4])})
c1086match_tier1 <- colnames(dat)[grepl("1mer.1086.match.tier1",colnames(dat))]
c1086matchNewName_tier1 <- plyr::laply(strsplit(c1086match_tier1, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
cTV1match_tier1 <- colnames(dat)[grepl("1mer.TV1.match.tier1",colnames(dat))]
cTV1matchNewName_tier1 <- plyr::laply(strsplit(cTV1match_tier1, split = ".", fixed = TRUE) , function(list){paste0(list[2])})

sequon_tier1 <- colnames(dat)[grepl("is.sequon.tier1",colnames(dat))]
sequonNewName_tier1 <- plyr::laply(strsplit(sequon_tier1, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
length_tier1 <- colnames(dat)[grepl("length.hypervariable",colnames(dat))]
numSequon_tier1 <- colnames(dat)[grepl("num.sequons",colnames(dat))]
charge_tier1 <- colnames(dat)[grepl("charge",colnames(dat))]
numCysteines_tier1 <- colnames(dat)[grepl("cysteine",colnames(dat)) & grepl("tier1",colnames(dat))]


alvacMatch_tier2 <-  colnames(dat)[grepl("1mer.96ZM651.match.tier2",colnames(dat))]
alvacMatchNewName_tier2 <- plyr::laply(strsplit(alvacMatch_tier2, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
c1086match_tier2 <- colnames(dat)[grepl("1mer.1086.match.tier2",colnames(dat))]
c1086matchNewName_tier2 <- plyr::laply(strsplit(c1086match_tier2, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
cTV1match_tier2 <- colnames(dat)[grepl("1mer.TV1.match.tier2",colnames(dat))]
cTV1matchNewName_tier2 <- plyr::laply(strsplit(cTV1match_tier2, split = ".", fixed = TRUE) , function(list){paste0(list[2])})
posIsAA_tier2 <- colnames(dat)[grepl(".is.",colnames(dat)) & grepl("hxb2",colnames(dat)) & grepl("tier2",colnames(dat)) &(!grepl("is.sequon.tier2",colnames(dat)))]
posIsAAnewName_tier2 <- plyr::laply(strsplit(posIsAA_tier2, split = ".", fixed = TRUE) , function(list){paste0(list[2],list[4])})
hd_tier2 <- colnames(dat)[grepl("hdist",colnames(dat)) & !colnames(dat) %in% c("hdist.zspace.1086.gp120.tier2", "hdist.zspace.TV1.gp120.tier2")]


tab <- table(dat_c$hiv1event, dat_c$TRT01P)
n = dim(dat_c)[1]
n.vx = tab[1,"T1"]
#tier1 binary
tier1.binary.cutoff <- get.minvar (n, n.vx, p=(length(alvacMatch_tier1)+ length(posIsAA_tier1) + length(sequon_tier1)), fwer=0.05) 

#tier 1 quantiative
tier1.quantitative.cutoff <- get.minvar (n, n.vx, p=(length(length_tier1)+ length(numSequon_tier1) + length(charge_tier1) + length(numCysteines_tier1)), fwer=0.05) 

#tier2 binary
tier2.match.cutoff <- get.minvar (n, n.vx, p=(length(alvacMatch_tier2)), fwer=0.05) 
tier2.residue.cutoff <- get.minvar (n, n.vx, p=length(posIsAA_tier2), fwer=0.05)

#tier 2 quantitative
tier1.quantitative.cutoff <- get.minvar (n, n.vx, p=(length(hd_tier2)+ 1), fwer=0.05) 


#screen binary tier 1 marks 
featureScreen.MatchvsMismatch(dat_c, alvacMatch_tier1, alvacMatchNewName_tier1, cutoff = tier1.binary.cutoff, tableDir, figureDir, "primaryRmatchvsMismatch_tier1", "Primary Reference Match\nvs. Mismatch at position")
featureScreen.MatchvsMismatch(dat_c, c1086match_tier1, c1086matchNewName_tier1, cutoff = tier1.binary.cutoff, tableDir, figureDir, "c1086RmatchvsMismatch_tier1", "1086.C Match vs. Mismatch\nat Position")
featureScreen.MatchvsMismatch(dat_c, cTV1match_tier1, cTV1matchNewName_tier1, cutoff = tier1.binary.cutoff, tableDir, figureDir, "cTV1RmatchvsMismatch_tier1", "TV1.C Match vs. Mismatch\nat Position")
featureScreen.presentVsAbsent(dat_c, posIsAA_tier1, posIsAAnewName_tier1, cutoff = tier1.binary.cutoff, tableDir, figureDir, "residuePresentAbsent_tier1", "Amino Acid Present vs.\nAbsent at Position")
featureScreen.presentVsAbsent(dat_c, sequon_tier1, sequonNewName_tier1, cutoff = tier1.binary.cutoff, tableDir, figureDir, "sequonPresentAbsent_tier1", "Presence vs. Absence of a Sequon\nStarting at Position")

#screen binary tier 2 marks 
featureScreen.MatchvsMismatch(dat_c, alvacMatch_tier2, alvacMatchNewName_tier2, cutoff = tier2.match.cutoff, tableDir, figureDir, "primaryRmatchvsMismatch_tier2", "Primary Reference Match\nvs. Mismatch at position")
featureScreen.MatchvsMismatch(dat_c, c1086match_tier2, c1086matchNewName_tier2, cutoff = tier2.match.cutoff, tableDir, figureDir, "c1086RmatchvsMismatch_tier2", "1086.C Match vs. Mismatch\nat Position")
featureScreen.MatchvsMismatch(dat_c, cTV1match_tier2, cTV1matchNewName_tier2, cutoff = tier2.match.cutoff, tableDir, figureDir, "cTV1RmatchvsMismatch_tier2", "TV1.C Match vs. Mismatch\nat Position")
featureScreen.presentVsAbsent(dat_c, posIsAA_tier2, posIsAAnewName_tier2, cutoff =tier2.residue.cutoff, tableDir, figureDir, "residuePresentAbsent_tier2", "Amino Acid Present vs.\nAbsent at Position")


alvacMatchScreenedIn <- read.csv(file.path(tableDir, paste0("primaryRmatchvsMismatch_tier1",".csv")))
c1086matchScreenedIn <- read.csv(file.path(tableDir, paste0("c1086RmatchvsMismatch_tier1",".csv")))
cTV1matchScreenedIn <- read.csv(file.path(tableDir, paste0("cTV1RmatchvsMismatch_tier1",".csv")))
matchScreenedIn <- rbind(alvacMatchScreenedIn, c1086matchScreenedIn, cTV1matchScreenedIn)
matchScreenedIn <- data.frame(insert = c("96ZM651.C + LAI.B ", "1086.C ", "TV1.C" ), matchScreenedIn)
write.csv(matchScreenedIn, file.path(tableDir, paste0("matchVsMismatchScreenedIn_tier1",".csv")),row.names=FALSE)  

alvacMatchScreenedIn <- read.csv(file.path(tableDir, paste0("primaryRmatchvsMismatch_tier2",".csv")))
c1086matchScreenedIn <- read.csv(file.path(tableDir, paste0("c1086RmatchvsMismatch_tier2",".csv")))
cTV1matchScreenedIn <- read.csv(file.path(tableDir, paste0("cTV1RmatchvsMismatch_tier2",".csv")))
matchScreenedIn <- rbind(alvacMatchScreenedIn, c1086matchScreenedIn, cTV1matchScreenedIn)
matchScreenedIn <- data.frame(insert = c("96ZM651.C + LAI.B ", "1086.C ", "TV1.C" ), matchScreenedIn)
write.csv(matchScreenedIn, file.path(tableDir, paste0("matchVsMismatchScreenedIn_tier2",".csv")),row.names=FALSE)  

#screening tier 1 quantitative variables
quantitativeMarks_tier1 <- c(length_tier1, numSequon_tier1, charge_tier1, numCysteines_tier1)
quantitativeMarks_tier1_screened_in <- lapply(quantitativeMarks_tier1, function(x){
  data_mark <- dat_c[, x]
  tab <- table( data_mark)
  #calculate the total of cases that have the mark value that is different from the most frequent mark value
  y <- sum(tab) - max(tab)
  if(y >= tier1.quantitative.cutoff)
    return(TRUE)
  else return(FALSE)
})
