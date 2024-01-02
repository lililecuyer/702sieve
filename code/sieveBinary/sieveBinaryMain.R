# Purpose: Estimation of hazard ratio-based VE by each Tier 1 and Tier 2 binary features (in SAP Table 5); 
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Lunn and McNeil (1995)
# Input:   702 primary endpoints and binary features
# Output:  csv and PDF files
# Author:  Li Li
# Date:    Dec 20, 2023



rm(list=ls(all=TRUE))
here::i_am("702sieve/rootDir.R")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "702sieve/code/sieveBinary")
outputDir <- file.path(repoDir, "702sieve/code/sieveBinary/output")
figureDir <- file.path(repoDir, "702sieve/figures/sieveBinary")
tableDir <- file.path(repoDir, "702sieve/tables/sieveBinary")


source(file.path(codeDir,"lunnMcneil.R"))
source(file.path(codeDir,"sieveBinaryUtils.R"))
source(file.path(codeDir,"forest.R"))
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(gtable)

loadFile = TRUE

source(file.path(repoDir, "702sieve/code/common.R"))

dat <- dat %>%
  mutate(tx=as.numeric(TRT01P=="T1"), 
         lineageLabel = as.numeric(transmitted.founder.status.tier2=="Multi"),
         eventTime = HIV24fu, 
         eventInd = HIV24)



alvacMatchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier1",".csv")))
residueScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier1",".csv")))
sequonScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varsequonPresentAbsent_tier1",".csv")))
c1086matchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varc1086RmatchvsMismatch_tier1",".csv")))
cTV1matchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varcTV1RmatchvsMismatch_tier1",".csv")))

alvacMatchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier2",".csv")))
residueScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier2",".csv")))
c1086matchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varc1086RmatchvsMismatch_tier2",".csv")))
cTV1matchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varcTV1RmatchvsMismatch_tier2",".csv")))

quantitativeMark_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier1",".csv"))) 
quantitativeMark_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier2",".csv"))) 



adjustList <- list(tier1Type1to3 = c(alvacMatchScreenedIn_tier1$x, residueScreenedIn_tier1$x, sequonScreenedIn_tier1$x),
                   tier1Type4to7 = quantitativeMark_tier1$x,
                   tier2Type1 = residueScreenedIn_tier2$x,
                   tier2Type3 = alvacMatchScreenedIn_tier2$x,
                   tier2Type2and4and5 = quantitativeMark_tier2$x
)
####################################################################################




binaryMarks <- c(alvacMatchScreenedIn_tier1$x, residueScreenedIn_tier1$x, sequonScreenedIn_tier1$x,  c1086matchScreenedIn_tier1$x,
                 cTV1matchScreenedIn_tier1$x, alvacMatchScreenedIn_tier2$x, residueScreenedIn_tier2$x, c1086matchScreenedIn_tier2$x,
                 cTV1matchScreenedIn_tier2$x, "lineageLabel")

continousMarks<- c(quantitativeMark_tier1$x, quantitativeMark_tier2$x)


#use competing risks survival models to obtain VE for all binary marks
if(loadFile == TRUE){
  load(file.path(outputDir, paste0("binaryMarkVE.RData")))
}else{
  binaryResultList <- list()
  for(mark in binaryMarks){
    fitData <- dplyr::select(dat, all_of(c("eventTime","eventInd",mark,"tx")))
    colnames(fitData) <- c("eventTime","eventInd","mark","tx")
    #convert mark values for non-primary endpoints through tau to NA
    fitData$mark[fitData$eventInd==0] <- NA
    
    # complete-case analysis, i.e., discard cases with a missing mark
    fitData <- filter(fitData, !(eventInd==1 & is.na(mark)))
    counts <- data.summary(fitData, "mark")
    #mark-specific cox regression
    fitTable <- tibble(mark = character(), inc = character(), VE = numeric(), LB = numeric(), UB = numeric(), p = numeric())
    
    for(i in c(1,0)){
      
      fitDataMS <- mutate(fitData, newStatus = as.numeric((!is.na(mark)) & mark==i))
      cox<- coxph(Surv(eventTime, newStatus) ~ tx, data=fitDataMS)
      scox <- summary(cox)$coef
      VE <- 1 - scox[2]
      LB <- 1 - exp(scox[1] + qnorm(0.975)*scox[3])
      UB <- 1 - exp(scox[1] - qnorm(0.975)*scox[3])
      p <- scox[5]
      vaccineCases <- counts$vaccineCases[paste(i)]
      placeboCases <- counts$placeboCases[paste(i)]
      vaccineIncRate <- counts$vaccineIncRate[paste(i)]
      placeboIncRate <- counts$placeboIncRate[paste(i)]
      
      fitTable <- add_row(. = fitTable, 
                          mark = ifelse(i==0, "0", "1"),
                          inc = paste0(vaccineCases," (",vaccineIncRate,")", " vs. ",
                                       placeboCases," (",placeboIncRate,")"),
                          VE = VE, LB = LB, UB = UB, p = p
      )
      
    }
    
    lunnMcneilP <- with(fitData, lunnMcneilTest(eventTime,eventInd,mark+1,tx))$coef[3,5]
    binaryResultList[[mark]] <- list(VEtable = fitTable, diffP = lunnMcneilP)
  }
save(binaryResultList, file=file.path(outputDir, paste0("binaryMarkVE.RData")))
}




#Table and graphic summary for tier 2 type 1 and 3
WestfallYoungAdjPvalues_tier2Type3 <- read.csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", paste0("WestfallYoungAdjPvalues_tier2Type3.csv")))
colnames(WestfallYoungAdjPvalues_tier2Type3) <- c("X", "p.unadj","p.FWER","p.FDR")
WestfallYoungAdjPvalues_tier2Type1 <- read.csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", paste0("WestfallYoungAdjPvalues_tier2Type1.csv")))
colnames(WestfallYoungAdjPvalues_tier2Type1 ) <- c("X", "p.unadj","p.FWER","p.FDR")


primaryVE_tier2 = table.VE.match.w.adjP(binaryResultList,marks = alvacMatchScreenedIn_tier2$x,WestfallYoungAdjPvalues_tier2Type3, forestplot = FALSE, Insert = "C.96ZM651")
c1086matchVE_tier2 = table.VE.match.wo.adjP (binaryResultList,marks = c1086matchScreenedIn_tier2$x, forestplot = FALSE, Insert = "C.1086")
cTV1matchVE_tier2 = table.VE.match.wo.adjP (binaryResultList,marks = cTV1matchScreenedIn_tier2$x, forestplot = FALSE, Insert = "C.TV1")
residue_scanning_tier2 <- table.VE.residue(binaryResultList, residueScreenedIn_tier2$x, WestfallYoungAdjPvalues_tier2Type1, forestplot = FALSE)

aggpvaluesPlot (WestfallYoungAdjPvalues_tier2Type3, alvacMatchScreenedIn_tier2$x, adj = TRUE, 
                topText = "Tier 2: Sieve Effect Test P-values (Hazard-based VE)\nMatch/Mismatch C.96ZM651", 
                bottomText = "Env HXB2 Position", sizeText = 16, panelHeight = 2.2,
                yLab = "Differential VE 2-Sided P-value", 
                figureDir = figureDir, figWidth = 8.5, figHeight = 8, 
                fileName = "aggPvaluePlots_matchVsMismatchPrimaryInsert_tier2.pdf")

aggpvaluesPlot (WestfallYoungAdjPvalues_tier2Type1, residueScreenedIn_tier2$x, adj = TRUE, 
                topText = "Tier 2: Sieve Effect Test P-values (Hazard-based VE)\nResidue Present/Absent", 
                bottomText = "Env HXB2 Position", sizeText = 16, panelHeight = 2.2,
                yLab = "Differential VE 2-Sided P-value", 
                figureDir = figureDir, figWidth = 8.5, figHeight = 8, 
                fileName = "aggPvaluePlots_residueScanning_tier2.pdf")


write.csv(primaryVE_tier2, file.path(tableDir, paste0("primaryVE_tier2", ".csv")),row.names = FALSE)
write.csv(c1086matchVE_tier2, file.path(tableDir, paste0("c1086matchVE_tier2", ".csv")),row.names = FALSE)
write.csv(cTV1matchVE_tier2, file.path(tableDir, paste0("cTV1matchVE_tier2", ".csv")),row.names = FALSE)
write.csv(residue_scanning_tier2, file.path(tableDir, paste0("residue_scanning_tier2", ".csv")),row.names = FALSE)

#Forest plot
residue_scanning_tier2_fp <- table.VE.residue(binaryResultList, residueScreenedIn_tier2$x, WestfallYoungAdjPvalues_tier2Type1, forestplot = TRUE)
residue_scanning_tier2_fp$mark <- laply(residue_scanning_tier2_fp$mark, function (x)strsplit(x, split = "\\.")[[1]][2])
#Item-specific summary
aaPos <- unique(laply(colnames(dat)[(grepl("\\.is.",colnames(dat))) &(! (grepl("sequon",colnames(dat))))],function(list)strsplit(list,split="\\.")[[1]][[2]]))
loc <- function(string,pattern){
  seq = seq(1:length(string))
  output = NULL
  for(i in 1:length(pattern)){
    output[i] = seq[string == pattern[i]]
  }
  return(output)
}
tier2TypeA <- aaPos[loc(aaPos,"160"):loc(aaPos,"180")] 
tier2TypeB <- sort(c(aaPos[loc(aaPos,"364"):loc(aaPos,"374")], aaPos[loc(aaPos,"421"):loc(aaPos,"435")], 
                     aaPos[loc(aaPos,"455"):loc(aaPos,"459")], aaPos[loc(aaPos,"466"):loc(aaPos,"477")],
                     aaPos[loc(aaPos,"197"):loc(aaPos,"207")], aaPos[loc(aaPos,"425"):loc(aaPos,"437")],
                     aaPos[loc(aaPos,"274"):loc(aaPos,"283")]))
tier2TypeC <- aaPos[loc(aaPos,"295"):loc(aaPos,"322")] 
tier2TypeD <- c("459", aaPos[loc(aaPos,"466"):loc(aaPos,"470")]) 

forestplotBinary (filter(residue_scanning_tier2_fp, mark %in% tier2TypeA)[,-1], "VE_residue_scanning_tier2_typeA.pdf", figureDir)
forestplotBinary (filter(residue_scanning_tier2_fp, mark %in% tier2TypeB)[,-1], "VE_residue_scanning_tier2_typeB.pdf", figureDir)
forestplotBinary (filter(residue_scanning_tier2_fp, mark %in% tier2TypeC)[,-1], "VE_residue_scanning_tier2_typeC.pdf", figureDir)
forestplotBinary (filter(residue_scanning_tier2_fp, mark %in% tier2TypeD)[,-1], "VE_residue_scanning_tier2_typeD.pdf", figureDir)

# Tier 1 
WestfallYoungAdjPvalues_tier1Type1to3 <- read.csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", paste0("WestfallYoungAdjPvalues_tier1Type1to3.csv")))
colnames(WestfallYoungAdjPvalues_tier1Type1to3) <- c("X", "p.unadj","p.FWER","p.FDR")
primaryVE_tier1 = table.VE.match.w.adjP(binaryResultList,marks = alvacMatchScreenedIn_tier1$x,WestfallYoungAdjPvalues_tier1Type1to3, forestplot = TRUE, Insert = "C.96ZM651")
c1086matchVE_tier1 = table.VE.match.wo.adjP (binaryResultList,marks = c1086matchScreenedIn_tier1$x, forestplot = TRUE, Insert = "C.1086")
cTV1matchVE_tier1 = table.VE.match.wo.adjP (binaryResultList,marks = cTV1matchScreenedIn_tier1$x, forestplot = TRUE, Insert = "C.TV1")
residue_scanning_tier1 <- table.VE.residue(binaryResultList, residueScreenedIn_tier1$x, WestfallYoungAdjPvalues_tier1Type1to3, forestplot = TRUE)

sequon_tier1 <- table.VE.sequon(binaryResultList, sequonScreenedIn_tier1$x, WestfallYoungAdjPvalues_tier1Type1to3, forestplot = TRUE)
aggpvaluesPlot (WestfallYoungAdjPvalues_tier1Type1to3, alvacMatchScreenedIn_tier1$x, adj = TRUE, 
                topText = "Tier 1: Sieve Effect Test P-values (Hazard-based VE)\nMatch/Mismatch C.96ZM651", 
                bottomText = "Env HXB2 Position", sizeText = 16, panelHeight = 2.2,
                yLab = "Differential VE 2-Sided P-value", 
                figureDir = figureDir, figWidth = 8.5, figHeight = 8, 
                fileName = "aggPvaluePlots_tier1.pdf")

forestplotBinary (primaryVE_tier1[,-1], "VE_match_mismatch_tier1.pdf", figureDir, width = 11)
forestplotBinary (residue_scanning_tier1[,-1], "VE_residue_scanning_tier1.pdf", figureDir)
forestplotBinary (sequon_tier1[,-1], "VE_sequon_tier1.pdf", figureDir,width = 10, height1 = 0.23, height2 = 0)
write.csv(c1086matchVE_tier1, file.path(tableDir, paste0("c1086matchVE_tier1", ".csv")),row.names = FALSE)
write.csv(cTV1matchVE_tier1, file.path(tableDir, paste0("cTV1matchVE_tier1", ".csv")),row.names = FALSE)

#single vs multiple lineage
fit <- binaryResultList$lineageLabel
VEtable <- fit$VEtable
VE <- tibble(feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
             P = character(), diffP = character())
VE <- add_row(VE, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
              diffP = format.p(fit$diffP))

VE <- add_row(VE, feature  = "Multiple Founders", cases = VEtable$inc[1],
              VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
              
              mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
              P = format.p(VEtable$p[1]), diffP = "")
VE <- add_row(VE, feature = "Single Founder", cases = VEtable$inc[2],
              VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
              mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
              P = format.p(VEtable$p[2]), diffP = "")
header.plot <- c("Seq\nFeature", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE(%)(95% CI)", "Two-sided\nP-value", "Differential VE\nUnadjusted P-value")
p <- forestWrap ( VE, c(4,5,6), header.plot, xlim = c(-30, 100), 
                  x_ticks_at = c(-20, 0, 20, 40, 60, 80), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black","royalblue", "darkred"),dim(VE)[1]/3),
                  nrows_underline = NULL,
                  insertHeaderText = NULL)

ggplot2::ggsave(filename = "VE_singleVsMultipleFounder.pdf", 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 9, height = 2,units = "in")  

