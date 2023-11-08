# Purpose: Westfall and Young permutation-based multiplicity adjustment for sieve test p-values
#          The SAP states that p-values are calculated for the vaccine vs. placebo comparison.
# Method:  Westfall and Young (1993)
#          Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.1 on CRAN
#          Lunn and McNeil (1995)
# Input:   702 primary endpoints,  sequence features specified in the sieve SAP that are to be adjusted in multiple comparison
# Output:  A list with each component being a vector of sieve test p-values for the analyzed marks from a single permutation of the mark variable.
# Author:  Li Li

rm(list=ls(all=TRUE))

#p.adj.perm2 needs update
# Setting directory paths -------------------------------------------------
library(here)

rm(list=ls(all=TRUE))
here::i_am("702sieve/rootDir.R")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "702sieve/code/westfallYoung")
outputDir <- file.path(repoDir, "702sieve/code/westfallYoung/output")
figureDir <- file.path(repoDir, "702sieve/figures/westfallYoung")
tableDir <- file.path(repoDir, "702sieve/tables/westfallYoung")

source(file.path(codeDir,"lunnMcneil.R"))
source(file.path(codeDir,"p.adj.perm2.R"))

library(sievePH)
library(tidyverse)
library(plyr)

####################################################################################
#this section of code is the same with sieveBinaryMain.R
source(file.path(repoDir, "702sieve/code/common.R"))

dat <- read.csv(file.path(dataDir, datFile)) %>%
  mutate(tx=as.numeric(TRT01P=="T1"), 
         lineageLabel = as.numeric(transmitted.founder.status.tier2=="Multi"),
         eventTime = HIV24fu, 
         eventInd = HIV24)


alvacMatch_tier1 <-  colnames(dat)[grepl("1mer.96ZM651.match.tier1",colnames(dat))]
c1086match_tier1 <- colnames(dat)[grepl("1mer.1086.match.tier1",colnames(dat))]
cTV1match_tier1 <- colnames(dat)[grepl("1mer.TV1.match.tier1",colnames(dat))]
posIsAA_tier1 <- colnames(dat)[grepl(".is.",colnames(dat)) & grepl("hxb2",colnames(dat)) & grepl("tier1",colnames(dat)) &(!grepl("is.sequon.tier1",colnames(dat)))]
alvacMatch_tier2 <-  colnames(dat)[grepl("1mer.96ZM651.match.tier2",colnames(dat))]
c1086match_tier2 <- colnames(dat)[grepl("1mer.1086.match.tier2",colnames(dat))]
cTV1match_tier2 <- colnames(dat)[grepl("1mer.TV1.match.tier2",colnames(dat))]
posIsAA_tier2 <- colnames(dat)[grepl(".is.",colnames(dat)) & grepl("hxb2",colnames(dat)) & grepl("tier2",colnames(dat)) &(!grepl("is.sequon.tier2",colnames(dat)))]

sequon_tier1 <- colnames(dat)[grepl("is.sequon.tier1",colnames(dat))]
numSequon_tier1 <- colnames(dat)[grepl("num.sequons",colnames(dat))]
lengthHypervariable_tier1 <- colnames(dat)[grepl("length.hypervariable",colnames(dat))]
electrochemicalCharge_tier1 <- colnames(dat)[grepl("charge",colnames(dat))]
cysteineCount_tier1 <- colnames(dat)[grepl("cysteine.count",colnames(dat))]

hdist.zspace_tier1 <- colnames(dat)[grepl("hdist.zspace",colnames(dat)) & grepl("tier1",colnames(dat))]
hdist.zspace_tier2 <- colnames(dat)[grepl("hdist.zspace",colnames(dat)) & grepl("tier2",colnames(dat))]


alvacMatchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier1",".csv")))
c1086matchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varc1086RmatchvsMismatch_tier1",".csv")))
cTV1matchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varcTV1RmatchvsMismatch_tier1",".csv")))
residueScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier1",".csv")))

alvacMatchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier2",".csv")))
c1086matchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varc1086RmatchvsMismatch_tier2",".csv")))
cTV1matchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varcTV1RmatchvsMismatch_tier2",".csv")))
residueScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier2",".csv")))

sequonScreenedIn <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varsequonPresentAbsent_tier1",".csv")))




tier2type1and2 <- c(mos1match, mos2match, c97zaMatch, posIsAA)
tier2posScreenedInList <- list(mos1 = mos1PosScreenedIn$x, mos2 = mos2PosScreenedIn$x, 
                               c97za = c97zaPosScreenedIn$x, 
                               posIsAA = posIsAAPosScreenedIn$x)

#tier 1 match/mismatch
tier1binary <- colnames(dat)[grepl("tier1",colnames(dat))]
tier1posScreenedInList <-  list(mos1 = mos1PosScreenedIn$x[mos1PosScreenedIn$x%in%tier1binary], 
                                mos2 = mos2PosScreenedIn$x[mos2PosScreenedIn$x%in%tier1binary], 
                                c97za = c97zaPosScreenedIn$x[c97zaPosScreenedIn$x%in%tier1binary],
                                v1v21428 = v1v21428PosScreenedIn$x[v1v21428PosScreenedIn$x %in% tier1binary],
                                posIsAA = posIsAAPosScreenedIn$x[posIsAAPosScreenedIn$x %in% tier1binary])
####################################################################################




binaryMarksList <- list( tier1 = c(tier1posScreenedInList$mos1, tier1posScreenedInList$mos2,
                                   tier1posScreenedInList$c97za, tier1posScreenedInList$v1v21428,
                                   tier1posScreenedInList$posIsAA,
                                   sequonPosScreenedIn2,
                                   "lineageLabel"),
                         tier2 = c(tier2posScreenedInList$mos1,tier2posScreenedInList$mos2,
                                   tier2posScreenedInList$c97za, 
                                   tier2posScreenedInList$posIsAA,
                                   "numSequonV5")
)  
continousMarksList <- list(tier1 = c("hdist.zspace.mos1.v2","hdist.zspace.mos2.v2","hdist.zspace.c97za.v2" ,
                                     "hdist.zspace.mos1.v2_ab","hdist.zspace.mos2.v2_ab","hdist.zspace.c97za.v2_ab",
                                     "hdist.zspace.mos1.adcc","hdist.zspace.mos2.adcc","hdist.zspace.c97za.adcc",
                                     "hdist.zspace.mos1.c_ab","hdist.zspace.mos2.c_ab","hdist.zspace.c97za.c_ab",
                                     "hdist.zspace.c97za.hvtn505.cd4bs.antibody", "hdist.zspace.c97za.hvtn505.cd4bs.kmer" ,
                                     "hdist.zspace.1428v1v2",
                                     "length.v1v2",
                                     "num.sequons.v1v2",
                                     "charge.v2",
                                     "cysteine.count"),
                           tier2 = c("hdist.zspace.mos1.gp120","hdist.zspace.mos2.gp120","hdist.zspace.c97za.gp120",
                                     "hdist.zspace.mos1.gp41","hdist.zspace.mos2.gp41","hdist.zspace.c97za.gp41",
                                     "hdist.zspace.mos1.v5","hdist.zspace.mos2.v5","hdist.zspace.c97za.v5",
                                     "length.v5"                                        )
)


#multiple testing adjust groups
adjustList <- list(tier2Type1and2 = c(c97zaPosScreenedIn$x, posIsAAPosScreenedIn$x),
                   tier2Type3and4 = c("hdist.zspace.c97za.gp120","hdist.zspace.c97za.v5", "hdist.zspace.c97za.gp41",
                                      "numSequonV5","length.v5"),
                   tier1Type1to4 = c(c97zaPosScreenedIn$x[c97zaPosScreenedIn$x%in% tier1binary],
                                     v1v21428PosScreenedIn$x,
                                     posIsAAPosScreenedIn$x[posIsAAPosScreenedIn$x%in% tier1binary],
                                     sequonPosScreenedIn2),
                   tier1Type5to7= c("hdist.zspace.c97za.v2", "hdist.zspace.c97za.v2_ab","hdist.zspace.c97za.adcc","hdist.zspace.c97za.c_ab",
                                    "hdist.zspace.c97za.hvtn505.cd4bs.antibody", "hdist.zspace.c97za.hvtn505.cd4bs.kmer" ,
                                    "hdist.zspace.1428v1v2",
                                    "length.v1v2","num.sequons.v1v2","charge.v2","cysteine.count")
)

nPerm <- 1000


# Compute p-values from data sets with resampled marks --------------------
for(adjustGroup in c("tier1Type1to4","tier1Type5to7","tier2Type1and2","tier2Type3and4")){
  marks <- adjustList[[adjustGroup]]
  datSub <- select(dat, all_of(c("eventTime","eventInd","tx","ind_sa", marks)))
  # get the p-values for individual permutations of the observed marks
  # the first vector in 'pvals' stores unadjusted p-values based on original data
  pvals <- lapply(1:(nPerm + 1), function(seed){
    set.seed(seed)
    
    # permute the marks observed or missing in cases
    idx <- sample(1:sum(datSub$eventInd))
    
    # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
    pvals1 <- sapply(1:length(marks), function(i){
      data1 <- select(datSub,all_of(c("eventTime","eventInd","tx","ind_sa", marks[i])))
      colnames(data1) <- c("eventTime","eventInd","tx","ind_sa","mark")
      
      # convert mark values for non-primary endpoints to NA
      data1$mark <- ifelse(data1$eventInd==0, NA, data1$mark)
      
      # apply the permutation
      if (seed>1){
        data1$mark[data1$eventInd==1] <- data1$mark[data1$eventInd==1][idx]  
      }
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
      
      # if the mark is quantitative
      if (marks[i] %in% unlist(continousMarksList)){
        # fit the mark-specific HR model
        markRng <- range(data1$mark, na.rm=TRUE)
        markGrid <- seq(markRng[1], markRng[2], length.out=200)
        fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=ind_sa))
        if(marks[i] %in% hdist.zspace){
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
          # double 1-sided Wald test of {H0: PE(v) constant for all v}
          return(min(2*sfit$pWald.HRconstant.1sided,1))
        }else{
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
          # 2-sided Wald test of {H0: PE(v) constant for all v}
          return(sfit$pWald.HRconstant.2sided) 
        }
        
        
      } else {
        sfit <- with(data1, lunnMcneilTestS(eventTime, eventInd, mark + 1, tx, stratVar=ind_sa))
        return(sfit$coef[3, 5])
      }
    })
    
    names(pvals1) <- marks
    return(pvals1)
  })
  
  pvals <- do.call(rbind, pvals)
  # Apply Westfall and Young (1993) to obtain adjusted p-values -------------
  pvals.adj <- p.adj.perm2(p.unadj=pvals[1, ], p.perms=pvals[-1, ], alpha=1)
  mark <- as.vector(row.names(pvals.adj))
  npvals.adj <- data.frame(mark,pvals.adj )
  write.csv(npvals.adj, file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_",adjustGroup,fileTag,".csv")), row.names=FALSE)
  save(pvals, file=file.path(tableDir, paste0("WestfallYoungPermPvalues_",adjustGroup, fileTag,".RData")))
  
}}







#QQ plots 
for(alignment in c("original","alternative")){
  if(alignment == "original"){
    source(file.path(repoDir, "finalAnalysis/code/common.R"))
  }else{
    source(file.path(repoDir, "finalAnalysis/code/common_alt.R"))
  }
  
  dat <- read.csv(file.path(dataDir, datFile)) %>%
    filter(cohort=="Per-Protocol") %>%
    mutate(tx=as.numeric(armdesc=="Vaccine"), 
           lineageLabel = as.numeric(transmitted.founder.status=="multiple"),
           eventTime = hiv1fposday, 
           eventInd = hiv1event,
           numSequonV5 = as.numeric(num.sequons.v5 >= 2))
  
  mos1match <- colnames(dat)[grepl("1mer.mos1.match",colnames(dat))]
  mos2match <- colnames(dat)[grepl("1mer.mos2.match",colnames(dat))]
  c97zaMatch <- colnames(dat)[grepl("1mer.c97za.match",colnames(dat))]
  v1v21428match <- colnames(dat)[grepl("1mer.1428v1v2.match",colnames(dat))]
  posIsAA <- colnames(dat)[grepl(".is.",colnames(dat))&grepl("hxb2",colnames(dat)) & (!grepl("is.sequon.tier1",colnames(dat)))]
  sequonPos <- colnames(dat)[grepl("is.sequon.tier1",colnames(dat))]
  hdist.zspace <- colnames(dat)[grepl("hdist.zspace",colnames(dat))]
  hdist_ab_no364 <- grep("no364", hdist.zspace, value = TRUE)
  hdist.zspace <- hdist.zspace[!hdist.zspace%in% hdist_ab_no364 ]
  
  
  paste.p <- function(p){if(p<0.001){return("< 0.001")}else{return(as.character(format(p,digits=2)))}}
  tier1Type1to4 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier1Type1to4", fileTag,".csv")))
  tier1Type5to7 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier1Type5to7", fileTag,".csv")))
  tier2Type1and2 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier2Type1and2", fileTag,".csv")))
  tier2Type3and4 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier2Type3and4", fileTag,".csv")))
  adjustedPvalues <- list("tier1Type1to4" = tier1Type1to4,
                          "tier1Type5to7" = tier1Type5to7,
                          "tier2Type1and2" = tier2Type1and2,
                          "tier2Type3and4" = tier2Type3and4)
  
  tier1Positions <- colnames(dat)[grepl("tier1",colnames(dat))]
  tier1aaPos <- laply(tier1Positions[1:67],function(list)strsplit(list,split="\\.")[[1]][[2]])
  
  loc <- function(string,pattern){
    seq = seq(1:length(string))
    output = NULL
    for(i in 1:length(pattern)){
      output[i] = seq[string == pattern[i]]
    }
    return(output)
  }
  
  subTypeA <- tier1aaPos[loc(tier1aaPos,"157"):loc(tier1aaPos,"184")]
  subTypeB <- c("130", "155","160","161","164","165","166","170","172","173","178","179","181","186","200")
  subTypeC <- tier1aaPos[c(loc(tier1aaPos,"51"):loc(tier1aaPos,"61"), loc(tier1aaPos,"69"):loc(tier1aaPos,"78"))]
  subTypeD <- c("169", "181")
  subTypeE <- c("230", tier1aaPos[c(loc(tier1aaPos,"279"):loc(tier1aaPos,"282"))], "350","364","429","432","456","471")
  
  for(adjustGroup in c("tier1Type1to4","tier1Type5to7","tier2Type3and4")){
    
    plotdata <- adjustedPvalues[[adjustGroup]]
    library(EnvStats)
    
    qq.p.unadj = qqPlot(plotdata$p.unadj, distribution = "unif", param.list = list(min = 0, max = 1))
    names(qq.p.unadj) <- c("Uniform", "p.unadj")
    p.unadj.uniformQuantile <- NULL
    
    for(i in 1:length(plotdata$p.unadj)){
      p.unadj.uniformQuantile[i] <- sample(qq.p.unadj$Uniform[qq.p.unadj$p.unadj == plotdata$p.unadj[i]],1)
    }
    plotdata$p.unadj.uniformQuantile <- p.unadj.uniformQuantile
    plotdata$FWERSig <- 1*(plotdata$p.FWER<=0.05)
    plotdata$FDRSig <- 1*(plotdata$p.FDR <= 0.2 & plotdata$p.unadj <= 0.05)
    plotdata$unadjSig <- 1*(plotdata$p.unadj <= 0.05)
    
    zoomx = 0.15
    for(zoom in c(zoomx, 1)){
      pdf(file=file.path(figureDir,paste0("qqplotUnadjustedPvalues", 
                                          adjustGroup,ifelse(zoom == zoomx, "zoom",""),fileTag, ".pdf")),  width=6, height =5)
      
      p <- ggplot(plotdata)+
        geom_point(aes(p.unadj.uniformQuantile, p.unadj), shape = 19,alpha = 0.7) 
      
      p <- p + geom_abline(slope = 1, intercept = 0, size = 0.3, alpha = 0.3)
      
      
      if(zoom == zoomx){
        p <- p + 
          scale_x_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))+
          scale_y_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))
      }else{
        p <- p +
          scale_x_continuous(lim = c(-0.05,1), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))+
          scale_y_continuous(lim = c(-0.05,1), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))
      }
      #yposadd = ifelse(zoom == zoomx, 0.025, 0.5)   
      yposadd = min(plotdata$p.unadj)+ifelse(zoom == zoomx, 0.02, 0.1) 
      for(i in 1:length(plotdata$mark)){
        mark = plotdata$mark[i]
        if(adjustGroup == "tier2Type1and2"){
          if(mark %in% posIsAA){
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markName = paste0("HXB2 site: " ,markSplit [[2]],"; Residue is ", markSplit [[4]],"; ")
            
          }else{
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markName = paste0("HXB2 site: ", markSplit [[2]],"; Residue Match/Mismatch Clade C C97ZA Insert; ")
          }
          
        }
        if(adjustGroup == "tier1Type1to4"){
          if(mark %in% c97zaMatch){
            markSplit = strsplit(mark, split = "\\.")[[1]]
            subType = case_when(markSplit[[2]] %in% subTypeA ~ "subTypeA",
                                markSplit[[2]] %in% subTypeB ~ "subTypeB",
                                markSplit[[2]] %in% subTypeC ~ "subTypeC",
                                markSplit[[2]] %in% subTypeD ~ "subTypeD",
                                markSplit[[2]] %in% subTypeE ~ "subTypeE")
            addDesc = case_when(subType == "subTypeA" ~ "V1V2-hotspot Correlates in NHP Challange Model\n Residue Match/Mismatch Clade C C97ZA Insert",
                                subType == "subTypeB" ~ "IgG3 V1V2 Correlates in Imbokodo\n Residue Match/Mismatch Clade C C97ZA Insert",
                                subType == "subTypeC" ~ "RV144 V2 C1 Sites of Vaccine Pressure in ADCC Epitope\n Residue Match/Mismatch Clade C C97ZA Insert",
                                subType == "subTypeD" ~ "RV144 Primary Sieve Sites\n Residue Match/Mismatch Clade C C97ZA Insert",
                                subType == "subTypeE" ~ "bNab Resistance Signature Sites\n Residue Match/Mismatch Clade C C97ZA Insert")
            
            markName = paste0(addDesc,"\n","HXB2 site: " ,markSplit [[2]],"; ")
          }else if(mark%in% v1v21428match){
            addDesc = "Env gp120 1428 V1V2 Antigen Sequence Residue Match/Mismatch"
            markName = paste0(addDesc,"\n","HXB2 site: " ,markSplit [[2]],"; ")
          }else if (mark %in% posIsAA){
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markSplitnum = gsub("[^[:digit:]]", "", markSplit[[2]])
            subType = case_when(markSplit[[2]] %in% subTypeA ~ "subTypeA",
                                markSplit[[2]] %in% subTypeB ~ "subTypeB",
                                markSplit[[2]] %in% subTypeC ~ "subTypeC",
                                markSplit[[2]] %in% subTypeD ~ "subTypeD",
                                markSplit[[2]] %in% subTypeE ~ "subTypeE")
            addDesc = case_when(subType == "subTypeA" ~ "V1V2-hotspot Correlates in NHP Challange Model",
                                subType == "subTypeB" ~ "IgG3 V1V2 Correlates in Imbokodo",
                                subType == "subTypeC" ~ "RV144 V2 C1 Sites of Vaccine Pressure in ADCC Epitope",
                                subType == "subTypeD" ~ "RV144 Primary Sieve Sites",
                                subType == "subTypeE" ~ "bNab Resistance Signature Sites")
            
            markName = paste0(addDesc,"\n","HXB2 site: " ,markSplit [[2]],"; Residue is ", markSplit [[4]],"; ")
          }else{
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markName = paste0("Sequon with Starting Position at ", markSplit[[2]],";\n")
          }
          
        }
        if(adjustGroup %in% c("tier2Type3and4", "tier1Type5to7")){
          markName <- case_when(mark=="hdist.zspace.mos1.v2" ~ "V1V2 Hotspot Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.v2" ~ "V1V2 Hotspot Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.v2" ~ "V1V2 Hotspot Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.v2_ab" ~ "IgG3 V1V2 CoR Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.v2_ab" ~ "IgG3 V1V2 CoR Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.v2_ab" ~ "IgG3 V1V2 CoR Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.adcc" ~ "V2 C1 ADCC Epitope Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.adcc" ~ "V2 C1 ADCC Epitope Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.adcc" ~ "V2 C1 ADCC Epitope Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.c_ab" ~ "Clade C bNAb Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.c_ab" ~ "Clade C bNAb Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.c_ab" ~ "Clade C bNAb Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.gp120" ~ "gp120 Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.gp120" ~ "gp120 Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.gp120" ~ "gp120 Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.v5" ~ "V5 Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.v5" ~ "V5 Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.v5" ~ "V5 Hamming distance to C97ZA",
                                mark=="hdist.zspace.mos1.gp41" ~ "gp41 Hamming distance to Mos1",
                                mark=="hdist.zspace.mos2.gp41" ~ "gp41 Hamming distance to Mos2S",
                                mark=="hdist.zspace.c97za.gp41" ~ "gp41 Hamming distance to C97ZA",
                                mark=="hdist.zspace.1428v1v2" ~ "Hamming distance from 1428 V1V2 antigen sequence",
                                mark=="hdist.zspace.c97za.hvtn505.cd4bs.antibody" ~ "CD4bs Antibody Contact Sites Hamming distance to C97ZA",
                                mark=="hdist.zspace.c97za.hvtn505.cd4bs.kmer" ~ "CD4bs overlapping k-mers Hamming distance to C97ZA",
                                mark=="length.v1v2" ~ "Length of V1V2",
                                mark=="length.v5" ~ "Length of V5 Loop",
                                mark=="num.sequons.v1v2" ~ "Number of Sequons in V1V2",
                                mark=="charge.v2" ~ "Electrochemical Charge of V2",
                                mark=="numSequonV5" ~ "Number of Sequons in V5",
                                mark == "cysteine.count" ~ "Cysteine count in gp160")
          markName <- paste0(markName,";\n")
        }
        
        
        #arrowbegin.x = ifelse(zoom == zoomx, 0.05+0.001, 0.05+0.06)
        arrowbegin.x = max(plotdata$p.unadj.uniformQuantile[plotdata$p.unadj < 0.05]) +
          ifelse(zoom == zoomx, 0.02, 0.2)
        arrowend.x = plotdata$p.unadj.uniformQuantile[i]
        arrowbegin.y = yposadd
        arrowend.y = ifelse(zoom == zoomx,plotdata$p.unadj[i]+0.001, plotdata$p.unadj[i]+0.01)
        
        
        if(plotdata$unadjSig[i] == 1){
          p <- p +  annotate(
            geom = "curve", x = arrowbegin.x, y = arrowbegin.y,
            xend = arrowend.x, yend = arrowend.y, 
            curvature = .2, arrow = arrow(length = unit(1, "mm")), size = 0.15)
          
          p <- p + annotate(geom = "text", x = arrowbegin.x +0.001, y = arrowbegin.y, 
                            label = paste0(markName,
                                           "FWER P-value ",ifelse(plotdata$p.FWER[i]<0.001,"","= "), 
                                           paste.p(plotdata$p.FWER[i] ),"; ",
                                           "FDR P-value ",ifelse(plotdata$p.FDR[i]<0.001,"","= "), 
                                           paste.p(plotdata$p.FDR[i] )), 
                            hjust = "left", size = 2.8)
          yposadd <- yposadd + ifelse(zoom == zoomx, 0.02, 0.2)
        }
      } 
      
      adjustType <- case_when(adjustGroup=="tier2Type1and2" ~ "Tier 2: Type 1 and 2",
                              adjustGroup=="tier2Type3and4" ~ "Tier 2: Type 3 and 4",
                              adjustGroup=="tier1Type5to7" ~ "Tier 1: Type 5 to Type 7",
                              adjustGroup=="tier1Type1to4" ~ "Tier 1: Type 1 to Type 4")
      
      p <- p+ xlab ("Uniform (0,1) Quantiles")+
        ylab( "Differential VE Unadjusted P-value")+
        ggtitle(paste0("Q-Q Plot of Differential VE Unadjusted P-values \n", adjustType))+
        theme(
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 13),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11, vjust = 0.2),
          axis.title = element_text(size=13)
        )
      print(p)
      dev.off()
    }}
  
  
  
  
  for(adjustGroup in c("tier2Type1and2")){
    
    plotdata <- adjustedPvalues[[adjustGroup]]
    library(EnvStats)
    
    qq.p.unadj = qqPlot(plotdata$p.unadj, distribution = "unif", param.list = list(min = 0, max = 1))
    names(qq.p.unadj) <- c("Uniform", "p.unadj")
    p.unadj.uniformQuantile <- NULL
    
    for(i in 1:length(plotdata$p.unadj)){
      p.unadj.uniformQuantile[i] <- sample(qq.p.unadj$Uniform[qq.p.unadj$p.unadj == plotdata$p.unadj[i]],1)
    }
    plotdata$p.unadj.uniformQuantile <- p.unadj.uniformQuantile
    plotdata$FWERSig <- 1*(plotdata$p.FWER<=0.05)
    plotdata$FDRSig <- 1*(plotdata$p.FDR <= 0.2 & plotdata$p.unadj <= 0.05)
    plotdata$unadjSig <- 1*(plotdata$p.unadj <= 0.05)
    
    zoomx = 0.15
    for(zoom in c(zoomx, 1)){
      pdf(file=file.path(figureDir,paste0("qqplotUnadjustedPvalues", 
                                          adjustGroup,ifelse(zoom == zoomx, "zoom",""),fileTag, ".pdf")),  width=6, height =5)
      
      p <- ggplot(plotdata)+
        geom_point(aes(p.unadj.uniformQuantile, p.unadj), shape = 19,alpha = ifelse(zoom==zoomx, 0.5, 0.05)) 
      
      p <- p + geom_abline(slope = 1, intercept = 0, size = 0.3, alpha = 0.3)
      
      if(zoom == zoomx){
        p <- p + 
          scale_x_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))+
          scale_y_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))
      }else{
        p <- p +
          scale_x_continuous(lim = c(-0.05,1.4), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))+
          scale_y_continuous(lim = c(-0.05,1.4), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))
      }
      
      #yposadd = ifelse(zoom == zoomx, 0.025, 0.5)   
      yposadd = min(plotdata$p.unadj)+ifelse(zoom == zoomx, 0.015, 0.1) 
      for(i in 1:length(plotdata$mark)){
        mark = plotdata$mark[i]
        if(adjustGroup == "tier2Type1and2"){
          if(mark %in% posIsAA){
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markName = paste0("HXB2 site: " ,markSplit [[2]],"; Residue is ", markSplit [[4]],"; ")
            
          }else{
            markSplit = strsplit(mark, split = "\\.")[[1]]
            markName = paste0("HXB2 site: ", markSplit [[2]],"; Residue Match/Mismatch Clade C C97ZA Insert; ")
          }
          
        }
        
        #arrowbegin.x = ifelse(zoom == zoomx, 0.05+0.001, 0.05+0.06)
        if(adjustGroup == "tier2Type1and2" & !zoom == zoomx){
          arrowbegin.x = 0.1
          arrowend.x = plotdata$p.unadj.uniformQuantile[i]
          
          arrowbegin.y = 0.05 + yposadd/2
          
          
          arrowend.y = ifelse(zoom == zoomx,plotdata$p.unadj[i]+0.001, plotdata$p.unadj[i]+0.01)
          
        }else{
          arrowbegin.x = max(plotdata$p.unadj.uniformQuantile[plotdata$p.unadj < 0.05]) -0.015
          arrowend.x = plotdata$p.unadj.uniformQuantile[i]
          arrowbegin.y = yposadd
          arrowend.y = ifelse(zoom == zoomx,plotdata$p.unadj[i]+0.001, plotdata$p.unadj[i]+0.01)
          
        }
        
        if(plotdata$unadjSig[i] == 1){
          p <- p +  annotate(
            geom = "curve", x = arrowbegin.x, y = arrowbegin.y,
            xend = arrowend.x, yend = arrowend.y, 
            curvature = .2, arrow = arrow(length = unit(1, "mm")), size = 0.15)
          
          p <- p + annotate(geom = "text", x = arrowbegin.x +0.001, y = arrowbegin.y, 
                            label = paste0(markName,
                                           "FWER P-value ",ifelse(plotdata$p.FWER[i]<0.001,"","= "), 
                                           paste.p(plotdata$p.FWER[i] ),"; ",
                                           "FDR P-value ",ifelse(plotdata$p.FDR[i]<0.001,"","= "), 
                                           paste.p(plotdata$p.FDR[i] )), 
                            hjust = "left", size = 2.2)
          yposadd <- yposadd + ifelse(zoom == zoomx, 0.01, 0.1)
        }
      } 
      
      adjustType <- case_when(adjustGroup=="tier2Type1and2" ~ "Tier 2: Type 1 and 2",
                              adjustGroup=="tier2Type3and4" ~ "Tier 2: Type 3 and 4",
                              adjustGroup=="tier1Type5to7" ~ "Tier 1: Type 5 to Type 7",
                              adjustGroup=="tier1Type1to4" ~ "Tier 1: Type 1 to Type 4")
      
      p <- p+ xlab ("Uniform (0,1) Quantiles")+
        ylab( "Differential VE Unadjusted P-value")+
        ggtitle(paste0("Q-Q Plot of Differential VE Unadjusted P-values \n", adjustType))+
        theme(
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 13),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11, vjust = 0.2),
          axis.title = element_text(size=13)
        )
      print(p)
      dev.off()
    }
  }  
