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

dat <- dat %>%
  mutate(tx=as.numeric(TRT01P=="T1"), 
         eventTime = HIV24fu, 
         eventInd = HIV24)



alvacMatchScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier1",".csv")))
residueScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier1",".csv")))
sequonScreenedIn_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varsequonPresentAbsent_tier1",".csv")))

alvacMatchScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varprimaryRmatchvsMismatch_tier2",".csv")))
residueScreenedIn_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varresiduePresentAbsent_tier2",".csv")))

quantitativeMark_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier1",".csv"))) 
quantitativeMark_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier2",".csv"))) 





binaryMarks <- c(alvacMatchScreenedIn_tier1$x, residueScreenedIn_tier1$x, sequonScreenedIn_tier1$x, alvacMatchScreenedIn_tier2$x, residueScreenedIn_tier2$x)

continousMarks<- c(quantitativeMark_tier1$x, quantitativeMark_tier2$x)



#multiple testing adjust groups
adjustList <- list(tier1Type1to3 = c(alvacMatchScreenedIn_tier1$x, residueScreenedIn_tier1$x, sequonScreenedIn_tier1$x),
                   tier1Type4to7 = quantitativeMark_tier1$x,
                   tier2Type1 = residueScreenedIn_tier2$x,
                   tier2Type3 = alvacMatchScreenedIn_tier2$x,
                   tier2Type2and4and5 = quantitativeMark_tier2$x
)

nPerm <- 1000


# Compute p-values from data sets with resampled marks --------------------
for(adjustGroup in c("tier1Type1to3","tier1Type4to7","tier2Type1", "tier2Type3","tier2Type2and4and5")){
  marks <- adjustList[[adjustGroup]]
  datSub <- select(dat, all_of(c("eventTime","eventInd","tx", marks)))
  # get the p-values for individual permutations of the observed marks
  # the first vector in 'pvals' stores unadjusted p-values based on original data
  pvals <- lapply(1:(nPerm + 1), function(seed){
    set.seed(seed)
    
    # permute the marks observed or missing in cases
    idx <- sample(1:sum(datSub$eventInd))
    
    # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
    pvals1 <- sapply(1:length(marks), function(i){
      data1 <- select(datSub,all_of(c("eventTime","eventInd","tx", marks[i])))
      colnames(data1) <- c("eventTime","eventInd","tx","mark")
      
      # convert mark values for non-primary endpoints to NA
      data1$mark <- ifelse(data1$eventInd==0, NA, data1$mark)
      
      # apply the permutation
      if (seed>1){
        data1$mark[data1$eventInd==1] <- data1$mark[data1$eventInd==1][idx]  
      }
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
      
      # if the mark is quantitative
      if (marks[i] %in% unlist(continousMarks)){
        # fit the mark-specific HR model
        markRng <- range(data1$mark, na.rm=TRUE)
        markGrid <- seq(markRng[1], markRng[2], length.out=200)
        fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx))
        if(marks[i] %in% hd_tier2){
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
          # double 1-sided Wald test of {H0: PE(v) constant for all v}
          return(min(2*sfit$pWald.HRconstant.1sided,1))
        }else{
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
          # 2-sided Wald test of {H0: PE(v) constant for all v}
          return(sfit$pWald.HRconstant.2sided) 
        }
        
        
      } else {
        sfit <- with(data1, lunnMcneilTest(eventTime, eventInd, mark + 1, tx))
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
  write.csv(npvals.adj, file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_",adjustGroup,".csv")), row.names=FALSE)
  save(pvals, file=file.path(tableDir, paste0("WestfallYoungPermPvalues_",adjustGroup,".RData")))
  
}




#QQ plots 
paste.p <- function(p){if(p<0.001){return("< 0.001")}else{return(as.character(format(p,digits=2)))}}
tier1Type1to3 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier1Type1to3",".csv")))
tier1Type4to7 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier1Type4to7",".csv")))
tier2Type1 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier2Type1",".csv")))
tier2Type3 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier2Type3",".csv")))
tier2Type2and4and5 <- read.csv(file.path(tableDir,paste0("WestfallYoungAdjPvalues_tier2Type2and4and5",".csv")))
adjustedPvalues <- list("tier1Type1to3" = tier1Type1to3,
                        "tier1Type4to7" = tier1Type4to7,
                        "tier2Type1" = tier2Type1,
                        "tier2Type3" = tier2Type3,
                        "tier2Type2and4and5" = tier2Type2and4and5)


aaPos <- unique(laply(colnames(dat)[(grepl("\\.is.",colnames(dat))) &(! (grepl("sequon",colnames(dat))))],function(list)strsplit(list,split="\\.")[[1]][[2]]))

loc <- function(string,pattern){
  seq = seq(1:length(string))
  output = NULL
  for(i in 1:length(pattern)){
    output[i] = seq[string == pattern[i]]
  }
  return(output)
}

# Tier 2 sites that overlap with tier 1 will be excluded in tier 2 analyses
tier1TypeA <- sort(c("6", "19", "169", "181", "307", "268", "307", "317", "343", "353", "369", "379"))
tier1TypeB <- c("97", "122", "364", "369", "424")
tier2TypeA <- aaPos[loc(aaPos,"160"):loc(aaPos,"180")] 
tier2TypeB <- sort(c(aaPos[loc(aaPos,"364"):loc(aaPos,"374")], aaPos[loc(aaPos,"421"):loc(aaPos,"435")], 
                aaPos[loc(aaPos,"455"):loc(aaPos,"459")], aaPos[loc(aaPos,"466"):loc(aaPos,"477")],
                aaPos[loc(aaPos,"197"):loc(aaPos,"207")], aaPos[loc(aaPos,"425"):loc(aaPos,"437")],
                aaPos[loc(aaPos,"274"):loc(aaPos,"283")]))
tier2TypeC <- aaPos[loc(aaPos,"295"):loc(aaPos,"322")] 
tier2TypeD <- c("459", aaPos[loc(aaPos,"466"):loc(aaPos,"470")]) 

for(adjustGroup in c("tier1Type1to3","tier1Type4to7","tier2Type1", "tier2Type3","tier2Type2and4and5")){
  
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
  
  zoomx = 0.22
  for(zoom in c(zoomx, 1)){
    pdf(file=file.path(figureDir,paste0("qqplotUnadjustedPvalues", 
                                        adjustGroup,ifelse(zoom == zoomx, "zoom",""), ".pdf")),  width=6, height =5)
    
    p <- ggplot(plotdata)+
      geom_point(aes(p.unadj.uniformQuantile, p.unadj), shape = 19,alpha = 0.7) 
    
    p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.3, alpha = 0.3)
    
    
    if(zoom == zoomx){
      p <- p + 
        scale_x_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))+
        scale_y_continuous(lim = c(0,zoomx), breaks = c(0,0.05,0.1,0.15), minor_breaks = c(0,0.05,0.1,0.15))
    }else{
      p <- p +
        scale_x_continuous(lim = c(-0.05,1.1), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))+
        scale_y_continuous(lim = c(-0.05,1.2), breaks = c(0,0.25,0.5,0.75,1), minor_breaks = c(0,0.25,0.5,0.75,1))
    }
    #yposadd = ifelse(zoom == zoomx, 0.025, 0.5)   
    yposadd = min(plotdata$p.unadj)+ifelse(zoom == zoomx, 0.02, 0.1) 
    for(i in 1:length(plotdata$mark)){
      mark = plotdata$mark[i]
      if(adjustGroup %in% c("tier1Type1to3", "tier2Type1", "tier2Type3")){
        if(mark %in% c(posIsAA_tier1, posIsAA_tier2)){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("HXB2 site: " ,markSplit [[2]],"; Residue is ", markSplit [[4]],"; ")
          
        }else if (mark %in% c(alvacMatch_tier1, alvacMatch_tier2)){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("HXB2 site: ", markSplit [[2]],"; Residue Match/Mismatch C.96ZM651; ")
        }else if (mark %in% c(sequon_tier1)){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("HXB2 site: " ,markSplit [[2]],"; Sequon Present/Absent; ")
        }
        
      }
      
      if(adjustGroup %in% c("tier1Type4to7")){
        if(mark %in% numSequon_tier1){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("Number of Sequons in " ,markSplit [[3]])
          
        }else if (mark %in% length_tier1){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("The Length of", markSplit [[3]])
        }else if (mark %in% charge_tier1 ){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("The Total Electrochemical Charge of " ,markSplit [[2]])
        }else if (mark %in% numCysteines_tier1){
          markSplit = strsplit(mark, split = "\\.")[[1]]
          markName = paste0("The Number of Cysteine in " ,markSplit [[3]])
        }
        
      }
      if(adjustGroup %in% c("tier2Type2and4and5")){
        markName <- case_when(mark=="cysteine.count.gp160.tier2" ~ "Cysteine count in gp160",
                              mark=="hdist.zspace.96ZM651.V2.linear.epitope.tier2" ~ "V2 Linear Epitope PC-weighted Hamming Distance to C.96ZM651",
                              mark=="hdist.zspace.96ZM651.CD4.binding.site.tier2" ~ "CD4 Contact Epitopes PC-weighted Hamming Distance to C.96ZM651",
                              mark=="hdist.zspace.96ZM651.V5.tier2" ~ "V5 PC-weighted Hamming Distance to C.96ZM651",
                              mark=="hdist.zspace.96ZM651.V3.tier2" ~ "V3 PC-weighted Hamming Distance to C.96ZM651",
                              mark=="hdist.zspace.96ZM651.tier2" ~ "PC-weighted Hamming Distance to C.96ZM651")
        markName <- paste0(markName,";\n")
      }
      
      
      #arrowbegin.x = ifelse(zoom == zoomx, 0.05+0.001, 0.05+0.06)
      arrowbegin.x = max(plotdata$p.unadj.uniformQuantile[plotdata$p.unadj < 0.05]) +
        ifelse(zoom == zoomx, 0.02, 0.15)
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
                                         "\nFWER P-value ",ifelse(plotdata$p.FWER[i]<0.001,"","= "), 
                                         paste.p(plotdata$p.FWER[i] ),"; ",
                                         "FDR P-value ",ifelse(plotdata$p.FDR[i]<0.001,"","= "), 
                                         paste.p(plotdata$p.FDR[i] )), 
                          hjust = "left", size = 2.8)
        yposadd <- yposadd + ifelse(zoom == zoomx, 0.02, 0.12)
      }
    } 
    
    adjustType <- case_when(adjustGroup=="tier2Type1" ~ "Tier 2: Set 1",
                            adjustGroup=="tier2Type3" ~ "Tier 2: Set 3",
                            adjustGroup=="tier2Type2and4and5" ~ "Tier 2: Sets 2, 4, and 5",
                            adjustGroup=="tier1Type1to3" ~ "Tier 1: Set 1 to 3",
                            adjustGroup=="tier1Type4to7" ~ "Tier 1: Set 4 to 7")
    
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


