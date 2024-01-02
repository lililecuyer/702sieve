# Purpose: Estimation of hazard ratio-based VE by each Tier 2 Hamming distance (in SAP Table 5); 
#          Considers primary endpoints, and the PP cohort is the analysis cohort
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska andn Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on Github
# Input:   702 primary endpoints and Hamming distances
# Output:  PDF files, each containing a single plot
# Author:  Li Li
# Date:    Jan 1, 2024

rm(list=ls(all=TRUE))
here::i_am("702sieve/rootDir.R")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "702sieve/code/sievePH")
outputDir <- file.path(repoDir, "702sieve/code/sievePH/output")
figureDir <- file.path(repoDir, "702sieve/figures/sievePH")
tableDir <- file.path(repoDir, "702sieve/tables/sievePH")

library(tidyverse)
library(sievePH)
source(file.path(codeDir, "ggplot.summary.sievePH.R"))
source(file.path(repoDir, "702sieve/code/common.R"))

dat <- dat %>%
  mutate(tx=as.numeric(TRT01P=="T1"), 
         eventTime = HIV24fu, 
         eventInd = HIV24) %>% 
  select(tx, eventTime, eventInd, starts_with("hdist"))


# Hamming distance variable names
dist <- grep("hdist", colnames(dat), value=TRUE)
# axis labels
xLab <- case_when(dist =="hdist.zspace.1086.gp120.tier2" ~ "PC-weighted Hamming Distance to C.1086",
                  dist =="hdist.zspace.TV1.gp120.tier2" ~ "PC-weighted Hamming Distance to C.TV1",
                  dist=="hdist.zspace.96ZM651.V2.linear.epitope.tier2" ~ "V2 Linear Epitope PC-weighted Hamming Distance to C.96ZM651",
                  dist=="hdist.zspace.96ZM651.CD4.binding.site.tier2" ~ "CD4 Contact Epitopes PC-weighted Hamming Distance to C.96ZM651",
                  dist=="hdist.zspace.96ZM651.V5.tier2" ~ "V5 PC-weighted Hamming Distance to C.96ZM651",
                  dist=="hdist.zspace.96ZM651.V3.tier2" ~ "V3 PC-weighted Hamming Distance to C.96ZM651",
                  dist=="hdist.zspace.96ZM651.tier2" ~ "PC-weighted Hamming Distance to C.96ZM651")

# strings included in file names
markFileString <- gsub(".", "_", dist, fixed=TRUE)

# initialize vectors of p-values
pHRconstancy <- pHRunity <- NULL

format.p <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "= 1"}
      else{pp[i] <-paste0(" = ",as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits)))
      if(pp[i]== " = 1.00") {pp[i]= " = 1"}}
    }
  }
  return (pp)
}

p.df <- read_csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", 
                                "WestfallYoungAdjPvalues_tier2Type2and4and5.csv"), show_col_types=FALSE)
p.df <- filter(p.df, mark%in%dist)

# for each Hamming distance
for (i in 1:length(dist)){
  
  dat1 <- select(dat, tx, eventTime, eventInd, mark=dist[i])
  
  # convert mark values for non-primary endpoints to NA
  dat1$mark <- ifelse(dat1$eventInd==0, NA, dat1$mark)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  dat1 <- filter(dat1, !(eventInd==1 & is.na(mark)))
  
  # fit the mark-specific HR model
  markRng <- range(dat1$mark, na.rm=TRUE)
  markGrid <- seq(markRng[1], markRng[2], length.out=200)
  fit <- with(dat1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx))
  sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
  
  # 2-sided Wald test of {H0: VE(v)=0 for all v}
  pHRunity <- c(pHRunity, sfit$pWald.HRunity.2sided)
  
  # Double 1-sided Wald test of {H0: VE(v)=VE for all v} against {H1: VE(v) varies with v}
  
  pHRconstancy <- c(pHRconstancy, min(2*sfit$pWald.HRconstant.1sided,1))
  
  pdf(file.path(figureDir, paste0("702_sievePH_VE_", markFileString[i], "_PP",".pdf")), width=0.9*7, height=0.9*6.3)
  xLim <- range(dat1$mark, na.rm=TRUE)
  plotHeights <- c(0.32, 0.68)
  ggsave.width <- 0.7 * 7.3
  ggsave.height <- 0.7 * 6.5
  title <- NULL
  
  
  p <- as.numeric(p.df[p.df$mark==dist[i], c("p.unadj", "p.FWER", "p.FDR")])
  if(is.na(p[1])){
    unajustedP = min(2*sfit$pWald.HRconstant.1sided,1)
    fmt.p <- format.p(unajustedP)
    subtitle <- paste0("Double One-Sided Unadjusted Sieve P ", fmt.p)
    
  }else{
    fmt.p <-  sapply(p, function(x){ format.p(x)})
    subtitle <- paste0("Double One-Sided Unadjusted Sieve P ", fmt.p[1], "\nFWER P ", fmt.p[2], ", FDR P ", fmt.p[3])
  }
  ylimL <- -0.4
  
  p <- ggplot(sfit, 
              mark=dat1$mark, 
              tx=dat1$tx, 
              xlim=xLim,
              ylim=c(ylimL, 1),
              xtickAt=NULL,
              xtickLab=NULL,
              ytickAt=seq(-1, 1, by=0.2),
              ytickLab=seq(-100, 100, by=20),
              xlab=xLab[i],
              ylab="Vaccine Efficacy (%)",
              axisLabSize=15.5,
              legendLabSize=12,
              txLab=c("Placebo", "Vaccine"),
              jitterFactor=0.1,
              title=title,
              subtitle=subtitle,
              subtitleSize=13,
              estLineSize=1.8,
              ciLineSize=1.4,
              pointSize=2.1,
              plotHeights=plotHeights)
  
  
  print(p)
  # mtext(paste0(doseTitleString[dose], " by ", markTitleString[mark], ": ", trialTitleString[trial]), side=3, font=2, line=0.5, cex=1.4)
  dev.off()
}

pVals <- tibble(dist, pHRunity, pHRconstancy)
save(pVals, file=file.path(tableDir, "702_sievePH_pValues_hdist.RData"))

# test validity of the method's assumptions ------------------------------------------------------------------

#  dat2 <- dat %>%
#    select(tx, eventTime, eventInd, mark=dist[1]) %>%
#    filter(!(eventInd==1 & is.na(mark)))
#  dat2 <- as.data.frame(dat2)
# 
# # test validity of the assumption that T and V are conditionally independent given Z
#  testIndepTimeMark(subset(dat2, tx==0, select=c("eventTime", "eventInd", "mark")), iter=1000)
# # 0.176
#  testIndepTimeMark(subset(dat2, tx==1, select=c("eventTime", "eventInd", "mark")), iter=1000)
# # 0.571
# # conclusion: we do not reject validity of the conditional independence assumption
# #
# # test validity of the specified mark density ratio model
#  with(dat2, testDensRatioGOF(eventInd, mark, tx))
# # 0.539
# # conclusion: we do not reject validity of the mark density ratio model
