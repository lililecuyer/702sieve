# Purpose: Estimation of hazard ratio-based VE by Tier 1 feature type 4-7 
#          Considers primary endpoints, and the PP cohort is the analysis cohort
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska andn Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on Github
# Input:   702 primary endpoints and mark variables specified above
# Output:  PDF files, each containing a single plot
# Author:  Li Li
# Date:    1/5, 2023

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
         eventInd = HIV24)
dat <- select(dat, all_of(c("tx", "eventTime", "eventInd", length_tier1, numSequon_tier1, charge_tier1, numCysteines_tier1,"cysteine.count.gp160.tier2")))

quantitativeMark_tier1 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier1",".csv"))) 
quantitativeMark_tier2 <- read.csv(file.path(repoDir, paste0("702sieve/tables/preScreen/varquantitativeMarksScreenedIn_tier2",".csv"))) 

# mark variable names
markName <- c(quantitativeMark_tier1$x, "cysteine.count.gp160.tier2")

# axis labels
xLab <- case_when(markName=="length.hypervariable.V1.tier1" ~ "Length of V1",
                  markName=="length.hypervariable.V2.tier1" ~ "Length of V2",
                  markName=="length.hypervariable.V4.tier1" ~ "Length of V4",
                  markName=="length.hypervariable.V5.tier1" ~ "Length of V5",
                  markName=="num.sequons.V1.tier1" ~ "Number of Sequons in V1",
                  markName=="num.sequons.V2.tier1" ~ "Number of Sequons in V2",
                  markName=="num.sequons.V3.tier1" ~ "Number of Sequons in V3",
                  markName=="num.sequons.V4.tier1" ~ "Number of Sequons in V4",
                  markName=="num.sequons.V5.tier1" ~ "Number of Sequons in V5",
                  markName=="num.sequons.CD4bs.Loop.D.tier1" ~ "Number of Sequons in CD4bs + Loop D",
                  markName=="charge.V1V2.tier1" ~ "Electrochemical Charge of V1V2",
                  markName=="charge.V3.tier1" ~ "Electrochemical Charge of V3",
                  markName=="charge.V4.tier1" ~ "Electrochemical Charge of V4",
                  markName=="charge.V5.tier1" ~ "Electrochemical Charge of V5",
                  markName == "cysteine.count.V1.tier1" ~ "Cysteine Count in V1",
                  markName == "cysteine.count.gp160.tier2" ~ "Cysteine Count in gp160")

# strings included in file names
markFileString <- gsub(".", "_", markName, fixed=TRUE)

# initialize vectors of p-values
pHRconstancy <- pHRunity <- NULL

# for each Hamming distance
for (i in 1:length(markName)){
  
  dat1 <- select(dat, tx, eventTime, eventInd, mark=markName[i])
  
  # convert mark values for non-primary endpoints to NA
  dat1$mark <- ifelse(dat1$eventInd==0, NA, dat1$mark)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  dat1 <- filter(dat1, !(eventInd==1 & is.na(mark)))
  
  # fit the mark-specific HR model
  if (markName[i] %in% c("length.v1v2")){  # plot a continuous VE curve
    markRng <- range(dat1$mark, na.rm=TRUE)
    markGrid <- seq(markRng[1], markRng[2], length.out=500)  
  } else {  # plot VE estimates as discrete points
    markGrid <- sort(unique(na.omit(dat1$mark)))
  }
  fit <- with(dat1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx))
  sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
  
  # 2-sided Wald test of {H0: VE(v)=0 for all v}
  pHRunity <- c(pHRunity, sfit$pWald.HRunity.2sided)
  
  # 2-sided Wald test of {H0: VE(v)=VE for all v} against {H1: VE(v) varies with v}
  pHRconstancy <- c(pHRconstancy, sfit$pWald.HRconstant.2sided)
  
  pdf(file.path(figureDir, paste0("702_sievePH_VE_", markFileString[i], "_PP",".pdf")), width=0.9*7, height=0.9*6.3)
  # plot(sfit, 
  #      mark=dat1$mark, 
  #      tx=dat1$tx, 
  #      xlim=NULL,
  #      ylim=c(-0.4, 1),
  #      xtickAt=NULL,
  #      xtickLab=NULL,
  #      ytickAt=seq(-0.4, 1, by=0.2),
  #      ytickLab=seq(-40, 100, by=20),
  #      xlab=xLab[i],
  #      ylab="Vaccine Efficacy (%)             ",
  #      txLab=c("Placebo", "Vaccine"),
  #      annotateP=TRUE)
  
  xLim <- range(dat1$mark, na.rm=TRUE)
  
  plotHeights <- c(0.32, 0.68)
  ggsave.width <- 0.7 * 7.3
  ggsave.height <- 0.7 * 6.5
  title <- NULL
  
  p.df1 <- read_csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", 
                                  "WestfallYoungAdjPvalues_tier2Type2and4and5.csv"), show_col_types=FALSE)
  p.df2 <- read_csv(file=file.path(repoDir,"702sieve/tables/westfallYoung", 
                                   "WestfallYoungAdjPvalues_tier1Type4to7.csv"), show_col_types=FALSE)
  p.df <- rbind(p.df1, p.df2)
  p.df <- filter(p.df, mark%in%markName)
  
  p <- as.numeric(p.df[p.df$mark==markName[i], c("p.unadj", "p.FWER", "p.FDR")])
  if(is.na(p[1])){
    unajustedP = sfit$pWald.HRconstant.2sided
    fmt.p <- ifelse(unajustedP<0.001, "< 0.001", paste0("= ", format(unajustedP, digits=2, nsmall=2)))
    subtitle <- paste0("Two-Sided Unadjusted Sieve P ", fmt.p)
  }else{
    fmt.p <-  sapply(p, function(x){ ifelse(x<0.001, "< 0.001", paste0("= ", format(x, digits=2, nsmall=2))) })
    subtitle <- paste0("Two-Sided Unadjusted Sieve P ", fmt.p[1], "\nFWER P ", fmt.p[2], ", FDR P ", fmt.p[3])
  }
  ylimL <- min(c(sfit$te$TE, -0.4))
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

pVals <- tibble(markName, pHRunity, pHRconstancy)
save(pVals, file=file.path(tableDir, paste0("702_sievePH_pValues_otherQuantMarks", ".RData")))

# test validity of the method's assumptions ------------------------------------------------------------------

# dat2 <- dat %>%
#   select(tx, eventTime, eventInd, mark=dist[1]) %>%
#   filter(!(eventInd==1 & is.na(mark)))
# dat2 <- as.data.frame(dat2)
# 
# # test validity of the assumption that T and V are conditionally independent given Z
# testIndepTimeMark(subset(dat2, tx==0, select=c("eventTime", "eventInd", "mark")), iter=1000)
# 0.552
# testIndepTimeMark(subset(dat2, tx==1, select=c("eventTime", "eventInd", "mark")), iter=1000)
# 0.692
# # conclusion: we do not reject validity of the conditional independence assumption
# 
# # test validity of the specified mark density ratio model
# with(dat2, testDensRatioGOF(eventInd, mark, tx))
# 0.758
# # conclusion: we do not reject validity of the mark density ratio model
