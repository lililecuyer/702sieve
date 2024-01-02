#data summary
data.summary <- function(markData, var){

  fit.data.cases <- subset(markData,eventInd==1)
  
  table.cases <- table(fit.data.cases$mark, fit.data.cases$tx)
  #the order is 0-1 for mark and 0-1 for armLabel
  #number of cases per person-years
  vaccinePersonYears <- sum(markData$eventTime[markData$tx==1])/365.25
  placeboPersonYears <- sum(markData$eventTime[markData$tx==0])/365.25
  
  vaccineIncRate <- round(table.cases[,"1"]/vaccinePersonYears*100,1)
  placeboIncRate <- round(table.cases[,"0"]/placeboPersonYears*100,1)
  
  vaccineCases <- table.cases[,"1"]
  placeboCases <- table.cases[,"0"]
  
  return(list(vaccineCases = vaccineCases, placeboCases = placeboCases, 
              vaccineIncRate = vaccineIncRate, 
              placeboIncRate = placeboIncRate))
}

#Format VE estimate or the CI of VE to show percentages with one digit after zero
format.VE <- function(VE){

  if(as.numeric(VE) < (-10)){return("-Inf")}
  else{
    return(as.character(format(round(as.numeric(VE)*100,1),nsmall=1)))
  }}

#' Convert a p-value to a string and if it is less than 0.001, return '<0.001'
#' @param p An p-value
#' @example paste.p(0.1)
format.p <- function(p){if(as.numeric(p)<0.001){return("< 0.001")}else{return(as.character(format(as.numeric(p),digits=2)))}}

#' locate the position of a pattern in a vector of string
loc <- function(stringVec, pattern){
  seq = 1:length(stringVec)
  return(seq[stringVec %in% pattern])
}


#VE summary table with adjusted p values
table.VE.match.w.adjP <- function(markResultList, marks, WestfallYoungAdjPvalues, forestplot = FALSE, Insert = "Primary"){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(), feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                 P = character(), diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      Qvalue = format.p(fdr)
      if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
      
     
      VE <- add_row(VE, mark = mark, position = pos, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p(unadjp), FWER = format.p(holm), Qvalue = Qvalue)
      
      VE <- add_row(VE, mark = mark, position = "", feature  = paste0(Insert, " Match"), cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                   
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                    P = format.p(VEtable$p[1]), diffP = "", FWER = "", Qvalue = "")
      VE <- add_row(VE,mark = mark,  position = "", feature = paste0(Insert, " Mismatch"), cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                    P = format.p(VEtable$p[2]), diffP = "", FWER = "", Qvalue = "")
    }
  }else{
    VE <- tibble(position = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      sig.fdr <- ifelse(fdr<=0.2 & unadjp <=0.05, 1, 0) 
      sig.holm <- ifelse(holm<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.fdr, paste(annotatedPosition,"$^{\\S}$",sep=""), annotatedPosition)
      annotatedPosition = ifelse(sig.holm, paste(annotatedPosition,"$^{\\dag}$",sep=""), annotatedPosition)
      
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = format.p(fdr)
      )
    }
  }
  
  return(VE) 
}


#VE summary table without adjusted p values
table.VE.match.wo.adjP <- function(markResultList, marks, forestplot = FALSE, Insert = ""){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(), feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                 P = character(), diffP = character())
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
     
      
      VE <- add_row(VE, mark = mark, position = pos, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p(unadjp))
      
      VE <- add_row(VE, mark = mark, position = "", feature  = paste0(Insert, " Match"), cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                    P = format.p(VEtable$p[1]), diffP = "")
      VE <- add_row(VE,mark = mark,  position = "", feature = paste0(Insert, " Mismatch"), cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                    P = format.p(VEtable$p[2]), diffP = "")
    }
  }else{
    VE <- tibble(position = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character())
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP)
      )
    }
  }
  
  return(VE) 
}

#VE summary table for residue present/absent
table.VE.residue <- function(markResultList, marks, WestfallYoungAdjPvalues, forestplot = FALSE){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(),residue = character(), cases = character(), VE = character(), mean = numeric(),
                 ll = numeric(), ul = numeric(), P = character(), 
                 diffP = character(), FWER = character(), Qvalue = character() )
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      y = strsplit(mark,split="\\.")[[1]]
      pos <- y[2]
      residue <- y[4]
      unadjp <- markResultList[[mark]]$diffP
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      Qvalue = format.p(fdr)
      if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
      
      VE <- add_row(VE, 
                    mark = mark,
                    position = pos,
                    residue = "",
                    cases = "",
                    VE = "",
                    mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = Qvalue)
      
      VE <- add_row(VE, 
                    mark = mark,
                    position = "",
                    residue = residue,
                    cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,P = format.p(VEtable$p[1]), 
                    diffP = "", 
                    FWER = "",Qvalue = "")
      VE <- add_row(VE, 
                    mark = mark,
                    position = "",
                    residue = paste0 ("Not ", residue),
                    cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,P = format.p(VEtable$p[2]), 
                    diffP = "", 
                    FWER = "", Qvalue = "")
  }}else{
    VE <- tibble(position = character(), residue = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character(), FWER = character(), Qvalue = character() )
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      y = strsplit(mark,split="\\.")[[1]]
      pos <- y[2]
      residue <- y[4]
      
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      sig.fdr <- ifelse(fdr<=0.2 & unadjp <=0.05, 1, 0) 
      sig.holm <- ifelse(holm<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.fdr, paste(annotatedPosition,"$^{\\S}$",sep=""), annotatedPosition)
      annotatedPosition = ifelse(sig.holm, paste(annotatedPosition,"$^{\\dag}$",sep=""), annotatedPosition)
      
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    residue = residue,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = format.p(fdr)
      )
    }
  }
  
  return(VE) 
}

#VE summary table for sequon present/absent
table.VE.sequon <- function(markResultList, marks, WestfallYoungAdjPvalues, forestplot = FALSE){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(), feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                 P = character(), diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      Qvalue = format.p(fdr)
      if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
      
      
      VE <- add_row(VE, mark = mark, position = pos, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p(unadjp), FWER = format.p(holm), Qvalue = Qvalue)
      
      VE <- add_row(VE, mark = mark, position = "", feature  = paste0("PNG Present"), cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                    P = format.p(VEtable$p[1]), diffP = "", FWER = "", Qvalue = "")
      VE <- add_row(VE,mark = mark,  position = "", feature = paste0("PNG Absent"), cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                    P = format.p(VEtable$p[2]), diffP = "", FWER = "", Qvalue = "")
    }
  }else{
    VE <- tibble(position = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      sig.fdr <- ifelse(fdr<=0.2 & unadjp <=0.05, 1, 0) 
      sig.holm <- ifelse(holm<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.fdr, paste(annotatedPosition,"$^{\\S}$",sep=""), annotatedPosition)
      annotatedPosition = ifelse(sig.holm, paste(annotatedPosition,"$^{\\dag}$",sep=""), annotatedPosition)
      
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = format.p(fdr)
      )
    }
  }
  
  return(VE) 
}


#'plot the unadjusted and adjusted P values 
#'subType, addDesc = "", tier = 1, 
aggpvaluesPlot <- function(pvalues, marks, adj = FALSE, topText = NULL, bottomText = NULL, sizeText, panelHeight = 2.2,
                           yLab = "Differential VE 2-Sided P-value", 
                           figureDir = figureDir, figWidth = 8.5, figHeight = 8, fileName){
  library(ggpubr)
  library(RColorBrewer)
  library(gridExtra)
  

  orderedpvalues <- pvalues
  rownames(orderedpvalues) <- pvalues$X
  pvalues <- orderedpvalues[marks,]
  pvalues$abbrevMarks <- laply(pvalues$X, function(list){strsplit(list,split="\\.")[[1]][2]})
  
  numPerPanel = ifelse(length(pvalues$X)<75, length(pvalues$X), 75)
  #browser()
  seqs = seq(1,length(pvalues$abbrevMarks),numPerPanel)
  seqNum = c(seqs, length(pvalues$abbrevMarks))
  p <- list()
  for(i in 1:length(seqs)){
    
    begin = seqNum[i]
    end = ifelse(i == length(seqs), seqNum[i+1], seqNum[i+1]-1)
    subpvalues <-  pvalues[begin:end,]
    
    
    pVal = subpvalues$p.unadj
    pVal[pVal <= 0.01] = 0.007
    subpvalues$pValInd <- as.numeric(subpvalues$p.unadj <= 0.05)
    if(adj){
      pValHolm = subpvalues$p.FWER
      pValFDR = subpvalues$p.FDR
      pValHolm[pValHolm <= 0.01] = 0.007
      pValFDR[pValFDR <= 0.01] = 0.007
      subpvalues$HolmInd <- as.numeric(subpvalues$p.FWER <= 0.05)
      subpvalues$FDRInd <- case_when(subpvalues$p.FDR <= 0.05 ~ 1,
                                     subpvalues$p.FDR > 0.05 & subpvalues$p.FDR <=0.2 ~ -1,
                                     subpvalues$p.FDR > 0.2 ~0)
      pvaluesTable <- data.frame(c(pVal, pValHolm, pValFDR),
                                 c(rep("Unadjusted P-value",length(pVal)), rep("FWER-adjusted P-value",length(pVal)), rep("FDR-adjusted P-value",length(pVal))))
      pvaluesTable$ind = c(subpvalues$pValInd, subpvalues$HolmInd, subpvalues$FDRInd )
      pvaluesTable$x <- c(1:length(pVal), 1:length(pVal), 1:length(pVal))
      colnames(pvaluesTable) <- c("p","type","ind","x")
    }else{
      pvaluesTable <- data.frame(c(pVal),c(rep("Unadjusted P-value",length(pVal))))
      pvaluesTable$ind = c(subpvalues$pValInd, )
      pvaluesTable$x <- 1:length(pVal)
      colnames(pvaluesTable) <- c("p","type","ind","x")
    }
    
    yticklabels <- expression("1", "0.2", "0.1", "0.05", ""<= 0.01)
    
    p[[i]] <- ggplot()+
      geom_segment(aes(x = x, y = rep(0,length(x)), 
                       xend = x, yend = -log10(p) ), data = pvaluesTable)+
      geom_point(aes(x = x, y = -log10(p), colour = factor(ind), shape = factor(type)),stroke = 1.3,data = pvaluesTable)+
      scale_x_continuous(labels = subpvalues$abbrevMarks, breaks = 1:length(pVal), lim = c(1,numPerPanel))+
      scale_y_continuous(breaks=-log10(c(1,0.2,0.1,0.05,0.01)), labels=yticklabels, lim = c(0,2.3))+
      geom_hline(yintercept = -log10(0.05) ,linetype=2,color = "gray40", size = 0.5)+
      geom_hline(yintercept = -log10(0.2) ,linetype=2,color = "gray40", size = 0.5)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95, size = 6.5, colour = "black"),
            plot.margin = margin(0, 0.2, 0, 0.2, "cm"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black", hjust = 1),
            legend.position="right",legend.direction="horizontal",
            legend.title=element_text(colour="black"),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(colour="black",size=8,face="bold"),
            legend.key.height=grid::unit(0.8,"cm"),
            legend.key.width=grid::unit(0.4,"cm"))+
      xlab("")+
      ylab("")+
      guides(color = "none")
  if(adj){
    p[[i]] <- p[[i]] + scale_colour_manual(values = c("1" = "red","0" = "black","-1" = "darkgoldenrod2"))+
      scale_shape_manual(name="",values = c("Unadjusted P-value" = 1, "FWER-adjusted P-value" = 4, "FDR-adjusted P-value" = 2))
      
  }else{
    p[[i]] <- p[[i]] + scale_colour_manual(values = c("1" = "red","0" = "black"))+
      scale_shape_manual(name="",values = c("Unadjusted P-value" = 1))
      
  }  
    
  }
  if(length(p)>0){
    
    multiPanelPlot <- ggpubr::ggarrange(plotlist=p, nrow=length(p), ncol=1, align="n", heights=panelHeight, 
                                        common.legend=TRUE, legend="top")
    multiPanelPlot <- annotate_figure(multiPanelPlot, left=text_grob(yLab, rot=90, size=sizeText),
                                      bottom=text_grob(bottomText, size=sizeText), 
                                      top = text_grob(topText, size = sizeText))
   
    ggpubr::ggexport(multiPanelPlot, filename= file.path(figureDir, fileName), width=figWidth, height=figHeight)
    
    
  }
  
}

           
format.p0 <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "1"}
      else{pp[i] <-as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits))
      if(pp[i] == "1.00"){pp[i]="1"}
      }
    }
  }
  return (pp)
}
format.p2 <- function(p, ndigits=2){
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

forestplotBinary <- function(forestTable, fileName, outputDir, width = 11, height1 = 0.25, height2 = 1){
  header.plot <- c("HXB2\nPosition", "Seq\nFeature", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE(%)(95% CI)", "Two-sided\nP-value", "P-value", "FWER\nP-value","Q-value")
  p <- forestWrap ( forestTable, c(5,6,7), header.plot, xlim = c(-30, 100), 
                    x_ticks_at = c(-20, 0, 20, 40, 60, 80), xlab = "VE (%) (95% CI)", 
                    ci_cols = rep(c("black","royalblue", "darkred"),dim(forestTable)[1]/3),
                    nrows_underline = NULL,
                    insertHeaderText = "Two-sided Differential VE",
                    insertHeaderTextCol = c(7:9))
  
  ggplot2::ggsave(filename = fileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = width, height = height1*(length(forestTable$mean)+height2),units = "in")  
}

