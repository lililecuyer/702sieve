
#'produce a character string with screened sites 
cellPaste <- function(screenedIn, celllength = 65){

  positions <- screenedIn$position
  digits <- nchar(as.character(positions))
  totalDigits <- sum(digits)
  cumSumDigits <- cumsum(digits)
  
  digitSeq <- seq(1, totalDigits, celllength)
  if(max(digitSeq)< totalDigits){digitSeq <- c(digitSeq, totalDigits)}
  seq <- 1
  for(i in 2:length(digitSeq)){
    seq <- c(seq, seq(1, length(positions),1)[cumSumDigits == max(cumSumDigits[cumSumDigits<=digitSeq[i]])])
  }
  outputString <- paste(positions[1:(seq[2])],collapse=",")
  if(length(seq)>2){
    for(i in 2:(length(seq)-1)){
      outputString <- paste0(outputString, ",\\newline", paste(positions[(seq[i]+1):(seq[i+1])],collapse=","))
    }
  }
  if(max(seq) < length(positions)){
    outputString <- paste0(outputString, ",\\newline", paste(positions[(seq[i]+1):(length(positions))],collapse=","))
    
  }
  return(outputString)
}

#produce screen tables and heat maps for match/mismatch
featureScreen.MatchvsMismatch <- function(sieveData, featureVars, newFeatureName, cutoff = 4, tableDir, figureDir, fileName, plotTitle){

  featureData <- subset(dplyr::select(sieveData, all_of(c("subjid","TRT01P", "hiv1event",featureVars))), hiv1event == 1)
  colnames(featureData ) <- c("subjid","TRT01P", "hiv1event",newFeatureName)
  
  featureLongFormat <- tidyr::pivot_longer(featureData, cols = all_of(newFeatureName), 
                                           names_to = "position", 
                                           values_to = "match",
                                           values_drop_na = TRUE)#there are NAs in the sequence

  freqsTable <- plyr::ddply(featureLongFormat, .(TRT01P, position, match), 
                                  function(df){
                                    length(unique(df$subjid))
                                  })

  table1 <- tibble(TRT01P = freqsTable$TRT01P,
                   position = freqsTable$position,
                   match = freqsTable$match,
                   "nCases" = freqsTable$V1)
  
  table2 <- tibble("position" = character(), "match" = numeric(), "nTxPoolCases" = numeric())
  
  d_ply(table1, .(position), function(df){
    for(pre.k in c(0,1)){
      nTxPoolCases <- ifelse(is.null(subset(df, match == pre.k)), 0,sum(df$nCases[df$match == pre.k]))
      table2 <<- add_row(.data = table2, "position" = df$position[1], "match" = pre.k,
                         "nTxPoolCases" = nTxPoolCases)
    }
  })
  

  table2$ind <- as.factor(ifelse(table2$nTxPoolCases>=cutoff, paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$ind <- factor(table2$ind, levels = c(paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$plotY <- factor(table2$position, levels = rev(newFeatureName))
  
  table2 <- mutate(table2,matchString = case_when(
                       match == 1 ~ "Match",
                       match == 0 ~ "Mismatch",
                     ))
  
  table2$matchString <- factor(table2$matchString, levels = c("Match", "Mismatch"))
  
 
  #Screened in comparison table
  posScreenedIn <- subset(ddply(table2,.(position), function(df){sum(df$nTxPoolCases>=cutoff)==2}), V1 == TRUE)
  posScreenedIn$position <- sort(factor(posScreenedIn$position, levels = newFeatureName))
  
  
  
  table3 <- tibble("m" = numeric(),"Positions Screened in" = character())
  
  table3 <- add_row(.data = table3, "m" = length(posScreenedIn$position), 
                    "Positions Screened in" = cellPaste(posScreenedIn))
 
  
  write.csv(table3, file.path(tableDir, paste0(fileName, ".csv")), row.names = FALSE)
  #browser()
  x = data.frame(vars = featureVars, newName = newFeatureName)
  write.csv(x$vars[x$newName %in% posScreenedIn$position], file.path(tableDir, paste0("var", fileName, ".csv")), row.names = FALSE)
  
  
  #plot 
  loc <- newFeatureName
  if(length(loc)>=80){
    cutloc <- seq(1, length(loc), 80)
    
  }else{
    cutloc <- c(1, length(loc))
    
  }

  locGroup <- list()
  locGroupName <- NULL
  for(i in 1:(length(cutloc)-1)){
    if(i==1){
      locGroup[[i]] <- loc[1:(cutloc[i+1])]
      locGroupName[i] <- paste(loc[cutloc[1]], loc[cutloc[i+1]],sep="-")
    }else{
      locGroup[[i]] <- loc[(cutloc[i]+1):(cutloc[i+1])]
      locGroupName[i] <- paste(loc[cutloc[i]+1], loc[cutloc[i+1]],sep="-")
    }
    
  }
  names(locGroup) <- locGroupName
  for(i in 1:length(locGroupName)){
    subTable <- subset(table2, plotY %in% locGroup[[i]])
    pdf(file = file.path(figureDir, paste0("plot_",fileName,"_",i, ".pdf")), 
        width = 3.7, height = 10.5)
    textcol <- "black"
    labelsize <- 10
    p <- ggplot(subTable, aes(x = matchString, y = plotY, fill = ind))+
      geom_tile(colour="white")+
      guides(fill=guide_legend(title="No. of\ncases"))+
      scale_y_discrete(expand=c(0,0))+
      scale_fill_manual(values=c("#008080", "grey"))+
      theme_grey(base_size=10)+
      #facet_grid(cols = vars(studySiteF), labeller = label_value)+
      xlab("")+
      ylab("")+
      #coord_equal()+
      ggtitle(paste(plotTitle, locGroupName[i]))+
      theme(
        legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol, size = labelsize),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=labelsize,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.4,"cm"),
        #plot.margin = grid::unit(c(0,0,marginb,0 ), "cm"), #top, right, bottom, left
        plot.title = element_text(colour=textcol, hjust = 0.5, vjust = 0.2, size = 12),
        axis.text.x = element_text(colour=textcol, angle = 0, hjust = 0.5,vjust = 0, size = labelsize),
        axis.text.y = element_text(colour=textcol,size = labelsize, vjust = 0.2),
        strip.text.x = element_text(colour=textcol,size = 16))
    print(p)
    dev.off()
    
  }
}

#produce screen tables and heat maps for amino acid present/absent
featureScreen.presentVsAbsent <- function(sieveData, featureVars, newFeatureName, cutoff = 4, tableDir, figureDir, fileName, plotTitle){
  #browser()
  featureData <- subset(dplyr::select(sieveData, all_of(c("subjid","TRT01P", "hiv1event",featureVars))), hiv1event == 1)
  colnames(featureData ) <- c("subjid","TRT01P", "hiv1event",newFeatureName)
  
  featureLongFormat <- tidyr::pivot_longer(featureData, cols = all_of(newFeatureName), 
                                           names_to = "position", 
                                           values_to = "present",
                                           values_drop_na = TRUE)#there are NAs in the sequence
  
  freqsTable <- plyr::ddply(featureLongFormat, .(TRT01P, position, present), 
                            function(df){
                              length(unique(df$subjid))
                            })
  
  table1 <- tibble(TRT01P = freqsTable$TRT01P,
                   position = freqsTable$position,
                   present = freqsTable$present,
                   "nCases" = freqsTable$V1)
  
  table2 <- tibble("position" = character(), "present" = numeric(), "nTxPoolCases" = numeric())
  
  d_ply(table1, .(position), function(df){
    for(pre.k in c(0,1)){
      nTxPoolCases <- ifelse(is.null(subset(df, present == pre.k)), 0,sum(df$nCases[df$present == pre.k]))
      table2 <<- add_row(.data = table2, "position" = df$position[1], "present" = pre.k,
                         "nTxPoolCases" = nTxPoolCases)
    }
  })
  
  
  table2$ind <- as.factor(ifelse(table2$nTxPoolCases>=cutoff, paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$ind <- factor(table2$ind, levels = c(paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$plotY <- factor(table2$position, levels = rev(newFeatureName))
  
  table2 <- mutate(table2,matchString = case_when(
    present == 1 ~ "Present",
    present == 0 ~ "Absent",
  ))
  #browser()
  table2$matchString <- factor(table2$matchString, levels = c("Present", "Absent"))
  
  #Screened in comparison table
  posScreenedIn <- subset(ddply(table2,.(position), function(df){sum(df$nTxPoolCases>=cutoff)==2}), V1 == TRUE)
  posScreenedIn$position <- sort(factor(posScreenedIn$position, levels = newFeatureName))
  
  
  
  table3 <- tibble("m" = numeric(),"Positions Screened in" = character())
  
  table3 <- add_row(.data = table3,  "m" = length(posScreenedIn$position), 
                    "Positions Screened in" = cellPaste(posScreenedIn,celllength = 80))
  
 
  write.csv(table3, file.path(tableDir, paste0(fileName, ".csv")), row.names = FALSE)
  #browser()
  x = data.frame(vars = featureVars, newName = newFeatureName)
  write.csv(x$vars[x$newName %in% posScreenedIn$position], file.path(tableDir, paste0("var", fileName, ".csv")), row.names = FALSE)
  
  
  #plot 
  loc <- newFeatureName
  cutloc <- seq(1, length(loc), 40)
  
  locGroup <- list()
  locGroupName <- NULL
  for(i in 1:(length(cutloc)-1)){
    if(i==1){
      locGroup[[i]] <- loc[1:(cutloc[i+1])]
      locGroupName[i] <- paste(loc[cutloc[1]], loc[cutloc[i+1]],sep="-")
    }else{
      locGroup[[i]] <- loc[(cutloc[i]+1):(cutloc[i+1])]
      locGroupName[i] <- paste(loc[cutloc[i]+1], loc[cutloc[i+1]],sep="-")
    }
    
  }
  names(locGroup) <- locGroupName
  for(i in 1:length(locGroupName)){
    subTable <- subset(table2, plotY %in% locGroup[[i]])
    pdf(file = file.path(figureDir, paste0("plot_",fileName,"_",i, ".pdf")), 
        width = 3.8, height = 8)
    textcol <- "black"
      labelsize <- 12
      p <- ggplot(subTable, aes(x = matchString, y = plotY, fill = ind))+
        geom_tile(colour="white")+
        guides(fill=guide_legend(title="No. of\ncases"))+
        scale_y_discrete(expand=c(0,0))+
        scale_fill_manual(values=c("#008080", "grey"))+
        theme_grey(base_size=10)+
        #facet_grid(cols = vars(studySiteF), labeller = label_value)+
        xlab("")+
        ylab("")+
        #coord_equal()+
        ggtitle(paste(plotTitle, locGroupName[i]))+
        theme(
          legend.position="right",legend.direction="vertical",
          legend.title=element_text(colour=textcol, size = labelsize),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(colour=textcol,size=labelsize,face="bold"),
          legend.key.height=grid::unit(0.8,"cm"),
          legend.key.width=grid::unit(0.4,"cm"),
          #plot.margin = grid::unit(c(0,0,marginb,0 ), "cm"), #top, right, bottom, left
          plot.title = element_text(colour=textcol, hjust = 0.5, vjust = 0.2, size = 12),
          axis.text.x = element_text(colour=textcol, angle = 0, hjust = 0.5,vjust = 0, size = labelsize),
          axis.text.y = element_text(colour=textcol,size = labelsize, vjust = 0.2),
          strip.text.x = element_text(colour=textcol,size = 16))
      print(p)
      dev.off()
      
  }
}

