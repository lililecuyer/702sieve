datFile <- "HVTN_702_sieve_data_final_v3.csv"
dat <- read_csv(file.path(dataDir, datFile))
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
numCysteines_tier2 <- colnames(dat)[grepl("cysteine",colnames(dat)) & grepl("tier2",colnames(dat))]


