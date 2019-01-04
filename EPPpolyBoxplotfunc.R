########
#Coded by Kate T. Snyder
#Last Modified 11-28-2018
#Built using RStudio Version 1.1.456
#R Version 3.5.1
#
#ape_5.2  phytools_0.6-60   maps_3.3.0   nlme_3.1-137
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#R.utils_2.6.0   nortest_1.0-4   MCMCglmm_2.26
#BayesTraitsV2
########
########
source("Supplement_subsetbirddata.R")


EPPpolyfunc <- function() {
pdf(file = "EPPPolyALLboxplotsOnePdf.pdf", height = 5, width = 4)
par(mfrow = c(2,2))
par(mar=c(4,3,1,1))
pairwiseresults <- set.seed(10)
for (i in 1:7) {
SongParams <- c("Syllrep", "Syllsong", "Song", "Duration","Interval","Rate","Continuity")
songparam <- SongParams[i]
subsetoutput <- subsetbirddata("EPP", SongParam = songparam, withsuboscines = TRUE, newdata = FALSE)

matingtree <- subsetoutput$ditree
EPPdf <- subsetoutput$df
songcol <- subsetoutput$songcol

EPPpolydf <- EPPdf[which(!is.na(EPPdf$Final.polygyny)),]

songdatavec <- EPPpolydf[,songcol]
names(songdatavec) <- EPPpolydf$BirdtreeFormat
removespecies <- EPPdf$BirdtreeFormat[which(is.na(EPPdf$Final.polygyny))]
dropspeciestips <- which(matingtree$tip.label %in% removespecies)
EPPpolytree <- drop.tip(matingtree, dropspeciestips)

EPPdf <- EPPpolydf

withinEPPgroupvec <- set.seed(10)
withinEPPgroupvec[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 1)] <- "polyHEPP"
withinEPPgroupvec[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 0)] <- "polyLEPP"
withinEPPgroupvec[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 1)] <- "monoHEPP"
withinEPPgroupvec[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 0)] <- "monoLEPP"
names(withinEPPgroupvec) <- EPPdf$BirdtreeFormat

matecolsvec <- set.seed(10)
matecolsvec[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 1)] <- 4 #"deepskyblue"
matecolsvec[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 0)] <- 2 # "white"
matecolsvec[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 1)] <- 3 #"blue"
matecolsvec[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 0)] <- 1 #"black"
py <- c("black", "white", "blue", "deepskyblue")
names(matecolsvec) <- EPPdf$BirdtreeFormat
matecolsvec <- as.factor(matecolsvec)

polyHEPPdf <- EPPdf[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 1),]
NpolyHEPP <- length(polyHEPPdf$BirdtreeFormat)
polyLEPPdf <- EPPdf[which(EPPdf$Final.polygyny == 1 & EPPdf$Final.EPP == 0),]
NpolyLEPP <- length(polyLEPPdf$BirdtreeFormat)
monoHEPPdf <- EPPdf[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 1),]
NmonoHEPP <- length(monoHEPPdf$BirdtreeFormat)
monoLEPPdf <- EPPdf[which(EPPdf$Final.polygyny == 0 & EPPdf$Final.EPP == 0),]
NmonoLEPP <- length(monoLEPPdf$BirdtreeFormat)

logsong <- log(songdatavec,10)

df <- as.data.frame(cbind(names(withinEPPgroupvec),withinEPPgroupvec,as.numeric(logsong)), stringsAsFactors = FALSE)

boxplot(logsong ~ withinEPPgroupvec, data = df, xaxt = "n",  ylab="", cex.axis = 0.6)
title(ylab=paste("log10", songparam), cex.lab = 0.6)

groups <- c(paste0("Monogamy/HighEPP N=",NmonoHEPP),paste0("Monogamy/LowEPP N=", NmonoLEPP),paste0("Polygyny/HighEPP N=",NpolyHEPP),paste0("Polygyny/LowEPP N=",NpolyLEPP))
text(seq_along(withinEPPgroupvec[1:4]), par("usr")[3] - 0.05, labels = groups, srt = 45, adj = 1, xpd = TRUE, cex.lab = 0.5, cex=0.55, cex.main = 0.55)

set.seed(10)
phylanova <- phylANOVA(EPPpolytree,withinEPPgroupvec,logsong, nsim=2000)
phylanovap <- phylanova$Pf
print(songparam)
print(phylanova)

songparam4x <- rep(songparam,4)
temppairwise <- cbind(songparam4x,rownames(phylanova$Pt),phylanova$Pt)

pairwiseresults <- rbind(pairwiseresults,temppairwise)

title(paste("Log10", songparam, "& mating strategy \nphylANOVA p = ",phylanovap), cex.main = 0.6)

} # end for loop
dev.off()
write.csv(pairwiseresults, file = "EPPpolyPhylANOVApairwise.csv")
} #end function