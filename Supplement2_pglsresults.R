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

pglsfunc <- function() {
source("GLMM_subsetbirddataEfficient.R")

require(ape)
require(geiger)
require(nlme)
require(phytools)

pglsresults <- set.seed(10)
subset <- subsetbirddata(c("Polygyny","EPP","Syllrep"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Syllable.rep.final,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
pglsresults <- tableout
pglsresults <- cbind(rep("Syllrep",4), rownames(tableout), pglsresults)


subset <- subsetbirddata(c("Polygyny","EPP","Syllsong"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Syll.song.final,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Syllsong",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

subset <- subsetbirddata(c("Polygyny","EPP","Song"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Song.rep.final,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Song",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

subset <- subsetbirddata(c("Polygyny","EPP","Duration"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Duration.final,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Duration",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

subset <- subsetbirddata(c("Polygyny","EPP","Interval"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Interval.final,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Interval",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

subset <- subsetbirddata(c("Polygyny","EPP","Rate"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Song.rate,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Rate",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

subset <- subsetbirddata(c("Polygyny","EPP","Continuity"))
subdf <- subset$df
rownames(subdf) <- subdf$BirdtreeFormat
subtree <- subset$bothtree
Nspecies <- length(subdf$BirdtreeFormat)
pglsModel <- gls(log(Continuity,10) ~ Final.polygyny*Final.EPP, correlation = corBrownian(1,phy = subtree), data = subdf, method = "ML")
pglssumm <- summary(pglsModel)
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep("Continuity",4), rownames(tableout), tableout)
pglsresults <- rbind(pglsresults,temppglsresults)

write.csv(pglsresults, file = "pglsresults.csv")
}