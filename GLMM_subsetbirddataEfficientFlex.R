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
require(R.utils)
require(phytools)
require(ape)
require(base)
require(mnormt)
require(plyr)
require(geiger)
require(btw)
require(nlme)
require(phytools)
require(nortest)

## One function to rule them all - subset bird data and generate associated trees
## Can subset by 3 groups
## 
## Outputs list of 7 items:
## bothtree - subsetted tree
## ditree (equal to bothtree)
## df - subsetted dataframe
## matevec - species-named vector of mating data values
## songcontdata - species-named vector of song data values
## matecol - column name for mating data
## songcol - column name for song data
subsetbirddata <- function(params = NA, withsuboscines = FALSE, newdata = FALSE, newtree = FALSE, islog = TRUE, minmax = NA) {
require(phytools)
require(ape)
require(maps)
  
output <- list()

if (newdata == FALSE) {
  alldatadfos <- as.data.frame(read.csv("SupplementalData_R.csv"))
} else {alldatadfos <- as.data.frame(read.csv(newdata))}

if (newtree == FALSE) {
  matezillatree <- read.nexus("matezillatreeHack.nex")
}

if (withsuboscines == FALSE) {
  suboscines <- c("Synallaxis_erythrothorax", "Synallaxis_brachyura", "Synallaxis_zimmeri", "Cranioleuca_erythrops", "Pseudoseisura_cristata", "Pseudoseisura_lophotes", "Phacellodomus_striaticollis", "Phacellodomus_maculipectus", "Leptasthenura_striolata", "Leptasthenura_yanacensis", "Cinclodes_fuscus", "Cinclodes_antarcticus", "Cinclodes_comechingonus", "Cinclodes_albiventris", "Furnarius_rufus", "Automolus_ochrolaemus", "Thripadectes_rufobrunneus", "Automolus_rubiginosus", "Xenops_minutus", "Sclerurus_guatemalensis", "Sclerurus_caudacutus", "Conopophaga_lineata", "Cercomacra_tyrannina", "Gymnopithys_leucaspis", "Empidonax_traillii", "Empidonax_alnorum", "Empidonax_minimus", "Empidonax_difficilis", "Empidonax_flaviventris", "Empidonax_virescens", "Contopus_sordidulus", "Sayornis_nigricans", "Sayornis_phoebe", "Pyrocephalus_rubinus", "Myiarchus_crinitus", "Myiarchus_cinerascens", "Tyrannus_tyrannus", "Elaenia_chiriquensis", "Oxyruncus_cristatus")
  suboscinetips <- match(suboscines,matezillatree$tip.label)
  oscinetree <- drop.tip(matezillatree, tip=suboscinetips)
  matingtree <- matezillatree <- oscinetree
  subosrows <- alldatadfos$BirdtreeFormat %in% suboscines
  matedf <- alldatadf <- alldatadfos[which(!subosrows),]
 } else if (withsuboscines == TRUE) {
   matedf <- alldatadf <- alldatadfos
 }
alldatadf <- alldatadf[which(alldatadf$BirdtreeFormat %in% matezillatree$tip.label),]  #two species are in dataset that aren't in overall tree, need to get them out of the dataset

#Go through each mating system 
matecolvec <- set.seed(10)
if ("Polygyny" %in% params) {
  matecolvec <- c(matecolvec,"Final.polygyny")
}
if ("EPP" %in% params) {
  matecolvec <- c(matecolvec,"Final.EPP")
} 
if ("OC" %in% params) {
  matecolvec <- c(matecolvec,"O.C")
} 
if ("Syllrep" %in% params) {
  matecolvec <- c(matecolvec,"Syllable.rep.final")
}
if ("Syllsong" %in% params) {
  matecolvec <- c(matecolvec,"Syll.song.final")
}
if ("Song" %in% params) {
  matecolvec <- c(matecolvec,"Song.rep.final")
}
if ("Duration" %in% params) {
  matecolvec <- c(matecolvec,"Duration.final")
}
if ("Interval" %in% params) {
  matecolvec <- c(matecolvec,"Interval.final")
}
if ("Rate" %in% params) {
  matecolvec <- c(matecolvec,"Song.rate")
}
if ("Continuity" %in% params) {
  matecolvec <- c(matecolvec,"Continuity")
}

if (length(params) == 1) {
  songdf <- matedf[which(!is.na(matedf[,matecolvec[1]])),]
} else if (length(params) == 2) {
  songdf <- matedf[which(!is.na(matedf[,matecolvec[1]]) & !is.na(matedf[,matecolvec[2]])),]
} else if (length(params) == 3) {
  songdf <- matedf[which(!is.na(matedf[,matecolvec[1]]) & !is.na(matedf[,matecolvec[2]]) & !is.na(matedf[,matecolvec[3]])),] 
  }
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)


ditree <- multi2di(bothtree)
ditree$edge.length[ditree$edge.length == 0] <- 0.0000000000000000001

output$matecolvec <- matecolvec
output$ditree <- ditree
output$bothtree <- bothtree
output$df <- songdf


return(output)
}


