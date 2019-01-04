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
## 
## Outputs list of 7 items:
## bothtree - subsetted tree
## ditree (equal to bothtree)
## df - subsetted dataframe
## matevec - species-named vector of mating data values
## songcontdata - species-named vector of song data values
## matecol - column name for mating data
## songcol - column name for song data
subsetbirddata <- function(MateParam=c("Polygyny","EPP","OC","none"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity","OC","none"), withsuboscines = FALSE, newdata = FALSE, newtree = FALSE, islog = TRUE, minmax = NA) {
require(phytools)
require(ape)
require(maps)
  
output <- list()

if (newdata == FALSE) {
  alldatadfos <- as.data.frame(read.csv("SupplementalData_R.csv"))
} else {alldatadfos <- as.data.frame(read.csv(newdata))}

if (newtree == FALSE) {
  matezillatree <- read.nexus("matezillatreeHack.nex")
} else {matezillatree <- read.nexus(newtree)}

if (withsuboscines == FALSE) {
  suboscines <- c("Synallaxis_erythrothorax", "Synallaxis_brachyura", "Synallaxis_zimmeri", "Cranioleuca_erythrops", "Pseudoseisura_cristata", "Pseudoseisura_lophotes", "Phacellodomus_striaticollis", "Phacellodomus_maculipectus", "Leptasthenura_striolata", "Leptasthenura_yanacensis", "Cinclodes_fuscus", "Cinclodes_antarcticus", "Cinclodes_comechingonus", "Cinclodes_albiventris", "Furnarius_rufus", "Automolus_ochrolaemus", "Thripadectes_rufobrunneus", "Automolus_rubiginosus", "Xenops_minutus", "Sclerurus_guatemalensis", "Sclerurus_caudacutus", "Conopophaga_lineata", "Cercomacra_tyrannina", "Gymnopithys_leucaspis", "Empidonax_traillii", "Empidonax_alnorum", "Empidonax_minimus", "Empidonax_difficilis", "Empidonax_flaviventris", "Empidonax_virescens", "Contopus_sordidulus", "Sayornis_nigricans", "Sayornis_phoebe", "Pyrocephalus_rubinus", "Myiarchus_crinitus", "Myiarchus_cinerascens", "Tyrannus_tyrannus", "Elaenia_chiriquensis", "Oxyruncus_cristatus")
  suboscinetips <- match(suboscines,matezillatree$tip.label)
  oscinetree <- drop.tip(matezillatree, tip=suboscinetips)
  matezillatree <- oscinetree
  subosrows <- alldatadfos$BirdtreeFormat %in% suboscines
  alldatadf <- alldatadfos[which(!subosrows),]
 } else if (withsuboscines == TRUE) {
   alldatadf <- alldatadfos
 }
alldatadf <- alldatadf[which(alldatadf$BirdtreeFormat %in% matezillatree$tip.label),]  #two species are in dataset that aren't in overall tree, need to get them out of the dataset

#Go through each mating system 
if (MateParam == "Polygyny") {
  matedf <- alldatadf[!is.na(alldatadf$Final.polygyny),]
  havedatavec <- matezillatree$tip.label %in% as.character(matedf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  matingtree <- drop.tip(matezillatree, tip = dropformating)
  matecol <- "Final.polygyny"
} else if (MateParam == "EPP") {
  matedf <- alldatadf[!is.na(alldatadf$Final.EPP),] 
  matecol <- "Final.EPP"
  havedatavec <- matezillatree$tip.label %in% as.character(matedf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  matingtree <- drop.tip(matezillatree, tip = dropformating)
} else if (MateParam=="OC") {
  matedf <- alldatadf[!is.na(alldatadf$O.C),] 
  matecol <- "O.C"
  havedatavec <- matezillatree$tip.label %in% as.character(matedf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  matingtree <- drop.tip(matezillatree, tip = dropformating)
} else if (MateParam=="none") {
  matedf <- alldatadf
  matingtree <- matezillatree
  matecol <- "None"
} else {
  return(print("Not a valid value for MateParam"))}


#Go through all song metrics to subset data and tree
if (SongParam == "Syllrep") {
  songdf <- matedf[!is.na(matedf$Syllable.rep.final),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Syllable.rep.final"
  
} else if (SongParam == "Syllsong") {
  songdf <- matedf[!is.na(matedf$Syll.song.final),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Syll.song.final"
  
} else if (SongParam == "Song") {
  songdf <- matedf[!is.na(matedf$Song.rep.final),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Song.rep.final"
  
} else if (SongParam == "Interval") {
  songdf <- matedf[!is.na(matedf$Interval.final),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Interval.final"
  
} else if (SongParam == "Duration") {
  songdf <- matedf[!is.na(matedf$Duration.final),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Duration.final"
  
} else if (SongParam == "Continuity") {
  songdf <- matedf[!is.na(matedf$Continuity),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Continuity"
  
} else if (SongParam == "Rate") {
  songdf <- matedf[!is.na(matedf$Song.rate),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "Song.rate"
} else if (SongParam == "OC") {
  songdf <- matedf[!is.na(matedf$O.C),]
  havedatavec <- matingtree$tip.label %in% as.character(songdf$BirdtreeFormat)
  matingtips <- which(havedatavec == TRUE)
  dropformating <- which(havedatavec == FALSE)
  bothtree <- drop.tip(matingtree, tip = dropformating)
  songcol <- "O.C"
} else if (SongParam == "none") {
  songdf <- matedf
  songcol <- NULL
  bothtree <- matingtree
} else {return(print("Not a valid value for SongParam"))}


if (MateParam != "none") {
  matevec <- as.character(songdf[,matecol])
  names(matevec) <- songdf$BirdtreeFormat
} else {
  matevec <- "none"
}

#Extract data vectors for song traits
if (is.na(minmax)) {
  if (SongParam == "Syllrep") {
    songcol <- "Syllable.rep.final"
  } else if (SongParam == "Syllsong") {
    songcol <- "Syll.song.final"
  } else if (SongParam == "Song") {
    songcol <- "Song.rep.final"
  } else if (SongParam == "Interval") {
    songcol <- "Interval.final"
  } else if (SongParam == "Duration") {
    songcol <- "Duration.final"
  } else if (SongParam == "Rate") {
    songcol <- "Song.rate"
  } else if (SongParam == "Continuity") {
    songcol <- "Continuity"
  } else if (SongParam == "OC") {
    songcol <- "O.C"
  } else if (SongParam == "none") {
    songdatavec <- matevec
  } 
  
} else if (minmax == "min") {
  if (SongParam == "Syllrep") {
    songcol <- "Syllable.rep.min"
  } else if (SongParam == "Syllsong") {
    songcol <- "Syll.song.min"
  } else if (SongParam == "Song") {
    songcol <- "Song.rep.min"
  } else if (SongParam == "Interval") {
    songcol <- "Interval.min"
  } else if (SongParam == "Duration") {
    songcol <- "Duration.min"
  } else if (SongParam == "Rate") {
    songcol <- "Song.rate.min"
  } else if (SongParam == "Continuity") {
    songcol <- "Continuity.min"
  } else if (SongParam == "none") {
    songdatavec <- matevec
  } 
} else if (minmax == "max") {
  if (SongParam == "Syllrep") {
    songcol <- "Syllable.rep.max"
  } else if (SongParam == "Syllsong") {
    songcol <- "Syll.song.max"
  } else if (SongParam == "Song") {
    songcol <- "Song.rep.max"
  } else if (SongParam == "Interval") {
    songcol <- "Interval.max"
  } else if (SongParam == "Duration") {
    songcol <- "Duration.max"
  } else if (SongParam == "Continuity") {
    songcol <- "Continuity.max"
  } else if (SongParam == "Rate") {
    songcol <- "Song.rate.max"  
  } else if (SongParam == "none") {
    songdatavec <- matevec
  } 
}


ditree <- multi2di(bothtree)
ditree$edge.length[ditree$edge.length == 0] <- 0.0000000000000000001

output$ditree <- ditree
output$bothtree <- bothtree
output$matecol <- matecol
output$matevec <- matevec
output$df <- songdf

if (SongParam != "none") {
songcontvec <- as.numeric(songdf[,songcol])

if (islog == TRUE) {
  songcontvec <- log(songcontvec,10)
}
names(songcontvec) <- songdf$BirdtreeFormat
output$songcol <- songcol
output$songcontvec <- songcontvec
} #end if SongParam != "none"

return(output)
}


