########
#Coded by Kate T. Snyder
#Last Modified 11-28-2018
#Built using RStudio Version 1.0.136
#R Version 3.3.1
#
#ape_4.1  phytools_0.5-38   maps_3.1.0  btw_V1.0
#BayesTraitsV2 
########
########

#Must set .BayesTraitsPath to location of program BayesTraitsV2, which must be located in your working directory. E.g.:
#.BayesTraitsPath <- "~/Documents/BayesTraits/BayesTraitsV2"



require(phytools)
require(btw)
source("Supplement2_btwfunction_GeneTree.R")

cyclematesongsGene <- function(jackknife, csvsout, nsim = 5, plotbt = TRUE) {
  # MateParams <- c("Polygyny","EPP")
  # SongParams <- c("Syllrep","Syllsong","Song","Duration","Interval","Rate","Continuity")
  for (n in 1:length(MateParams)) {
    for (m in 1:length(SongParams)) {
      MateParam = MateParams[n]
      SongParam = SongParams[m]
      btwfunction(MateParam,SongParam,plot = plotbt,jackknife = jackknife,csvsout = csvsout,nsim=nsim,newtreefile = "matezillaGeneTreeHack.nex")
      #    transitionplots(MateParams[n],SongParams[n],csvfiles[n]) 
    } #end for m in 1:songparams
  } #end for n in 1:MateParams
} #end function BTfunc
