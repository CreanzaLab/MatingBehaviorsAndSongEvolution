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

matesongfunc <- function(MateParam=c("Polygyny","EPP","OC","none"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity","OC","none"), matemodel = "ARD", heattree=FALSE, brownie=FALSE, minmax = NA, nsim = 1000, newmatingtreefile = FALSE) {
  require(R.utils)
  require(phytools)
  require(ape)
  require(base)
  require(mnormt)
  require(plyr)
  require(geiger)
  
  source(file = "Supplement_subsetbirddataNewTree.R")
  source(file = "Supplement2_browniefunctionNewTree.R")

  
  output <- list()
  
  subsetout <- subsetbirddata(MateParam, SongParam, islog = TRUE, newdata = FALSE, newtree = newmatingtreefile, withsuboscines = TRUE) #log10
  matingtree <- subsetout$bothtree
  songdf <- subsetout$df
  matingdatavec <- subsetout$matevec
  matingtree <- subsetout$ditree
  songdatavec <- subsetout$songcontvec
  


  
  if (SongParam != "none") {
    set.seed(10)
    phylanova <- phylANOVA(matingtree,matingdatavec,songdatavec, nsim=nsim)
    output$phylanovatitle <- paste(MateParam,SongParam,minmax,"N =", length(songdatavec))
    output$phylanova <- phylanova
    
    if (brownie == TRUE) {
      browniefunction(MateParam,SongParam, matensim = nsim, phylanova = phylanova$Pf)
    } # end if brownie=T
    
    if (heattree == TRUE) {
        plotheattree(MateParam,SongParam)
    }  #end if plot=T 
    
  }  #end != "none"
  
}