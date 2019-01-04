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

#Supplemental Code: matesongfunc
#Performs phylanova, brownie, plots brownie rate distributions

#
# matesongfunc("Polygyny","Syllrep",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Syllsong",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Song",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Duration",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Interval",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Rate",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("Polygyny","Continuity",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Syllrep",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Syllsong",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Song",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Duration",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Interval",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Rate",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)
# matesongfunc("EPP","Continuity",brownie = TRUE,matensim=20, ntrees=100, minmax = NULL)


matesongfunc <- function(MateParam=c("Polygyny","EPP","OC","none"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity","OC","none"), matemodel = "ARD", heattree=FALSE, brownie=FALSE, minmax = NA, matensim = 1000, ntrees = 10) {
  require(R.utils)
  require(phytools)
  require(ape)
  require(base)
  require(mnormt)
  require(plyr)
  require(geiger)
  
  source(file = "Supplement_subsetbirddataNewTree.R")
  source(file = "Supplement_findQratesNewTree.R")
  source(file = "Supplement_browniefunctionNewTree.R")

  
  output <- list()
  
  phylanovavec <- set.seed(10)
  browniedata <- data.frame(MatePar=character(),SongPar=character(),Pval=numeric(),ERRate=numeric(),ERloglik=numeric(),ERace=numeric(),ARDRate0=numeric(),ARDRate1=numeric(),ARDloglik=numeric(),ARDace=numeric(),k2=numeric(),convergence=character(),simmapnumber=integer(),phylanovaP=numeric(),stringsAsFactors = FALSE)

  brownie100df <- data.frame(MatePar=character(),SongPar=character(),Pval=numeric(),ERRate=numeric(),ERloglik=numeric(),ERace=numeric(),ARDRate0=numeric(),ARDRate1=numeric(),ARDloglik=numeric(),ARDace=numeric(),k2=numeric(),convergence=character(),simmapnumber=integer(),phylanovaP=numeric(),treenum=numeric(),stringsAsFactors = FALSE)
  
  for (i in 1:ntrees) {
  
  newtree <- paste("samplematezillaHack",i,".nex",sep="")
  
  subsetout <- subsetbirddata(MateParam, SongParam, islog = TRUE, newdata = FALSE, newtree = newtree) #log10
  matingtree <- subsetout$bothtree
  songdf <- subsetout$df
  matingdatavec <- subsetout$matevec
  matingtree <- subsetout$ditree
  songdatavec <- subsetout$songcontvec


  if (SongParam != "none") {
    set.seed(10)
    phylanova <- phylANOVA(matingtree,matingdatavec,songdatavec, nsim=matensim*20)
    output$phylanovatitle <- paste(MateParam,SongParam,minmax,"N =", length(songdatavec))
    output$phylanova <- phylanova
    phylanovavec[i] <- phylanova$Pf
     
    if (brownie == TRUE) {
     # pdf(file = paste(MateParam,SongParam,"10treebrownieplots.pdf"), width = 8, height = 11)
     # par(mfrow = c(3,6))
      brownietempdf <- browniefunction(MateParam,SongParam, matensim = matensim, phylanova = phylanova$Pf, newtree = newtree)
      treenum <- rep(i,times = length(brownietempdf$Pval))
      brownietempdf <- cbind(brownietempdf,treenum)
     # dev.off()
    } # end if brownie=T
    brownie100df <- rbind(brownie100df,brownietempdf)
    
  }  #end != "none"
  
}  # end for i in 1:100 (trees)
  write.csv(brownie100df,file = paste(MateParam,SongParam,ntrees,"treesBrownie.csv",sep=""))
  
  
  ###Plot overall Brownie distributions
  
  brownied <- brownie100df
  browniedf <- brownied[brownied$convergence == "Optimization has converged.",]
  browniedf <- browniedf[!is.na(browniedf$convergence),]
  ARDratio0to1 <- browniedf$ARDRate0/browniedf$ARDRate1
  #browniedf <- cbind(browniedf,ARDratio0to1)
  all0over1df <- browniedf[which(ARDratio0to1 > 1),]
  all1over0df <- browniedf[which(ARDratio0to1 < 1),]
  
  ERloglikmean <- mean(browniedf$ERloglik)
  ARDloglikmean <- mean(browniedf$ARDloglik)
  P.chisqAll=pchisq(2*(ARDloglikmean-as.numeric(ERloglikmean)),1,lower.tail=FALSE)
  ERloglikmean0over1 <- mean(all0over1df$ERloglik)
  ARDloglikmean0over1 <- mean(all0over1df$ARDloglik)
  P.chisq0over1=pchisq(2*(ARDloglikmean0over1-as.numeric(ERloglikmean0over1)),1,lower.tail=FALSE)
  ERloglikmean1over0 <- mean(all1over0df$ERloglik)
  ARDloglikmean1over0 <- mean(all1over0df$ARDloglik)
  P.chisq1over0=pchisq(2*(ARDloglikmean1over0-as.numeric(ERloglikmean1over0)),1,lower.tail=FALSE)
  
 # pdf(file=paste("ALLbrownie",ntrees,"treetest",matensim, "sims", MateParam,SongParam,".pdf",sep=""),width = 5, height = 5)
  par(mfrow=c(1,2))
  if (MateParam == "Polygyny") {
    state0 <- "Monogamy"
    state1 <- "Polygyny"
  } else if (MateParam == "EPP") {
    state0 <- "Low EPP"
    state1 <- "High EPP"
  }
  D0 <- density(browniedf$ARDRate0)
  D1 <- density(browniedf$ARDRate1)
  nsimlab = matensim
#  par(mar = c(1,1,2,1))
  plot(D0,col="blue",
       xlim=c(min(c(D0$x,D1$x)),
              max(c(D0$x,D1$x))),
       ylim=c(min(c(D0$y,D1$y)),
              max(c(D0$y,D1$y))),main="All converged runs, multiple trees",
       cex.main = 0.75, 
       xlab=paste("Rate of",SongParam,"evolution"),
       ylab="Frequency")
  lines(D1, col="red")
#  title("bottomright",labels = paste("#trees = ", ntrees, ", #sims/tree = ", matensim))
  abline(v=browniedf$ERRate[1], lty = 2)
  if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] > browniedf$ERRate[1]) {
    legend("topleft",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
  } else if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] < browniedf$ERRate[1]) {
    legend("topright",legend = c(paste(state0),paste(state1), "Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2)) 
  }
    #plot brownie pvalue 
  D0 <- density(browniedf$Pval)
  sdev <- sd(browniedf$Pval)
  meanphy <- mean(browniedf$Pval)
  plot(D0,col="black",
       xlim=c(min(D0$x),
              max(D0$x)),
       ylim=c(min(D0$y),
              max(D0$y)),
       main=paste(MateParam, SongParam, "Brownie pvals", ntrees, "trees", ", # sims per tree = ", nsimlab, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4)), cex.main = 0.75, xlab="Pval" ,ylab="Frequency") 
  abline(v=0.05, col = "gray")
  dev.off()
  
  
  ###Plot overall phylanova distribution
  
  D0 <- density(phylanovavec)
  sdev <- sd(phylanovavec)
  meanphy <- mean(phylanovavec)
  nsimlab = matensim*20
  pdf(file = paste(MateParam,SongParam,ntrees,"treephylanovadist.pdf",sep=""))
  plot(D0,col="blue",
       xlim=c(min(D0$x),
              max(D0$x)),
       ylim=c(min(D0$y),
              max(D0$y)),
       main=paste(MateParam, SongParam, "PhylANOVA pvals", ntrees, "trees", ", # sims per tree = ", nsimlab, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4)), cex.main = 0.75, xlab="Pval" ,ylab="Frequency") 
  abline(v=0.05, col = "gray")
  dev.off()
} 