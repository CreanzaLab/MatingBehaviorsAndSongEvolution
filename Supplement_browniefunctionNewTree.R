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
browniefunction <- function(MateParam = MateParam,SongParam = SongParam, matensim = 500, phylanova = "not evaluated", minmax = NULL, newtree = newtree) {
  source(file = "Supplement_findQratesNewTree.R")
  print("Beginning simmap for brownie") 
  starttimebrownie <- Sys.time()
  
  
  
  Qoutput <- findQrates(MateParam = MateParam, SongParam = "none", withsuboscines = TRUE, newdata = FALSE, newtree = newtree)
  qrates <- Qoutput$qrates
  
  
  subsetoutput <- subsetbirddata(MateParam, SongParam, withsuboscines = TRUE, newdata = FALSE, newtree = newtree) 
  
  matingtree <- subsetoutput$ditree  
  matecol <- subsetoutput$matecol
  matingdatavec <- subsetoutput$matevec
  songdatavec <- subsetoutput$songcontvec
  
  simmappy <- make.simmap(matingtree,matingdatavec,nsim=matensim,Q=qrates) 
  write.simmap(simmappy, file=paste(Sys.Date(),MateParam,SongParam,matensim,"simmaps",".txt",sep=""))
  simmapsdone <- Sys.time()
  simmaptime <- simmapsdone - starttimebrownie
  
  print(paste("Simmaps generated. That step took this much time: ", simmaptime, ".  Starting for loop with ", matensim, " loops."))
  
  browniedata <- data.frame(MatePar=character(matensim),SongPar=character(matensim),Pval=numeric(matensim),ERRate=numeric(matensim),ERloglik=numeric(matensim),ERace=numeric(matensim),ARDRate0=numeric(matensim),ARDRate1=numeric(matensim),ARDloglik=numeric(matensim),ARDace=numeric(matensim),k2=numeric(matensim),convergence=character(matensim),simmapnumber=integer(matensim),phylanovaP=numeric(matensim),stringsAsFactors = FALSE)
  brownied <- set.seed(10)
  browniedf <- set.seed(10)
  
  for (i in 1:matensim) {
    simmapfor <- simmappy[[i]]
    brownieliteresults <- set.seed(10)
    
    tryCatch(
      expr = {
        withTimeout(expr={
          
          brownieliteresults <- brownie.lite(simmapfor,songdatavec,maxit=75000)
          browniedata[i,3] <- brownieliteresults$P.chisq
          browniedata[i,4] <- brownieliteresults$sig2.single
          browniedata[i,5] <- brownieliteresults$logL1
          browniedata[i,6] <- brownieliteresults$a.single
          browniedata[i,7] <- brownieliteresults$sig2.multiple[1]
          browniedata[i,8] <- brownieliteresults$sig2.multiple[2]
          browniedata[i,9] <- brownieliteresults$logL.multiple
          browniedata[i,10] <- brownieliteresults$a.multiple
          browniedata[i,11] <- brownieliteresults$k2
          browniedata[i,12] <- as.character(brownieliteresults$convergence)
          browniedata[i,13] <- i}, timeout = 16, cpu=Inf, onTimeout = "error")
      },
      TimeoutException = function(ex) {browniedata[i,2:13]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"timeout",i);
      print(paste("timeout",i));
      })
    
    browniedata[i,1] <- MateParam
    browniedata[i,2] <- SongParam
    browniedata[i,14] <- phylanova
    
    # if (i %in% c(20,100,160,200,400,600,800,1000,1200,1400)) {
    #   write.csv(file = paste(Sys.Date(),"brownie",matensim,"sim",MateParam,SongParam,minmax,".csv",sep=""), browniedata) #cumulative brownie data results, saved during long process
    #   print("Saved data")
    # }
    if (i %in% seq(0,2000,by=50)) {
      print(paste("End brownie loop iteration",i,Sys.time()))
    }
  }  #end for loop 1:matensim
  print("End brownie loop")
  endtimebrownie <- Sys.time()
  looptime = endtimebrownie-starttimebrownie
  #write.csv(file = paste(Sys.Date(),"brownie",matensim,"sim",MateParam,SongParam,minmax,".csv",sep=""), browniedata)
  
  brownied <- browniedata
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
  
  
  #plot brownie distribution
  # if (MateParam == "Polygyny") {
  #   state0 <- "Monogamy"
  #   state1 <- "Polygyny"
  # } else if (MateParam == "EPP") {
  #   state0 <- "Low EPP"
  #   state1 <- "High EPP"
  # }
  # D0 <- density(browniedf$ARDRate0)
  # D1 <- density(browniedf$ARDRate1)
  # 
  # par(mar = c(1,1,2,1))
  # plot(D0,col="blue",
  #      xlim=c(min(c(D0$x,D1$x)),
  #             max(c(D0$x,D1$x))),
  #      ylim=c(min(c(D0$y,D1$y)),
  #             max(c(D0$y,D1$y))),
  #      main=paste(MateParam, SongParam, "all converged runs", " \nPval =", round(P.chisqAll,3), "/ Pval", state0, "higher =", round(P.chisq0over1,3), "/ Pval", state1, "higher =", round(P.chisq1over0,3)), cex.main = 0.75, xlab=paste("Rate of",SongParam,"evolution"),ylab="Frequency") 
  # lines(D1, col="red")
  # abline(v=browniedf$ERRate[1], lty = 2)
  # if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] > browniedf$ERRate[1]) {
  #   legend("topleft",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
  # } else if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] < browniedf$ERRate[1]) {
  #   legend("topright",legend = c(paste(state0),paste(state1), "Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2)) 
  # }
  
  return(browniedf[,1:14])
}