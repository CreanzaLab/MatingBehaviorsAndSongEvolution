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


btwfunction <- function(MateParam=c("Polygyny","EPP","OC"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Continuity","Rate","OC"), plot=TRUE, jackknife = FALSE, csvsout = FALSE, nsim = 10, newtreefile = FALSE) {
  output <- list()
  output$start <- Sys.time()
  print(Sys.time())
  require(btw)
  require(phytools)
  require(base)
  #require(mnormt)
  require(geiger)
  require(R.utils)
  
  source(file = "Supplement2_subsetbirddata_newtree.R")
 # source(file = "Supplement_BayesPlots.R")
  subset <- subsetbirddata(MateParam,SongParam, newtree = newtreefile)
  bothtree <- subset$ditree
  matecol <- subset$matecol
  songdf <- subset$df
  songcol <- subset$songcol
  
  output$lengthsongdf <- length(songdf$BirdtreeFormat)
  
  
  ##Jackknife test removing each family
  familyvecNoNone <- unique(songdf$Family)
  familyvec <- c("None",as.character(familyvecNoNone))

  if (plot == TRUE) {
    pdf(paste(Sys.Date(),"Bayes",MateParam,SongParam,"Pval",nsim,"repsJacks.pdf", sep=""), width = 7, height = 10)
    par(mfrow=c(round(length(familyvec)/4)+1,4))
#    layout(mat = matrix(nrow=20, ncol = 2, 1:40, byrow = TRUE))
    par(oma = c(1,3,3,1))
  }
    
    if (jackknife == TRUE) {

    transitions10reps <- data.frame()
     for (k in 1:length(familyvec)) {
   # for (k in 1:4) {
      jacksongdf <- songdf  #have to do this so songdf doesn't get whittled down every time the for loop loops
      jackedsongdf <- jacksongdf[which(jacksongdf$Family != familyvec[k]),]
      removedfamilydf <- jacksongdf[which(jacksongdf$Family == familyvec[k]),]
      
      notjackedvec <- bothtree$tip.label %in% as.character(jackedsongdf$BirdtreeFormat)
      dropforjack <- which(notjackedvec == FALSE)
      
      jacktree <- drop.tip(bothtree,tip = dropforjack)
      removedfamily <- familyvec[k]
      matevec <- as.character(jackedsongdf[,matecol])
      names(matevec) <- jackedsongdf$BirdtreeFormat
      songcontvec <- jackedsongdf[,songcol]
      songcontvec <- unique(songcontvec)
      print(paste("Without",familyvec[k],length(songcontvec)))
      LRstatall <- set.seed(10)
      LRpvalall <- set.seed(10)
      songcontvecall <- set.seed(10)
      transitions10reps <- data.frame(matrix(nrow=nsim*length(songcontvec), ncol=0))
      for (j in 1:nsim) {
        if (j == 1) {
          print(paste("Rep:",j, "Number of corrD points: 0", Sys.time()))
        } else if (j %in% c(10,20,40,60,80,100)) {
          print(paste("Rep:",j,"Number of corrD points",length(transitions10reps[,3]), Sys.time()))
        }
        LRstat <- set.seed(10)
        LRpval <- set.seed(10)
        transitions <- data.frame(matrix(nrow=length(songcontvec), ncol=0)) 
        for (i in 1:length(songcontvec)) {
          transandp <- set.seed(10)
          songdiscvec <- jackedsongdf[,songcol]
          thresh = songcontvec[i]
          songdiscvec[songdiscvec <= thresh] <- 0
          songdiscvec[songdiscvec > thresh] <- 1
          names(songdiscvec) <- jackedsongdf$BirdtreeFormat
          btwdfjack <- as.data.frame(cbind(as.character(jackedsongdf$BirdtreeFormat),matevec,as.character(songdiscvec)))
          nocorrD <- Discrete(jacktree, btwdfjack)
          corrD <- Discrete(jacktree, btwdfjack, dependent=TRUE)
          lrtestresults <- lrtest(corrD, nocorrD)
          LRstat[i] <- lrtestresults$LRstat
          LRpval[i] <- lrtestresults$pval
          transandp <- cbind(corrD,LRstat[i],LRpval[i],songcontvec[i])
          transitions <- rbind(transitions, transandp)
        } #end for i in 1:length(songcontvec)
        LRstatall <- c(LRstatall,LRstat)
        LRpvalall <- c(LRpvalall,LRpval)
        songcontvecall <- c(songcontvecall,songcontvec)
        transitions10reps <- rbind(transitions10reps,transitions)
      } #end for j in 1:nsim
      
      print(length(transitions10reps[,2]))
      
      if (csvsout == TRUE) {
        write.csv(transitions10reps,file=paste(Sys.Date(),"Bayes",MateParam,SongParam,nsim,"repsJack",familyvec[k],".csv",sep=""))
        csvfile <- paste(Sys.Date(),"Bayes",MateParam,SongParam,nsim,"repsJack",familyvec[k],".csv",sep="")
      }
      jackdf <- transitions10reps
      colnames(jackdf[,c(15:17)]) <- c("LRstat","LRpval","songcontvec")
      
      
      d = data.frame(songcontvecall,LRstatall, LRpvalall)
      
      if (plot == TRUE) {
        # par(mar = c(3,2,2,1))
        #   with(d, plot(songcontvecall, LRpvalall, pch=16, col="red3", 
        #                xlab=paste(SongParam, "threshold for binary (low/high) categorization"),ylab="p-value", main = paste(MateParam,"&",SongParam, "sign. of dependent correlation,","\nN =",length(songcontvec),"- Removed", paste(familyvec[k])), cex.main = 0.5, log = "x"
        #   ))
        plotBTjacks(MateParam,SongParam, d = jackdf, familysplit = removedfamily, nsim = nsim)
        
      } #end if plot ==TRUE
    } #end for loop covering family names
    #dev.off()
  } else { # if jackknife is anything but TRUE
    
    matevec <- as.character(songdf[,matecol])
    names(matevec) <- songdf$BirdtreeFormat
    songcontvec <- songdf[,songcol]
    songcontvec <- unique(songcontvec)
    #output$songcontvec <- songcontvec
    print(length(songcontvec))
    LRstatall <- set.seed(100)
    LRpvalall <- set.seed(100)
    songcontvecall <- set.seed(100)
    transitions100reps <- data.frame()
    for (j in 1:nsim) {
      if (j == 1) {
        print(paste("Rep:",j, "Number of corrD points: 0", Sys.time()))
      } else if (j %in% c(10,20,30,40,50,60,70,80,90)) {
        print(paste("Rep:",j,"Number of corrD points",length(transitions100reps[,3]),Sys.time()))
      }
      LRstat <- set.seed(10)
      LRpval <- set.seed(10)
      transitions <- data.frame(matrix(nrow=length(songcontvec), ncol=0))
      for (i in 1:length(songcontvec)) {
        transandp <- set.seed(10)
        songdiscvec <- songdf[,songcol]
        thresh = songcontvec[i]
        songdiscvec[songdiscvec <= thresh] <- 0
        songdiscvec[songdiscvec > thresh] <- 1
        names(songdiscvec) <- songdf$BirdtreeFormat
        btwdf <- as.data.frame(cbind(as.character(songdf$BirdtreeFormat),matevec,as.character(songdiscvec)))
        nocorrD <- Discrete(bothtree, btwdf)
        corrD <- Discrete(bothtree, btwdf, dependent=TRUE)
        lrtestresults <- lrtest(corrD, nocorrD)
        LRstat[i] <- lrtestresults$LRstat
        LRpval[i] <- lrtestresults$pval
        transandp <- cbind(corrD,LRstat[i],LRpval[i],songcontvec[i])
        transitions <- rbind(transitions, transandp)
      }
      LRstatall <- c(LRstatall,LRstat)
      LRpvalall <- c(LRpvalall,LRpval)
      songcontvecall <- c(songcontvecall,songcontvec)
      transitions100reps <- rbind(transitions100reps,transitions)
    }
    
    print(length(transitions100reps[,2]))
    if (csvsout == TRUE) {
      write.csv(transitions100reps,file=paste("Bayes",MateParam,SongParam,nsim,"reps.csv",sep=""))
      csvfile <- paste(Sys.Date(),"Bayes",MateParam,SongParam,nsim,"reps.csv",sep="")
      output$nojackdf<- transitions100reps
    }
    df <-  transitions100reps
    d = data.frame(songcontvecall,LRstatall, LRpvalall)
    colnames(df[,c(15:17)]) <- c("LRstat","LRpval","songcontvec")
    
    if (plot == TRUE) {
      transitionplots(MateParam,SongParam,df = df,newpdf = TRUE, nsim = nsim)
    } #end if plot == True 
    
  } #end for loop (jackknife == FALSE)
  dev.off()
  output$finish <- Sys.time()
  output$runtime <- output$finish - output$start
  return(output)
} #end btwfunction


cyclematesongs <- function(jackknife = FALSE, csvsout = FALSE, plotbt = TRUE, nsim = 10, treevec = trees) {
  for (n in 1:length(MateParams)) {
    for (m in 1:length(SongParams)) {
      MateParam = MateParams[n]
      SongParam = SongParams[m]
      for (i in treevec) {
        temptree <- paste0("samplematezillaHack",i,".nex")
      btwfunction(MateParam,SongParam,plot = plotbt,jackknife = jackknife,csvsout = csvsout, newtreefile = temptree, nsim = nsim)
      } #end for i in 10:29 going through trees
      #    transitionplots(MateParams[n],SongParams[n],csvfiles[n]) 
    } #end for m in 1:songparams
  } #end for n in 1:MateParams
} #end function BTfunc

