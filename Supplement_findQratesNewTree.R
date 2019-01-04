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

#findQrates - to be used within matingfunction to set the rates of transition between states for the building of simmaps for brownie
#Can also be used to generate simmaps with the computed rates with plot = TRUE

findQrates <- function(MateParam=c("Polygyny","EPP","OC","none"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity","OC","none"), matemodel = "ARD", plot=FALSE, matensim = 1500, newtree = newtree, withsuboscines = FALSE, newdata = FALSE) {
require(phytools)
source(file = "Supplement_subsetbirddataNewTree.R")
  
    output <- list()
  

subsetoutput <- subsetbirddata(MateParam, SongParam, withsuboscines = TRUE, newdata = newdata, newtree=newtree)  

  matingtree <- subsetoutput$ditree  
  matecol <- subsetoutput$matecol
  matingdatavec <- subsetoutput$matevec

  if (withsuboscines == TRUE) {  # for file naming
    subosc = "subosc"
  } else {subosc = NULL}
  
  ##########
  if (MateParam != "OC") { #this does the test for justification for ARD model for mating system
    ERmodel <- ace(matingdatavec,matingtree, type="discrete",model = "ER")
    ARDmodel <- ace(matingdatavec,matingtree, type="discrete",model = "ARD")
    anovaERARD <- anova(ERmodel,ARDmodel)
    testrates=c(ERmodel$rates,ERmodel$lik.anc[1,], ARDmodel$rates,ARDmodel$lik.anc[1,])
    output$anovaERARD <- anovaERARD
    output$ERrates <- paste("ERrates",ERmodel$rates)
    output$ERlikanc <- paste("ERlik.anc",ERmodel$lik.anc[1,])
    output$ARDrates <- paste("ARDrates",ARDmodel$rates)
    output$ARDlikanc <- paste("ARDlik.anc",ARDmodel$lik.anc[1,])
  
    
    matsize <- length(ARDmodel$rates)
    if(matsize == 1){matsize<- matsize+1}
    #create matrix of correct size; make diag negative; multiply by rates
    qrates <- matrix(rep(1,matsize^2),matsize,matsize)
    diag(qrates) <- -1
    qrates <- qrates*ARDmodel$rates  
    #pulls the state names and sets as col and row names
    rownames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    colnames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    
    if (plot == TRUE) {
    mateonlysimmap <- make.simmap(matingtree,matingdatavec,model = "ARD", nsim = 3) #makes six simmaps for viewing purposes
    mateonlysimmapQset <- make.simmap(matingtree,matingdatavec,model = "ARD", nsim = 3, Q = qrates)
    pdf(file = paste(Sys.Date(),MateParam,SongParam,subosc,"egSimmaps.pdf",sep=""),height=10,width=6)
    layout(matrix(1:6,nrow = 2,ncol=3))
    for (i in 1:3) {
      simmap <- mateonlysimmap[[i]]
      plotSimmap(simmap,fsize=0.2, lwd = 0.8)
      numrates1 <- lapply(simmap[[7]][,1],round,digits=9)
      title(paste("make.simmap"," \nQrates:", numrates1[1],numrates1[2]),cex.main = 0.5)
      simmap2 <- mateonlysimmapQset[[i]]
      plotSimmap(simmap2,fsize=0.2,lwd=0.8)
      numrates <- lapply(ARDmodel$rates,round,digits=9)
      title(main=paste("ARDmodel"," \nQrates:",numrates[1],numrates[2]),cex.main = 0.5)
    }
    dev.off()
    }

  } else {
    ERmodel <- ace(matingdatavec,matingtree, type="discrete",model = "ER")
    ARDmodel <- ace(matingdatavec,matingtree, type="discrete",model = "ARD")
    anovaERARD <- anova(ERmodel,ARDmodel)
    
    matsize <- length(ERmodel$rates)
    if(matsize == 1){matsize<- matsize+1}
    #create matrix of correct size; make diag negative; multiple by rates
    qrates <- matrix(rep(1,matsize^2),matsize,matsize)
    diag(qrates) <- -1
    qrates <- qrates*ERmodel$rates  
    #pulls the state names and sets as col and row names
    rownames(qrates) <- dimnames(ERmodel$lik.anc)[[2]]
    colnames(qrates) <- dimnames(ERmodel$lik.anc)[[2]]
  }
  
  
  output$qrates <- qrates
  return(output)

}