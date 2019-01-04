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

GLMMfunc <- function() {
  source('GLMM_subsetbirddataEfficientFlex.R')
  glmmdata <- "DataGLMMformatted_ratecontR.csv"
  
songlist <- c("Syllrep","Syllsong","Song","Duration","Interval","Rate","Continuity")
resultslist <- list()
for (i in 1:7) {
  songparam <- songlist[i]
  print(songparam)
  subsetglmm <- subsetbirddata(c("Polygyny","EPP",songparam), newdata = glmmdata)
  glmmdf <- subsetglmm$df
  subtree <- subsetglmm$bothtree
  matecolvec <- subsetglmm$matecolvec
  songcol <- matecolvec[length(matecolvec)]
  colnames(glmmdf)[1] <- "species"
  Tree<- subtree
  Tree <- multi2di(Tree,0.0001)
  Tree <- force.ultrametric(Tree)
  Tree$edge.length[which(Tree$edge.length == 0)] <- .0001
  Tree$tip.label <- gsub(" ", "_", Tree$tip.label)
  Tree <- force.ultrametric(Tree)
  TreeAinv<-inverseA(Tree, nodes="TIPS")$Ainv
  animal <- glmmdf$species
  glmmdf <- cbind(animal,glmmdf)
  glmmdf[,songcol] <- log(glmmdf[,songcol],10)
  Fix <-  paste(songcol ," ~ Final.polygyny*Final.EPP")
  Ran <- "~animal+species"

  if (songparam %in% c("Rate","Continuity")) {
    Ran <- "~animal"
  }

  mcmcglmmout <- MCMCglmm(as.formula(Fix), random = as.formula(Ran), data = glmmdf, ginverse = list(animal = TreeAinv), pr = TRUE)
  glmmsumm <- summary(mcmcglmmout)
  solutions <- glmmsumm$solutions
  resultslist[[i]] <- list(Fix,Ran,solutions)
  
  
  ### posterior predictive checks
  ### 

 # ConvergePlot(mcmcglmmout, songparam)
  
  predictions1 <- predict(mcmcglmmout, marginal=NULL, posterior="all", type="response", interval = "confidence")
  ##plot histogram of mean of predictions and line of mean of actual data
  pdf(file = paste0(Sys.Date(),"GLMMpredictmeansLine",songparam,".pdf"))
  df <- as.data.frame(cbind(glmmdf[,songcol], predictions1))
  print(ggplot(df, aes(x = df$V1, y = df$fit)) +
          geom_point(size = 4) +
          geom_errorbar(aes(ymax = df$upr, ymin = df$lwr)) + labs(x = paste("log10",songcol), y = "predictions"))
  dev.off()
  ### trying to do posterior predictive checks
  set.seed(NULL)
    pred_means <- vector()
    predictions2 <- simulate.MCMCglmm(mcmcglmmout, nsim = 1000)
    #print(predictions)
    pred_means <- colMeans(predictions2)
    print(pred_means)
    #plot histogram of mean of predictions and line of mean of actual data
    pdf(file = paste0(Sys.Date(),"GLMMpredictmeansHist",songparam,".pdf"))
    hist(pred_means, main=mcmcglmmout$Fixed[[1]][[2]], xlab="Predicted means")
    abline(v=mean(glmmdf[,songcol]))
    dev.off()
} #end for loop going through song features

sink(file = "GLMMout.txt")
resultslist
sink(file = NULL)
} # end function


# Check convergence
# ConvergePlot <- function(model, Name="", Alpha.Cex=2, PDF=TRUE){
#   if(Name != ""){
#     Dash <- "-log"
#   }else{Dash <- ""}
#   
#   if(PDF){
#     pdf(paste0(Sys.Date(),"GLMMConvergence",Dash,Name,".pdf"), height = 9, width = 8)
#     
#   }
#   ncols <- ceiling((ncol(model$Sol)+ncol(model$VCV))/2)
#   par(mfrow=c(6,4), mgp=c(1.5,.5,0), mar=c(1,1,1,1))
#   Count <- 1
#   for(i in 1:ncol(model$Sol)){
#     CheckPlot(model$Sol[,i], colnames(model$Sol)[i], Count, Alpha.cex=Alpha.Cex)
#     Count <- Count+1
#   }
#   for(i in 1:ncol(model$VCV)){
#     CheckPlot(model$VCV[,i], colnames(model$VCV)[i], Count, Alpha.cex=Alpha.Cex)
#     Count <- Count+1
#   }
#   if(PDF){
#     dev.off()
#   }
# }
# CheckPlot <- function(modelVec, title="", pos=1, Alpha.cex=2.5){
#   par(mar = c(1.5,1.5,1.5,1.5))
#   title <- gsub("[.]units", "", title)
#   title <- gsub("[[:punct:]]", "", title)
#   plot.default(modelVec, type='l', ylab="", xlab= "Iterations",main=paste("Trace of ", title), cex.main=0.7)
# 
#   dense <- density(modelVec)
#   plot(dense, cex.main=0.7,
#        main=paste0("Density of ", title), font.main=2)
#   segments(modelVec, 0, modelVec, max(dense$y)*.03)
# }



