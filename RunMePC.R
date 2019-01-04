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

##Run Me
##Revision edition
##PC stuff
#


###EPP/Polygyny Boxplots and phylANOVA
source("EPPpolyBoxplotfunc.R") 
##sources: ("Supplement_subsetbirddata.R")
EPPpolyfunc()
#returns plots and csv "EPPpolyPhylANOVApairwise.csv"



###Multitree - generate csvs of brownie+phylANOVA analyses
##Outputs one csv file per mate/song combination and one pdf of the distribution of PhylANOVA pvalues
source("Supplement2_MateSongFuncNewTreeV2onepdf.R")
#source (file = "Supplement_subsetbirddataNewTree.R")
#source(file = "Supplement_findQratesNewTree.R")
#source(file = "Supplement_browniefunctionNewTree.R")
#mateparams <- c("Polygyny")
#songparams <- c("Syllrep")
 mateparams <- c(rep("Polygyny",6),rep("EPP",5))
 songparams <- c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Song","Syllsong","Syllrep","Interval","Duration")
 mateparams <- c("Polygyny","EPP")
 songparams <- c("Continuity","Continuity")
nsim = 20
treeN = 100
for (i in 1:length(mateparams)) {
  for (j in 1: length(songparams)) {
    matesongfunc(mateparams[i],songparams[j],brownie = TRUE,matensim=nsim, ntrees=treeN, minmax = NULL)
  }
}
##MultiTree - plot from already-made csv files
source("Supplement2_MultiTreephylANOVAbrownieonepdf.R")
pdf(file = "ALLbrowniephylanovaMultiTree.pdf", width = 9,height = 9)
par(mfrow=c(ncol = 4, nrow = 3))
for (j in 1:length(mateparams)) {
  mateparam = mateparams[j]
  songparam = songparams[j]
    matesongplots(mateparam,songparam, matensim = nsim, ntrees = treeN)
}
dev.off()

 


###GeneTree
source("Supplement2_MateSongFuncNewTree.R")
mateparams <- c("Polygyny")
songparams <- c("Syllrep")
#mateparams <- c("Polygyny","EPP")
#songparams <- c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity")
matesim = 20
for (i in 1:length(mateparams)) {
  for (j in 1: length(songparams)) {
    matesongfunc(MateParam = mateparams[i],SongParam = songparams[j], brownie = TRUE, nsim = matesim, newmatingtreefile = "GENEonlyConsensusTreeHack.nex")
  }
}#returns multiple csv's

##Plots from already-made csv files
source("Supplement2_MateSongFuncNewGENETreeOnePDFplot.R")
pdf(file = "ALLbrowniephylanovaGenetree.pdf", width = 11,height = 9)
par(mfrow=c(ncol = 4, nrow = 4))
for (i in 1:length(songparams)) {
  songparam <- songparams[i]
  for (j in 1:length(mateparams)) {
    mateparam <- mateparams[j]
    matesongGeneplot(mateparam,songparam,matensim = matesim)
  }
}
dev.off()
#returns single pdf


###PGLS
source("Supplement2_pglsresults.R")
  ##Sources GLMM_subsetbirddataEfficient.R
pglsfunc()
##Outputs csv file of results "pglsresults.csv


###GLMM
source("Supplement2_EPPpolyGLMM.R")
#sources:   source('GLMM_subsetbirddataEfficientFlex.R')
GLMMfunc()
#returns single text document GLMMout.txt with summaries of MCMCglmm tests for all song parameters; returns two pdf plots per song parameter
