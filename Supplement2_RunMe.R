#Run Me
#Revision edition
#BayesTraits
.BayesTraitsPath <- "~/Documents/VANDERBILT/CreanzaLab/BayesTraits/Sandbox/BayesTraitsV2"

setwd("~/Documents/VANDERBILT/CreanzaLab/BayesTraits/Sandbox")

#GeneTree
#Perform and plot BayesTraits on GeneTree
source("Supplement2_genetreebtwcycle.R")
  #Sources: "Supplement2_btwfunction_GeneTree.R"
    #source(file = "Supplement_subsetbirddata_newtree.R")
    #source(file = "Supplement2_BayesPlotsGeneTree.R") #contains plotBTjacks and transitionplots
MateParams <- c("Polygyny")
SongParams <- c("Syllsong")
cyclematesongsGene(jackknife = FALSE, csvsout = TRUE, nsim = 4)
#Returns one csv file & one plot per Mate/Song combination


#MultiPhy
#need to be in folder with "sampleMatezillaHack#.nex" # = 10-29
#cyclematesongs will create .csv file outputs from BayesTraits
#plotBTmultiphy will plot from these .csv files
source("Supplement2_btwfunction_newtree.R") #contains btwfunction and cyclematesongs
MateParams = c("Polygyny")
SongParams = c("Syllsong")
cyclematesongs(jackknife = FALSE, csvsout = TRUE, plotbt = FALSE,nsim = 2, treevec = 10:11)

source("Supplement2_BayesPlotsmultiphyAveraged.R")
pdf(file = "ALLBayesplotsMultiphyAveraged.pdf",height = 16, width = 12)
par(mfrow=c(nrow=6,ncol =4))
par(mar=c(2,1,3,2))
plotBTmultiphy(matelist = MateParams, songlist = SongParams, nsims = 2)  #plots to file "ALLBayesplotsMultiphyAveraged.pdf"
dev.off()

