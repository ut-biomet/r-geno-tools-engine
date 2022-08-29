# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Utility file to manually test the tools


# packages

# load the functions
sapply(FUN = source,
       X = list.files("src",
                      pattern = ".R$",
                      full.names = T))




# Progeny BLUP variance and expected value calculation `progenyBlupVarExp` ----

genoFile <- 'data/geno/breedGame_phasedGeno.vcf.gz'
crossTableFile <- 'data/crossingTable/breedGame_small_crossTable.csv'
SNPcoordFile <- 'data/SNPcoordinates/breedingGame_SNPcoord.csv'
markerEffectsFile <- 'data/markerEffects/breedGame_markerEffects.csv'
outFile <- 'data/results/progenyBlupEstimation.json'


progenyBlups <- calc_progenyBlupEstimation(
  genoFile = genoFile,
  crossTableFile = crossTableFile,
  SNPcoordFile = SNPcoordFile,
  markerEffectsFile = markerEffectsFile,
  outFile = outFile
)


plotBlup(progenyBlups)

graphFile <- tempfile(fileext = ".html")
draw_progBlupsPlot(outFile, sorting = 'alpha', outFile = graphFile)
