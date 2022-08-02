# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Utility file to manually test the tools



# load the functions
sapply(FUN = source,
       X = list.files("src",
                      pattern = ".R$",
                      full.names = T))




# Progeny BLUP variance and expected value calculation `progenyBlupVarExp` ----

genoFile <- 'data/geno/breedGame_phasedGeno.vcf.gz'
crossTableFile <- 'data/crossingTable/breedGame_small_crossTable.csv'
SNPcoordFile <- 'data/SNPcoordinates/breedingGame_SNPcoord.csv'
markerEffectsFiles <- 'data/markerEffects/breedGame_markerEffects.csv'
outFile <- NULL


progenyBlupVarExp <- calc_progenyBlupVarExp(genoFile = genoFile,
                                            crossTableFile = crossTableFile,
                                            SNPcoordFile = SNPcoordFile,
                                            markerEffectsFiles = markerEffectsFiles,
                                            outFile = outFile)


progenyBlupVarExp$cross <- paste0(progenyBlupVarExp$ind1,
                                  '_X_',
                                  progenyBlupVarExp$ind2)

library(plotly)

fig <- plot_ly(
  data = progenyBlupVarExp,
  x = ~ cross,
  y = ~ exp,
  type = 'scatter',
  mode = 'markers',
  error_y = ~ list(array = sqrt(var),
                   color = '#000000'),
  hoverinfo = 'text',
  text = apply(progenyBlupVarExp, 1, function(l) {
    paste(names(l), ":", l, collapse = "\n")
  })
)
fig

