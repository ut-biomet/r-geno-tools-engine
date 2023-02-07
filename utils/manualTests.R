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






# GWAS
genoF <- 'data/geno/testMarkerData01.vcf.gz'
phenoF <- 'data/pheno/testPhenoData01.csv'

res_gwas <- run_gwas(genoFile = genoF,
                 phenoFile = phenoF,
                 trait = 'Flowering.time.at.Arkansas',
                 thresh_maf = 0.05,
                 thresh_callrate = 0.9,
                 test = 'score',
                 response = 'quantitative')

run_resAdjustment(gwasFile = res_gwas$file)



manPlot(gwas = res_gwas,
        adj_method = 'bonferoni',
        thresh_p = 0.05,
        interactive = FALSE,
        title = )



# relmat
genoF <- 'data/geno/breedGame_geno.vcf.gz'
calc_genoRelMAt(genoFile = genoF)

g <- readGenoData(genoF)

outFile <- tempfile()
draw_relHeatmap(relMatFile = "/tmp/RtmpWSw4RM/file1e5d76636ac49.csv",
                outFile = outFile)

draw_relHeatmap(relMatFile = "/tmp/RtmpNDRRPq/file2873b580e67cb.csv",
                interactive = FALSE,
                outFile = outFile)




# ped network
ped <- "data/pedigree/breedGame_pedigree.csv"

p <- readPedData(ped)

outFile <- tempfile()
draw_pedNetwork(pedFile = ped, removeUnusedInds = TRUE,
                outFile = outFile)








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
p <- draw_progBlupsPlot(outFile, sorting = 'alpha', outFile = graphFile)


# crossing simulation ----
dir.create('tmp')
outFile <- 'tmp/crossSim.vcf.gz'
crosSim <- crossingSimulation(
  genoFile = genoFile,
  crossTableFile = crossTableFile,
  SNPcoordFile = SNPcoordFile,
  nCross = 100,
  outFile = 'tmp/crossSim.vcf.gz')


geno <- gaston::read.vcf('tmp/crossSim.vcf.gz')
geno <- gaston::as.matrix(geno)

markerEffects <- readMarkerEffects(markerEffectsFile)

geno <- geno[, row.names(markerEffects$SNPeffects)]

bvDta <- data.frame(bv = as.vector(geno %*% as.matrix(markerEffects$SNPeffects)),
                    ind = as.factor(sub('-\\d+$', '', row.names(geno))))

#plot
p <- plot_ly(type = "box",
        # mode = "markers",
        data = bvDta,
        x = ~ind,
        y = ~bv,
        boxpoints = "all",
        jitter = 1,
        pointpos = 0
) %>% add_markers(
    data = progenyBlups,
    x = paste0(progenyBlups$ind1,'_X_', progenyBlups$ind2),
    y = ~ blup_exp,
    # type = 'scatter',
    # mode = 'markers',
    error_y = ~ list(array = 2*sqrt(blup_var),
                     type = 'data',
                     color = '#000000'),
    hoverinfo = 'text',
    text = apply(progenyBlups, 1, function(l) {
      paste(names(l), ":", l, collapse = "\n")
    })
  )
outFile <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(plotly::partial_bundle(p),
                        outFile, selfcontained = TRUE)



# SNP coord with missing ids ----
SNPcoordFile <- 'tests/data/breedingGame_SNPcoord_missingIds.csv'
SNPcoordFile <- 'tests/data/breedingGame_SNPcoord_missingPhysPos.csv'
readSNPcoord(SNPcoordFile)

