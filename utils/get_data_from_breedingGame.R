# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Utility file to get data from a breeding game session (of https://github.com/timflutre/PlantBreedGame)

# The best way to do it is to have access to the game data of the session and run the breeding game locally.
# The be sure the "breeder" you are interedted have the game master status.





# get pedigree data ----
# 1st download all the "Plant material Data"
# of the breeder in the `baseFolder`
# (I used UT-previous sessions: Todai_session3_2021-11_training, Sakurai2)
baseFolder <- "~/Downloads"
pedFiles <- list.files(baseFolder, "IndList", full.names = T)

# read plant material data to create the pedigree table
ped <- do.call(rbind, lapply(pedFiles, read.table, sep = "\t", header = T))
ped <- ped[, c("child", "parent1", "parent2")]
colnames(ped)[1] <- "ind"

# create list of the founder population
allFounders <- sprintf(
  paste0(
    "Coll",
    "%0", floor(log10(1000)) + 1, "i"
  ),
  1:1000
)
usedFounders <- allFounders[allFounders %in% do.call(c, lapply(ped, identity))]

# add the founders to the pedigree data
ped <- dplyr::bind_rows(
  ped,
  data.frame(
    ind = allFounders, # usedFounders
    parent1 = NA,
    parent2 = NA
  )
)

# save the well formated pedigree file with the founders
newPedFile <- paste0(baseFolder, "/breedingGame_pedigree.csv")
write.csv(ped, newPedFile,
  na = "", row.names = F
)


# get geno data ----
# we will request a new genotyping containing all the individuals
# create geno request:
allInds <- unique(do.call(c, lapply(ped, identity)))
allInds <- allInds[!is.na(allInds)]
req <- data.frame(
  ind = allInds,
  task = "geno",
  details = "hd"
)

# for a smaller request:
usedInds <- c(allInds[!allInds %in% allFounders], usedFounders)
req <- data.frame(
  ind = usedInds,
  task = "geno",
  details = "hd"
)

newReq <- paste0(baseFolder, "/breedGame_genoReq.txt")
write.table(req, newReq,
  na = "", row.names = F, sep = "\t", quote = FALSE
)

# get pheno data ----
# we will request a new phenotyping of the founder (on only one year)

# create pheno request:
req <- data.frame(
  ind = allFounders,
  task = "pheno-field",
  details = "4"
)

newReq <- paste0(baseFolder, "/breedGame_phenoReq.txt")
write.table(req, newReq,
  na = "", row.names = F, sep = "\t", quote = FALSE
)

# you can now request the phenotyping in the breeding game and download the result
# be sure the pathogen column is `FALSE` (if not, the `trait3` will have effect on `trait1`, and a BLUP calculation in needed)

# clean the created pheno file
pheno <- read.table(
  paste0(
    baseFolder,
    "/Result_pheno-field_breedGame_phenoReq_2015-01-01.txt.gz"
  ),
  header = T
)

# remove unessesary coluns
pheno <- pheno[, c("ind", "plot", "trait1", "trait2")]

# save the clean phenotype file
phenoFile <- paste0(baseFolder, "/breedGame_pheno.csv")
write.csv(pheno, phenoFile,
  na = "", row.names = F
)

# get markerEffects ----

baseFolder <- "~/Desktop/data"
load(paste0(baseFolder, "/truth/g0.RData"))
load(paste0(baseFolder, "/truth/p0.RData"))
p0$mu
trueMarkerEffects <- p0$Beta


markerEffects <- readMarkerEffects(file = "./data/markerEffects/breedGame_markerEffects.csv")
markerEffects <- readMarkerEffects(file = "./data/markerEffects/breedGame_markerEffects.json")



head(trueMarkerEffects[rownames(markerEffects$SNPeffects), ])
head(markerEffects$SNPeffects)

markerEffects2traits <- list(
  intercept = p0$mu,
  SNPeffects = trueMarkerEffects[rownames(markerEffects$SNPeffects), ]
)

jsonList <- lapply(names(markerEffects2traits$intercept), function(trait) {
  coefList <- as.list(markerEffects2traits$SNPeffects[, trait, drop = FALSE])
  names(coefList) <- row.names(markerEffects2traits$SNPeffects)
  list(
    intercept = p0$mu[trait],
    coefficients = coefList
  )
})
names(jsonList) <- names(markerEffects2traits$intercept)

outFile <- "data/markerEffects/breedGame_markerEffects_2traits.json"
(jsondta <- jsonlite::toJSON(jsonList,
  auto_unbox = TRUE,
  pretty = T,
  digits = NA,
  na = "string"
))

writeLines(jsondta, con = outFile)
