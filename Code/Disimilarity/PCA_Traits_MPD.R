rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
require(snowfall)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
# Load the Aggregated North American Pollen Data Base
NPDB.Pollen <- read.csv("./Data/Pollen/Processed/NPDB_Agg.csv")
# Load the climate for the North American Pollen Data Base
NPDB.Clim <- read.csv("./Data/PaleoClimate/Processed/NPDB_Clim.csv")

# Select only sites with Pollen data and Climate data
NPDB.SiteUse <- table(c(NPDB.Clim$ID1,NPDB.Pollen$ID1))
NPDB.SiteUse <- as.numeric(names(NPDB.SiteUse)[(NPDB.SiteUse==2)])
# Climate
NPDB.Clim <- NPDB.Clim[c(NPDB.Clim$ID1%in%NPDB.SiteUse),]
# Aggregated composition
NPDB.Pollen <- NPDB.Pollen[c(NPDB.Pollen$ID1%in%NPDB.SiteUse),]
# Load the Aggregated NEATOMA Data Base
NEATOMA.Agg <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")

# Which taxa to use - ensure match between NPDS and NEATOMA
TaxaUse <- table(c(names(NEATOMA.Agg[,-c(1:6)]),names(NPDB.Pollen[,-c(1:10)])))
TaxaUse <- names(TaxaUse)[TaxaUse==2]
#####

#####
# Load the trait data
Traits <- read.csv2("./Data/Traits/Raw/TraitsSummary.csv")
Traits <- Traits[c(Traits$Name%in%TaxaUse),]
# Seeds 
SWT <- apply(Traits[,c("Min..of.Seed.Size..mg.","Max..of.Seed.Size..mg.")],
             1,
             function(x){
               ifelse(is.na(x[1]),
                      NA,
                      runif(1,log10(x[1]),log10(x[2])))
             })
SWT <- scale(SWT)

# Height
Hmax <- apply(Traits[,c("Min..of.Height..m.","Max..of.Height..m.")],
              1,
              function(x){
                ifelse(is.na(x[1]),
                       NA,
                       runif(1,x[1],x[2]))
              })
Hmax <- scale(Hmax)
# LMA
LMA <- apply(Traits[,c("Min..of.LMA..log...g.m2.","Max..of.LMA..log...g.m2.")],
             1,
             function(x){
               ifelse(is.na(x[1]),
                      NA,
                      runif(1,x[1],x[2]))
             })
LMA <- scale(LMA)

# Build a trait dataset
TraitsDBS <- data.frame(SWT = SWT,
                        Hmax = Hmax,
                        LMA = LMA)
row.names(TraitsDBS) <- Traits$Name
TraitsDBS <- TraitsDBS[complete.cases(TraitsDBS),]

#turn traits into PCoA axes
GowdisDfrm <- FD::gowdis(TraitsDBS)
PCoATraits <- ape::pcoa(GowdisDfrm)
PCoATraitsEgVal <- PCoATraits $ vectors

# Which taxa to use - ensure match between NPDS and NEATOMA and TRAITS
TaxaUse <- row.names(PCoATraitsEgVal)
# Estimate the Trait distance in PCoA Space                       
TraitsDist <- as.matrix(dist(PCoATraitsEgVal))
#####


#######################################################################################################################
#####
# Load the North American Pollen Data Base and its translation to the Paleo data
# raw NPDB
SrtTime1 <- Sys.time()
BootTraits <- lapply(1:100,
                     function(RunIt){
                       #####
                       # Build a species names per site for the NPDS
                       SppPerSiteNPDB <- apply(NPDB.Pollen[,TaxaUse],
                                               1,
                                               function(x){TaxaUse[x>0]})
                       #####
                       #####
                       # Trait Mean pairwise distance contrast
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("TraitsDist")
                       SrtTime<-Sys.time()
                       MPDBaseline <- lapply(SppPerSiteNPDB,
                                             function(i){
                                               sapply(SppPerSiteNPDB,
                                                      function(x){
                                                        mean(TraitsDist[i,x])
                                                      })
                                             })
                       MPDBaseline <- do.call("rbind",MPDBaseline)
                       sfStop()
                       Sys.time()-SrtTime
                       #####
                       
                       #####
                       # Define the analogue threshold
                       # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                       TraitMPDDist.ROC <- roc(MPDBaseline, # current time Taxon data.frame 
                                               groups = NPDB.Pollen$BIOME # vector of group memberships
                       )
                       TraitMPDDist.ROC
                       #####
                       
                       #####
                       # Build a species per site Data for the NEATOMA dataset
                       PollenSppSiteNames <- apply(NEATOMA.Agg[,TaxaUse],
                                                   1,
                                                   function(x){TaxaUse[x>0]
                                                   })
                       #####
                       
                       #####
                       # Estimate distance to each paleo assemblage
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("PollenSppSiteNames")
                       sfExport("TraitsDist")
                       SrtTime<-Sys.time()
                       PaleoToPresentTraitDist <- sfLapply(PollenSppSiteNames,
                                                           function(i){
                                                             out <- sapply(SppPerSiteNPDB,
                                                                           function(x){
                                                                             mean(TraitsDist[i,x],na.rm=T)
                                                                           })
                                                             min(out,na.rm=T)
                                                           })
                       sfStop()
                       Sys.time()-SrtTime
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(DataFrm = data.frame(Time=-21:-1,
                                                             MPDDist = rev(tapply(do.call("c",PaleoToPresentTraitDist),
                                                                                  NEATOMA.Agg$Time,
                                                                                  median,
                                                                                  na.rm=T))),
                                        RocTresh = TraitMPDDist.ROC)
                       return(Out.list)
                     })
Sys.time()-SrtTime1
saveRDS(BootTraits,"./Results/BootPCoATraitsMPD.rds")
