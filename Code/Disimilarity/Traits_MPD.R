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
TaxaUseObs <- table(c(names(NEATOMA.Agg[,-c(1:6)]),names(NPDB.Pollen[,-c(1:10)])))
TaxaUseObs <- names(TaxaUseObs)[TaxaUseObs==2]
#####

#######################################################################################################################

SrtTime1 <- Sys.time()
BootTraits <- lapply(1:100,
                     function(RunIt){
                       # Generate a possible trait realization based on the know distribution of values for a evaluated taxa
                       # Load the trait data
                       Traits <- read.csv2("./Data/Traits/Raw/TraitsSummary.csv")
                       Traits <- Traits[c(Traits$Name%in%TaxaUseObs),]
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
                       # Which taxa to use - ensure match between NPDS and NEATOMA and TRAITS
                       TaxaUse <- table(c(TaxaUseObs,row.names(TraitsDBS)))
                       TaxaUse <- names(TaxaUse)[TaxaUse==2]
                       # Estimate the Trait distance in RAW Space                       
                       TraitsDist <- as.matrix(dist(TraitsDBS[TaxaUse,]))
                       #####
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
                       #SrtTime<-Sys.time()
                       MPDBaseline <- sfLapply(SppPerSiteNPDB,
                                               function(i){
                                                 sapply(SppPerSiteNPDB,
                                                        function(x){
                                                          mean(TraitsDist[i,x])
                                                        })
                                               })
                       MPDBaseline <- do.call("rbind",MPDBaseline)
                       sfStop()
                       #Sys.time()-SrtTime
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
                       # Estimate distance of paleo assemblage to the NPDB
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("PollenSppSiteNames")
                       sfExport("TraitsDist")
                       #SrtTime<-Sys.time()
                       PaleoToPresentTraitDist <- sfLapply(PollenSppSiteNames,
                                                           function(i){
                                                             out <- sapply(SppPerSiteNPDB,
                                                                           function(x){
                                                                             mean(TraitsDist[i,x],na.rm=T)
                                                                           })
                                                             min(out,na.rm=T) # Only keep the Min Distance
                                                           })
                       NEATOMA.Agg$minSCDToPres <- do.call("c",PaleoToPresentTraitDist) # Add the Min distance to the Agg NEATOMA DATA
                       sfStop()
                       #Sys.time()-SrtTime
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(NEATOMA.CommDis = NEATOMA.Agg[,-c(1,7:45)],
                                        ROC.Cutoff = TraitMPDDist.ROC)
                       return(Out.list)
                     })

saveRDS(BootTraits,"./Results/BootTraitsMPD.rds")



# Load and plot the MPD for traits
rm(list=ls());gc()
require(snowfall)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")

BootTraits <- readRDS("./Results/BootTraitsMPD.rds")


# Estimate the Variability in estimates per Trait space iteration 
sfInit( parallel=TRUE, cpus=10)
sfExport("BootTraits")

TraitsMPDVar <- sfLapply(BootTraits,
                       function(x){#x<-BootTraits[[3]]
                         TraitsMPD <- x$NEATOMA.CommDis
                         MeanMPD <- lapply(1:1000,
                                           function(j){
                                           tapply(TraitsMPD$minSCDToPres[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                  TraitsMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                  median)})
                         MeanMPD <- do.call("cbind",
                                            MeanMPD)
                         Out <- data.frame(Time = 1:21,
                                           t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975))))
                         return(Out)
                       })
sfStop()

TraitsMPD <- data.frame(Time = 1:21,
                        X50 = apply(sapply(TraitsMPDVar,function(x){x$X50.}),1,median),
                        X2.5 = apply(sapply(TraitsMPDVar,function(x){x$X2.5.}),1,median),
                        X97.5 = apply(sapply(TraitsMPDVar,function(x){x$X97.5.}),1,median))

plot(y = rev(TraitsMPD$X50),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Trait Change", 
     xlab ="Time (kyrBP)",
     ylab = "MPD (Euclidean Dist)",
     ylim = range(TraitsMPD[-1]))
lines(y = rev(TraitsMPD$X2.5),
      x = -21:-1)
lines(y = rev(TraitsMPD$X97.5),
      x = -21:-1)

abline(h=median(sapply(BootTraits,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                       })))
