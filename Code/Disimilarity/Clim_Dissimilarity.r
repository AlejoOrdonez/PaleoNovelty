rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
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
# Biome
NPDB.Pollen <- NPDB.Pollen[c(NPDB.Pollen$ID1%in%NPDB.SiteUse),]
#####

#####
# Estimate the StdEucDis for Temperature over NBPD sites
TempStdEucDisList <- lapply(c("Sprg","Summ","Fall","Wint"),#c("Summ","Wint")
                            function(x){#(x<-"Summ")
                              Dist1 <- as.matrix(dist(NPDB.Clim[,paste0("BslTempMn.",x)])^2)
                              StdMtx <- matrix(NPDB.Clim[,paste0("BslTempSd.",x)],
                                               nrow = dim(Dist1)[1],
                                               ncol = dim(Dist1)[2],
                                               byrow = TRUE)
                              Dist2 <- Dist1/StdMtx # row wise 
                              return(Dist2)
                            })
TempStdEucDis <- Reduce("+",TempStdEucDisList)
# Estimate the StdEucDis for Precipitation over NBPD sites
PrecStdEucDisList <- lapply(c("Sprg","Summ","Fall","Wint"),#c("Summ","Wint")
                            function(x){
                              Dist1 <- as.matrix(dist(NPDB.Clim[,paste0("BslPrecMn.",x)])^2)
                              StdMtx <- matrix(NPDB.Clim[,paste0("BslPrecSd.",x)],
                                               nrow = dim(Dist1)[1],
                                               ncol = dim(Dist1)[2],
                                               byrow = TRUE)
                              Dist2 <- Dist1/StdMtx # row wise 
                              return(Dist2)
                            })
PrecStdEucDis <- Reduce("+",PrecStdEucDisList)
# Estimate the StdEucDis over over NBPD sites 
FullStdEucDis <- (TempStdEucDis+PrecStdEucDis)^0.5
#####

#####
# Define the Analogue Threshold using a ROC criteria
EnvTreshold.ROC <- roc(FullStdEucDis, # current time Taxon data.frame 
                       groups = NPDB.Pollen$BIOME # vector of group memberships
                       )
EnvTreshold.ROC
#####

#####
# Estimate the dissimilarity to present (TRaCE21 defined) of paleo-climates of cored sites.
# Load the Aggregated NEATOMA Data Base
NEATOMA.Agg <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")
# Load the climate of the NEATOMA Data Base
NEATOMA.Clim <- read.csv("./Data/PaleoClimate/Processed/NEATOMA_Clim.csv")


# Stdz Euclidean Distance for Temperature variables
PaleoTempStdEucDisList <- lapply(c("Sprg","Summ","Fall","Wint"),#c("Summ","Wint")
                                 function(x){#x <- "Sprg"
                                   Dist1 <- analogue::distance(data.frame(x=NPDB.Clim[,paste0("BslTempMn.",x)]),
                                                               data.frame(x=NEATOMA.Clim[,paste0("Temp.Mn.",x)]),
                                                               method = "SQeuclidean")
                                   StdMtx <- matrix(NPDB.Clim[,paste0("BslTempSd.",x)],
                                                    nrow = dim(Dist1)[1],
                                                    ncol = dim(Dist1)[2],
                                                    byrow = F)
                                   Dist2 <- Dist1/StdMtx # Col wise 
                                   return(Dist2)
                                 })
PaleoTempStdEucDis <- Reduce("+",PaleoTempStdEucDisList)
# Stdz Euclidean Distance for Precipitation variables
PaleoPrecStdEucDisList <- lapply(c("Sprg","Summ","Fall","Wint"),#c("Summ","Wint")
                                 function(x){
                                   Dist1 <- analogue::distance(data.frame(x=NPDB.Clim[,paste0("BslPrecMn.",x)]),
                                                               data.frame(x=NEATOMA.Clim[,paste0("Prec.Mn.",x)]),
                                                               method = "SQeuclidean")
                                   StdMtx <- matrix(NPDB.Clim[,paste0("BslPrecSd.",x)],
                                                    nrow = dim(Dist1)[1],
                                                    ncol = dim(Dist1)[2],
                                                    byrow = F)
                                   Dist2 <- Dist1/StdMtx # Col wise 
                                   return(Dist2)
                                 })
PaleoTempStdEucDis <- Reduce("+",PaleoPrecStdEucDisList)
# Stdz Euclidean Distance 
FullPaleoStdEucDis <- (PaleoTempStdEucDis+PaleoTempStdEucDis)^0.5
colnames(FullPaleoStdEucDis) <- (NEATOMA.Clim$sites)
rownames(FullPaleoStdEucDis) <- NPDB.Clim$SITECODE

# Estimate the Min Distance to present for each Neatoma site
NEATOMA.Clim$minSESClim <- as.numeric(apply(FullPaleoStdEucDis,2,min))


# Save the climate dissimilarity output
NEATOMA.ClimDis <- list(FullDisMatrix = FullPaleoStdEucDis,
                        NEATOMA.ClimDis = NEATOMA.Clim[,-c(1,8:15)],
                        NEATOMA.Clim = NEATOMA.Clim,
                        NPDB.Clim = NPDB.Clim,
                        ROC.Cutoff = EnvTreshold.ROC)
saveRDS(NEATOMA.ClimDis,"./Results/StdEuDis_ClimDis.rds")


# Dummy plot (taking a even sub sample of sites across periods)
ClimDisBoot <- lapply(1:1000,
            function(i){
              SamplTmp <- do.call("c",lapply(1:21,
                                             function(x){
                                               sample(which(NEATOMA.Clim$Time==x), 10)
                                             }))
              tapply(NEATOMA.Clim$minSESClim[SamplTmp],
                     NEATOMA.Clim$Time[SamplTmp],
                     median)
            })

ClimDisBoot2 <- do.call("rbind",ClimDisBoot)
ClimDisBootQuant <- apply(ClimDisBoot2,2,quantile,c(0.0275,0.5,0.975))
plot(y = rev(ClimDisBootQuant[2,]),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Env Change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimularity (stndz Euc Dist)",
     ylim=range(ClimDisBootQuant))
lines(y = rev(ClimDisBootQuant[1,]),
      x = -21:-1,)
lines(y = rev(ClimDisBootQuant[3,]),
      x = -21:-1,)
abline(h=EnvTreshold.ROC$roc$ Combined$ optima)



