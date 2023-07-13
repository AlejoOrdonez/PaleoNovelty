rm(list=ls());gc()
require(ape)
require(snowfall)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")

#####
# Super tree data
#Load a super tree
MegaTree <- read.tree("./Data/Phylogeny/ALLMB.tre")
# Spp Translation
SppTrans <- read.csv2("./Data/Phylogeny/SppUse_Pylo_PLANTS.csv")
# Remove taxa not in the Magatree
SppTrans <- SppTrans[c(SppTrans$name.for.phylo%in%MegaTree$tip.label),]
## Prune the Phylogeny
SppPhylo <- keep.tip(MegaTree,MegaTree$tip.label[c(MegaTree$tip.label%in%SppTrans$name.for.phylo)])
# Estimate the Pairwise distance between species
SppDist <- cophenetic(SppPhylo)
#####

#####
# Community data NPDB and NEATOMA 
# Load the Aggregated North American Pollen Data Base
NPDB.Pollen <- read.csv("./Data/Pollen/Processed/NPDB_Agg.csv")
# Load the climate for the North American Pollen Data Base
NPDB.Clim <- read.csv("./Data/PaleoClimate/Processed/NPDB_Clim.csv")
# Select only sites with Pollen data and Climate data
NPDB.SiteUse <- table(c(NPDB.Clim$ID1,NPDB.Pollen$ID1))
NPDB.SiteUse <- as.numeric(names(NPDB.SiteUse)[(NPDB.SiteUse==2)])
# Aggregated composition
NPDB.Pollen <- NPDB.Pollen[c(NPDB.Pollen$ID1%in%NPDB.SiteUse),]
# Load the Aggregated NEATOMA Data Base
NEATOMA.Agg <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")
#####

#####
# List of Spp to use in all three data sets
TaxaUse <- table(c(names(NEATOMA.Agg)[-c(1:7)],
                   names(NPDB.Pollen)[-c(1:10)],
                   unique(SppTrans$Translation)))
TaxaUse <- names(TaxaUse)[c(TaxaUse==3)]
#####


BootPhylo <- lapply(1:100,
                     function(RunIt){
                       SrtTime1 <- Sys.time()
                       #####
                       # Random selection of pylogeny tips
                       SppTrasTbl<- SppTrans[SppTrans$Translation%in%TaxaUse,] # Ensure the matching spp are used
                       SppTrasTbl <- table(SppTrasTbl$name.for.phylo,SppTrasTbl$Translation)
                       SppTrasTbl <- apply(SppTrasTbl,2,function(x){sample(rownames(SppTrasTbl)[x==1],1)})
                       # Generate a sample soecific phylogenetic distance matrix
                       SppDistTmp <- SppDist[as.character(SppTrasTbl),as.character(SppTrasTbl)]
                       dimnames(SppDistTmp)[[1]] <- dimnames(SppDistTmp)[[2]] <- names(SppTrasTbl)
                       #####
                       #####
                       # Build a species names per site for the NPDS
                       SppPerSiteNPDB <- apply(NPDB.Pollen[,TaxaUse],
                                               1,
                                               function(x){TaxaUse[x>0]})
                       SppPerSiteNPDB <- SppPerSiteNPDB[sapply(SppPerSiteNPDB,length)>0]
                       #####
                       
                       #####
                       # Phylogenetic Mean pairwise distance contrast
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("SppDistTmp")
                       #SrtTime<-Sys.time()
                       MPDBaselineList <- sfLapply(SppPerSiteNPDB,
                                               function(i){#(i<-SppPerSiteNPDB[[1]])
                                                 sapply(SppPerSiteNPDB,
                                                        function(x){#(x<-SppPerSiteNPDB[[2]])
                                                          mean(SppDistTmp[i,x])
                                                        })
                                               })
                       MPDBaseline <- do.call("rbind",MPDBaselineList)
                       
                       sfStop()
                       #Sys.time()-SrtTime
                       #####
                       #####
                       # Define the analogue threshold
                       # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                       BiomeUse <- NPDB.Pollen$BIOME[sapply(apply(NPDB.Pollen[,TaxaUse],1,function(x){TaxaUse[x>0]}),length)>0]
                       PhyloMPDDist.ROC <- roc(MPDBaseline[c(BiomeUse!=names(table(BiomeUse))[table(BiomeUse)==1]),
                                                           c(BiomeUse!=names(table(BiomeUse))[table(BiomeUse)==1])], # current time Taxon data.frame 
                                               groups =  BiomeUse[c(BiomeUse!=names(table(BiomeUse))[table(BiomeUse)==1])]# vector of group memberships
                       )
                       PhyloMPDDist.ROC
                       #####
                       
                       #####
                       # Build a species per site Data for the NEATOMA dataset
                       NEATOMA.Agg2 <- NEATOMA.Agg
                       PollenSppSiteNames <- apply(NEATOMA.Agg2[,TaxaUse],
                                                   1,
                                                   function(x){TaxaUse[x>0]
                                                   })
                       #####
                       
                       #####
                       # Estimate distance of paleo assemblage to the NPDB
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("PollenSppSiteNames")
                       sfExport("SppDistTmp")
                       #SrtTime<-Sys.time()
                       PaleoToPresentTraitDist <- sfLapply(PollenSppSiteNames,
                                                           function(i){#i<-PollenSppSiteNames[[16]]
                                                             out <- sapply(SppPerSiteNPDB,
                                                                           function(x){
                                                                             mean(SppDistTmp[i,x],na.rm=T)
                                                                           })
                                                             min(out,na.rm=T) # Only keep the Min Distance
                                                           })
                       NEATOMA.Agg2$minSCDToPres <- do.call("c",PaleoToPresentTraitDist) # Add the Min distance to the Agg NEATOMA DATA
                       NEATOMA.Agg2 <- NEATOMA.Agg2[NEATOMA.Agg2$minSCDToPres!=Inf,]
                       
                       sfStop()
                       #Sys.time()-SrtTime
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(NEATOMA.CommDis = NEATOMA.Agg2[,-c(1,7:45)],
                                        ROC.Cutoff = PhyloMPDDist.ROC)
                       Sys.time() -SrtTime1
                       return(Out.list)
                     })
saveRDS(BootPhylo,"./Results/BootPhyloMPD.rds")
                     



# Estimate the Variability in estimates per Trait space iteration 
rm(list=ls());gc()
BootPhylo <- readRDS("./Results/BootPhyloMPD.rds")

sfInit( parallel=TRUE, cpus=10)
sfExport("BootPhylo")

PhyloMPDVar <- sfLapply(BootPhylo,
                         function(x){#x<-BootPhylo[[3]]
                           PhyloMPD <- x$NEATOMA.CommDis
                           MeanMPD <- lapply(1:1000,
                                             function(j){
                                               tapply(PhyloMPD$minSCDToPres[do.call("c",lapply(21:1,function(i){sample(which(PhyloMPD$Time==i),10)}))],
                                                      PhyloMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(PhyloMPD$Time==i),10)}))],
                                                      median)})
                           MeanMPD <- do.call("cbind",
                                              MeanMPD)
                           Out <- data.frame(Time = 1:21,
                                             t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975))))
                           return(Out)
                         })
sfStop()

PhyloMPD <- data.frame(Time = 1:21,
                        X50 = apply(sapply(PhyloMPDVar,function(x){x$X50.}),1,median),
                        X2.5 = apply(sapply(PhyloMPDVar,function(x){x$X2.5.}),1,median),
                        X97.5 = apply(sapply(PhyloMPDVar,function(x){x$X97.5.}),1,median))

plot(y = rev(PhyloMPD$X50),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Phylogentic Change", 
     xlab ="Time (kyrBP)",
     ylab = "MPD (Cophenetic Dist)",
     ylim = range(PhyloMPD[-1]))
lines(y = rev(PhyloMPD$X2.5),
      x = -21:-1)
lines(y = rev(PhyloMPD$X97.5),
      x = -21:-1)

abline(h=median(sapply(BootPhylo,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                       })))

abline(h=quantile(sapply(BootPhylo,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                       }),c(0.025,0.5,0.975)))


#"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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"

