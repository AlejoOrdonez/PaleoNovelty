rm(list=ls());gc()
require(analogue)
require(FD)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/Data")
#######################################################################################################################
#####
# Load the North American Pollen Data Base and its translation to the Paleo data
# raw NPDB
BootTraits <- lapply(1:100,
                     function(RunIt){
                       SrtTime <- Sys.time()
                       NPDB <- read.csv("./Pollen/NorthAmerican Pollen DBS/whitmoreetal2005_v1-72.csv")
                       # Taxa translate between NPDB and Paleo Data
                       TaxtTrans <- read.csv("./Pollen/NorthAmerican Pollen DBS/Taxa transation.csv")
                       # Which are the summarizing taxa - note the TaxtTrans$With.Traits!="" filter -> only species with traits are used
                       TaxaUse <- sort(unique(TaxtTrans$Translated[TaxtTrans$With.Traits!=""]))
                       #TaxaUse <- sort(unique(TaxtTrans$Translated))[-1]
                       #####
                       
                       #####
                       # Build a translated NPDB dataset
                       SummNPDB <- lapply(TaxaUse, # only use species with trait data (Woody species)
                                          function(x){
                                            out <- NPDB[,TaxtTrans$Original[TaxtTrans$Translated==x]]
                                            if(class(out)=="data.frame"){out <- apply(out,1,sum)}
                                            out <- as.numeric(out)
                                            return(out)
                                          })
                       names(SummNPDB) <- TaxaUse
                       SummNPDB2 <- as.data.frame(do.call("cbind",SummNPDB))
                       SummNPDB3 <- as.data.frame(t(apply(SummNPDB2,1,function(x){x/sum(x,na.rm=T)})))
                       SummNPDB3[is.na(SummNPDB3)] <- 0
                       #####
                       
                       #####
                       # Load the trait data an turnit into a distance matrix
                       Traits <- read.csv("./Traits/TraitsSummary.csv")
                       Traits$TaxtTrans <- TaxtTrans$Translated[match(Traits$Genus,TaxtTrans$PLANTS)]
                       # Seeds 
                       SWT <- apply(Traits[,c("Min..of.Seed.Size..mg.","Max..of.Seed.Size..mg.")],
                                    1,
                                    function(x){
                                      runif(1,x[1],x[2])})
                       SWT <- scale(log10(SWT))
                       
                       # Height
                       Hmax <- apply(Traits[,c("Min..of.Height..m.","Max..of.Height..m.")],
                                     1,
                                     function(x){
                                       runif(1,x[1],x[2])})
                       Hmax <- scale(Hmax)
                       # LMA
                       LMA <- apply(Traits[,c("Min..of.LMA..log...g.m2.","Max..of.LMA..log...g.m2.")],
                                    1,
                                    function(x){
                                      runif(1,x[1],x[2])})
                       LMA <- scale(LMA)
                       
                       # Build a trait dataset
                       TraitsDBS <- data.frame(SWT = SWT,
                                               Hmax = Hmax,
                                               LMA = LMA)
                       row.names(TraitsDBS) <- Traits$TaxtTrans
                       TraitsDBS <- TraitsDBS[complete.cases(TraitsDBS),]
                       
                       #turn traits into PCoA axes
                       GowdisDfrm <- gowdis(TraitsDBS)
                       PCoATraits <- ape::pcoa(GowdisDfrm)
                       PCoATraits <- PCoATraits $ vectors
                       #####
                       
                       #####
                       # Build a occurrence dataset
                       SummNPDB4 <- SummNPDB3[,rownames(TraitsDBS)]
                       SummNPDB4 <- SummNPDB4[apply(SummNPDB4,1,sum)!=0,]
                       SummNPDB4 <- SummNPDB4[c(apply(SummNPDB4,1,function(x){sum(x!=0)})>3),]
                       SummNPDB4 <- as.data.frame(t(apply(SummNPDB4,1,function(x){x/sum(x,na.rm=T)})))
                       
                       # Build a species names per site data set
                       SppPerSite <- lapply(1:dim(SummNPDB4)[1],
                                            function(x){
                                              names(SummNPDB4)[SummNPDB4[x,]!=0]
                                            })
                       
                       #####
                       #####
                       # Estimate the Trait distance in PCoA Space                       
                       PCoATraitsDist <- as.matrix(dist(PCoATraits))
                       
                       # Trait Mean pairwise distance contrast
                       MPDBaseline <- lapply(1:length(SppPerSite),
                                             function(i){
                                               sapply(SppPerSite,
                                                      function(x){
                                                        mean(PCoATraitsDist[SppPerSite[[i]],x])
                                                      })
                                             })
                       
                       MPDBaseline <- do.call("rbind",MPDBaseline)
                       #####
                       
                       #####
                       # Define the analogue threshold
                       BIOMES <- vect("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/Data/WWF-Biomes/wwf_terr_ecos.shp")
                       BIomesNames <- data.frame(code = 1:14,
                                                 Name = c("Tropical and Subtropical Moist Broadleaf Forests",
                                                          "Tropical and Subtropical Dry Broadleaf Forests",
                                                          "Tropical and Subtropical Coniferous Forests",
                                                          "Temperate Broadleaf and Mixed Forests",
                                                          "Temperate Coniferous Forests",
                                                          "Boreal forests/Taiga",
                                                          "Tropical and Subtropical Grasslands, Savannas, and Shrublands",
                                                          "Temperate Grasslands, Savannas, and Shrublands",
                                                          "Flooded Grasslands and Savannas",
                                                          "Montane Grasslands and Shrublands",
                                                          "Tundra",
                                                          "Mediterranean Forests, Woodlands, and Scrub",
                                                          "Deserts and Xeric Shrublands",
                                                          "Mangroves"))
                       # Get the BIOME type for each NPDB location
                       NPDBLoc <- NPDB[rownames(SummNPDB4),c("LONDD","LATDD")]
                       names(NPDBLoc) <- c('x',"y")
                       NPDBBiome <- extract(BIOMES[,"BIOME"],
                                            NPDBLoc)
                       # Make NA values a category
                       NPDBBiome[is.na(NPDBBiome[,2]),2] <- 100
                       
                       
                       # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                       TraitMPDDist <- MPDBaseline[!NPDBBiome[,2]%in%c(3,98,99,100),!NPDBBiome[,2]%in%c(3,98,99,100)]
                       TraitMPDDist.ROC <- roc(TraitMPDDist, # current time Taxon data.frame 
                                               groups = BIomesNames$Name[NPDBBiome[!NPDBBiome[,2]%in%c(3,98,99,100),2]] # vector of group memberships
                       )
                       TraitMPDDist.ROC
                       #####
                       
                       #####
                       # Load the Paleo-Pollen DBS
                       Old <- read.csv("./Pollen/Ordonez & Williams 2013/abund_spp_env_21to1bp.csv")
                       Old2 <- Old[,rownames(PCoATraits)]
                       Old2 <- as.data.frame(t(apply(Old2,1,function(x){x/sum(x,na.rm=T)})))
                       Old2 <- data.frame(Old[, c("Time","Clim.ID","sites","Latitude","Longitude")],
                                          Old2)
                       Old2 <- Old2[c(apply(Old2[,rownames(PCoATraits)],1,function(x){sum(x!=0)})>3),]
                       Old2 <- na.omit(Old2)
                       # Build a species per site Data
                       PollenSppSiteNames <- lapply(1:dim(Old2)[1],
                                                    function(x){
                                                      names(SummNPDB4)[Old2[x,names(SummNPDB4)]!=0]
                                                    })
                       #####

                       #####
                       # Estimate distance to each paleo assemblage
                       PaleoToPresentTraitDist <- lapply(PollenSppSiteNames,
                                                         function(i){
                                                           out <- sapply(SppPerSite,
                                                                         function(x){
                                                                           mean(PCoATraitsDist[i,x],na.rm=T)
                                                                         })
                                                           min(out)
                                                         })
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(DataFrm = data.frame(Time=-21:-1,
                                                             MPDDist = rev(tapply(do.call("c",PaleoToPresentTraitDist),Old2$Time,median))),
                                        RocTresh = TraitMPDDist.ROC)
                       Sys.time()-SrtTime
                       return(Out.list)
                     })
saveRDS(BootTraits,"./Traits/BootPCoATraitsMPD.rds")
