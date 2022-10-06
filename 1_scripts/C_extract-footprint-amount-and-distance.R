library(data.table)
library(dismo)
library(ggplot2)
library(ggspatial)
library(raster)
library(rpart)
library(maptools)
library(rgdal)
theme_set(theme_bw())
library(sf)
library(multidplyr)
library(dplyr)
library(tidyr)
library(parallel)

pkey.sa<-read.csv("0_data/1_processed/point counts/allpkey.studyarea-noABMIgrids.csv", header=TRUE)
str(pkey.sa)
pkey.sa$X.1<-NULL
pkey.sa$X.2<-NULL
#create sf features in "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs" projection
pkey.sf <- st_as_sf(x = pkey.sa,                         
               coords = c("x.alb10tm", "y.alb10tm"),
               crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#human footprint shapefiles extracted from ArcGIS
c.seismic<-st_read("0_data/0_raw/human footprint 2018/ConventionalSeismic.shp")
c.seismic<-st_transform(c.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

pkey.sf.summary<-list()
pkey.sf.demo<-pkey.sf[1:100,]
for (i in 1:nrow(pkey.sf.demo)){
  pkey.sf.i<-pkey.sf.demo[i,]
  PKEY.i<-pkey.sf.i$PKEY
  SS.i<-pkey.sf.i$SS
  YEAR.i<-pkey.sf.i$YEAR
  #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
  c.seismic.y<-c.seismic%>%
    filter(YEAR<YEAR.i)
  #1b. also after time-filtering, create 150-m and 565-m buffers
  b150<-st_buffer(pkey.sf.i,150)
  p.150<-st_intersection(c.seismic.y, b150)
  area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Shape_Area))
  
  b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
  p.565<-st_intersection(c.seismic.y, b565)
  area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Shape_Area))
  # pkey.sf.area[[i]]<-data.frame(PKEY=PKEY.i,
  #                             SS=SS.i,
  #                             YEAR=YEAR.i,
  #                             AREA150=area.150,
  #                             PROP150=area.150/(3.14*150*150),
  #                             AREA565=area.565,
  #                             PROP565=area.565/(3.14*565*565))
  
  #to filter amount of footprint type within 150 and 565 m
  #2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
  maxdist<-1000
  b.maxdist<-st_buffer(pkey.sf.i,maxdist)
  p.maxdist<-st_intersection(c.seismic.y, b.maxdist)
  #3. estimate minimum distance among the time-and-buffer-filtered polygons
  #4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.
  nearest<-ifelse((nrow(p.maxdist)==0), maxdist, 
                  min(st_distance(pkey.sf.i,p.maxdist)))
  pkey.sf.summary[[i]]<-data.frame(PKEY=PKEY.i,
                              SS=SS.i,
                              YEAR=YEAR.i,
                              AREA150=area.150,
                              PROP150=area.150/(3.14*150*150),
                              AREA565=area.565,
                              PROP565=area.565/(3.14*565*565),
                              NEAR.DIST=nearest)
}
summaries.old<-do.call(rbind,pkey.sf.summary)

#wrap in a function
areadist<-function(pointfile, polyfile, maxdist){
  pkey.sf.summary<-list()
  pkey.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pkey.sf.demo)){
    pkey.sf.i<-pkey.sf.demo[i,]
    PKEY.i<-pkey.sf.i$PKEY
    SS.i<-pkey.sf.i$SS
    YEAR.i<-pkey.sf.i$YEAR
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(pkey.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Shape_Area))
    
    b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
    p.565<-st_intersection(polyfile.y, b565)
    area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Shape_Area))

    #to filter amount of footprint type within 150 and 565 m
    #2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
    b.maxdist<-st_buffer(pkey.sf.i,maxdist)
    p.maxdist<-st_intersection(polyfile.y, b.maxdist)
    #3. estimate minimum distance among the time-and-buffer-filtered polygons
    #4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.
    nearest<-ifelse((nrow(p.maxdist)==0), maxdist, 
                    min(st_distance(pkey.sf.i,p.maxdist)))
    pkey.sf.summary[[i]]<-data.frame(PKEY=PKEY.i,
                                     SS=SS.i,
                                     YEAR=YEAR.i,
                                     AREA150=area.150,
                                     PROP150=area.150/(3.14*150*150),
                                     AREA565=area.565,
                                     PROP565=area.565/(3.14*565*565),
                                     NEAR.DIST=nearest)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(rbind,pkey.sf.summary)
  return(summaries)
}

SEISMIC<-areadist(pkey.sf, c.seismic, 1000)#takes ~5 minutes
 
#Parallel processing: dplyr doesn't work with boot or parLapply

#Try parallel processing with tidyverse functions using multiplyr
#https://multidplyr.tidyverse.org/articles/multidplyr.html
#multidplyr doesn't work with sfc objects: needs a vector or ordinary dataframe
parallel::detectCores()#8
#use 1-2 cores less than maximum
cluster <- new_cluster(6)
cluster
pkey.sf.cl <- pkey.sa %>% 
  group_by(PKEY) %>% 
  partition(cluster) %>%
  st_as_sf(coords = c("x.alb10tm", "y.alb10tm"),
           crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# Parallelize any simple features analysis.
#https://stackoverflow.com/questions/48800375/parallelize-st-union-from-rs-sf-package
#the st_par() is someone's custom function: doesn't work for Windows, only Macs
st_par <- function(sf_df, sf_func, n_cores, ...){
  
  # Create a vector to split the data set up by.
  split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))
  
  # Perform GIS analysis
  split_results <- split(sf_df, split_vector) %>%
    mclapply(function(x) sf_func(x), mc.cores = n_cores)
  
  # Combine results back together. Method of combining depends on the output from the function.
  if ( length(class(split_results[[1]]))>1 | class(split_results[[1]])[1] == 'list' ){
    result <- do.call("c", split_results)
    names(result) <- NULL
  } else {
    result <- do.call("rbind", split_results)
  }
  
  # Return result
  return(result)
}

#and parLapply:
no_cores <- detectCores(logical = TRUE)  
cl <- makeCluster(no_cores-2)  
registerDoParallel(cl)
seq_id_all <- seq_along(1:nrow(pkey.sf.demo)) ### intersect_2021 was the file I created with st_intersection. 
areadist.p<-function(){#
  #pkey.sf.summary<-list()
  #pkey.sf.demo<-pointfile[1:100,]
  for (i in seq_id_all(pkey.sf.demo)){
  pkey.sf.i<-pkey.sf[i,]
  PKEY.i<-pkey.sf.i$PKEY
  SS.i<-pkey.sf.i$SS
  YEAR.i<-pkey.sf.i$YEAR
  #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
  polyfile.y<-polyfile%>%mutate(YEAR.i==pkey.sf.i$YEAR)%>%
    filter(YEAR<YEAR.i)
  #1b. also after time-filtering, create 150-m and 565-m buffers
  b150<-st_buffer(pkey.sf.i,150)
  p.150<-st_intersection(polyfile.y, b150)
  area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Shape_Area))
  
  b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
  p.565<-st_intersection(polyfile.y, b565)
  area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Shape_Area))
  
  #to filter amount of footprint type within 150 and 565 m
  #2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
  b.maxdist<-st_buffer(pkey.sf.i,maxdist)
  p.maxdist<-st_intersection(polyfile.y, b.maxdist)
  #3. estimate minimum distance among the time-and-buffer-filtered polygons
  #4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.
  nearest<-ifelse((nrow(p.maxdist)==0), maxdist, 
                  min(st_distance(pkey.sf.i,p.maxdist)))
  pkey.sf.summary<-data.frame(PKEY=PKEY.i,
                              SS=SS.i,
                              YEAR=YEAR.i,
                              AREA150=area.150,
                              PROP150=area.150/(3.14*150*150),
                              AREA565=area.565,
                              PROP565=area.565/(3.14*565*565),
                              NEAR.DIST=nearest)
  return(pkey.sf.summary)
  }
  #summaries<-do.call(rbind,pkey.sf.summary)
  #return(summaries)
}
clusterExport(cl,'pkey.sf.demo')
clusterExport(cl,'c.seismic')
clusterEvalQ(cl,'sf')
clusterEvalQ(cl,'dplyr')
clusterEvalQ(cl,maxdist<-1000)
pkey.sf.summary<-list()
results<- c(parLapply(cl,seq_id_all,fun=areadist.p))
stopCluster(cl)
#Error in checkForRemoteErrors(val) : 
#6 nodes produced errors; first error: unused argument (X[[i]])

SEISMIC2<-st_par(pkey.sf, areadist(pkey.sf, c.seismic, 1000), 4)

source("F:/OSM Good Cumulative Effects/new habitat X edge paper/3_packrat/00-areas-and-distances.R")
# initialising the list to pass to the apply functions
#x.list <- sapply(1:100, list)
nodeslist<-8
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")
clusterExport(cl, "pkey.sf")
clusterExport(cl, "c.seismic")
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, source("F:/OSM Good Cumulative Effects/new habitat X edge paper/3_packrat/00-areas-and-distances.R"))
res<-parLapply(cl, X=1:100, areadist.p(pkey.sf, c.seismic, 1000))
stopCluster(cl)

cat("OK\n* Loading data on master ... ")
pkey.sa<-read.csv("0_data/1_processed/point counts/allpkey.studyarea-noABMIgrids.csv", header=TRUE)
str(pkey.sa)
pkey.sa$X.1<-NULL
pkey.sa$X.2<-NULL
#create sf features in "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs" projection
pkey.sf <- st_as_sf(x = pkey.sa,                         
                    coords = c("x.alb10tm", "y.alb10tm"),
                    crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#human footprint shapefiles extracted from ArcGIS
c.seismic<-st_read("0_data/0_raw/human footprint 2018/ConventionalSeismic.shp")
c.seismic<-st_transform(c.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

cat("OK\nload packages on workers .. .")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, areadist)

cat("OK\n* Exporting and data loading on workers ... ")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
  tmpcl <- clusterEvalQ(cl, setwd("E:/OSM Ideas/new habitat X edge paper/0_data/1_processed/point counts"))
#tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))


cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(YY[,c("REVI","LISP","LEFL")])
SPP <- colnames(YY[,c("AMRE","BBWA","BRCR","BTNW","CAWA","CMWA","OVEN","RCKI","TEWA","WTSP","YRWA")])
DONE <- character(0)
if (interactive() | TEST)
  SPP <- SPP[1:2]

DONE <- substr(list.files(paste0("/out/", PROJ)), 1, 4)
TOGO <- setdiff(SPP, DONE)
mods<-mods_veg#added Feb. 24, 2022

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
while (length(TOGO) > 0) {
  SPP1 <- sample(TOGO, 1)#random sample
  cat("\n  -", length(DONE), "done,", length(TOGO), "more to go, doing", SPP1, "on", date(), "... ")
  if (interactive())
    flush.console()
  t0 <- proc.time()
  #z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegc", ssh_fit="Space")
  if (interactive()) {
    res <- pblapply(cl=cl, X=1:BBB, FUN=run_path1,
                    i=SPP1, mods=mods, CAICalpha=CAICalpha,
                    wcol=NULL, ssh_class=NULL, ssh_fit=NULL)
  } else {
    res <- parLapply(cl, 1:BBB, run_path1,
                     i=SPP1, mods=mods, CAICalpha=CAICalpha,
                     wcol=NULL, ssh_class=NULL, ssh_fit=NULL)
  }
  attr(res, "timing") <- proc.time() - t0
  attr(res, "proj") <- PROJ
  attr(res, "spp") <- SPP1
  attr(res, "CAICalpha") <- CAICalpha
  attr(res, "date") <- as.character(Sys.Date())
  attr(res, "ncl") <- length(cl)
  #attr(res, "conv.seis") <- conv.seismic.bp*max(DAT$`Conventional Seismic`)
  #attr(res, "paved.road") <- road.bp*max(DAT$`Paved Roads`)
  #attr(res, "active.well") <- well.bp*max(DAT$`Active Wells`)
  #attr(res, "pipeline") <- pipeline.bp*max(DAT$Pipelines)
  #attr(res, "facility") <- facility.bp*max(DAT$Industrial)
  save(res,
       file=paste0(getwd(),"/out/", PROJ, "/", SPP1, ".RData"))
  DONE <- substr(list.files(paste0(getwd(),"/out/", PROJ)), 1, 4)
  TOGO <- setdiff(SPP, DONE)
  cat("OK")
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")


seismicXharvest<-st_read("0_data/0_raw/human footprint 2018/seismic_in_harvest2018.shp")
seismicXharvest <- st_set_crs(seismicXharvest, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#get distances
start.time <- Sys.time()
bam.dist_seismicXharvest<-st_distance(bam.studyarea, seismicXharvest, by.element=FALSE)#get the distance matrix
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#Error: cannot allocate vector of size 18.5 Gb
#For now, it looks like I'll have to use the 
#nearest distances I extracted from ArcGIS

#Low Impact Seismic
bam_LIS<-read.csv("0_data/0_raw/distance to footprint/bampoints_distlowimpactseismic.csv", header=TRUE)
bam_LIS$SS<-bam_LIS$SS_V4
bam_LIS$SS_V4<-NULL
bbs_LIS<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distlowimpactseismic.csv", header=TRUE)
bbs_LIS$SS<-bbs_LIS$lctn_nm
bbs_LIS$lctn_nm<-NULL
wt_LIS<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distlowimpactseismic.csv", header=TRUE)
wt_LIS$SS<-wt_LIS$location
wt_LIS$location<-NULL

allpts_LIS<-rbind(bam_LIS,bbs_LIS,wt_LIS)
hist(allpts_LIS$NEAR_DIST)

#Conventional Seismic
bam_CS<-read.csv("0_data/0_raw/distance to footprint/bampoints_distconventionalseismic.csv", header=TRUE)
bam_CS$SS<-bam_CS$SS_V4
bam_CS$SS_V4<-NULL
bbs_CS<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distconventionalseismic.csv", header=TRUE)
bbs_CS$SS<-bbs_CS$lctn_nm
bbs_CS$lctn_nm<-NULL
wt_CS<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distconventionalseismic.csv", header=TRUE)
wt_CS$SS<-wt_CS$location
wt_CS$location<-NULL

allpts_CS<-rbind(bam_CS,bbs_CS,wt_CS)
hist(allpts_CS$NEAR_DIST)

#Pipelines
bam_PL<-read.csv("0_data/0_raw/distance to footprint/bampoints_distpipeline.csv", header=TRUE)
bam_PL$SS<-bam_PL$SS_V4
bam_PL$SS_V4<-NULL
bbs_PL<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distpipeline.csv", header=TRUE)
bbs_PL$SS<-bbs_PL$lctn_nm
bbs_PL$lctn_nm<-NULL
wt_PL<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distpipeline.csv", header=TRUE)
wt_PL$SS<-wt_PL$location
wt_PL$location<-NULL

allpts_PL<-rbind(bam_PL,bbs_PL,wt_PL)
hist(allpts_PL$NEAR_DIST)

#Transmission Lines
bam_TL<-read.csv("0_data/0_raw/distance to footprint/bampoints_distpowerline.csv", header=TRUE)
bam_TL$SS<-bam_TL$SS_V4
bam_TL$SS_V4<-NULL
bbs_TL<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distpowerline.csv", header=TRUE)
bbs_TL$SS<-bbs_TL$lctn_nm
bbs_TL$lctn_nm<-NULL
wt_TL<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distpowerline.csv", header=TRUE)
wt_TL$SS<-wt_TL$location
wt_TL$location<-NULL

allpts_TL<-rbind(bam_TL,bbs_TL,wt_TL)
hist(allpts_TL$NEAR_DIST)

#Abandoned Wells
bam_ABW<-read.csv("0_data/0_raw/distance to footprint/bampoints_distabandonedwell.csv", header=TRUE)
bam_ABW$SS<-bam_ABW$SS_V4
bam_ABW$SS_V4<-NULL
bbs_ABW<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distabandonedwell.csv", header=TRUE)
bbs_ABW$SS<-bbs_ABW$lctn_nm
bbs_ABW$lctn_nm<-NULL
wt_ABW<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distabandonedwell.csv", header=TRUE)
wt_ABW$SS<-wt_ABW$location
wt_ABW$location<-NULL

allpts_ABW<-rbind(bam_ABW,bbs_ABW,wt_ABW)
hist(allpts_ABW$NEAR_DIST)

#Active Wells
bam_ACW<-read.csv("0_data/0_raw/distance to footprint/bampoints_distactivewell.csv", header=TRUE)
bam_ACW$SS<-bam_ACW$SS_V4
bam_ACW$SS_V4<-NULL
bbs_ACW<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distactivewell.csv", header=TRUE)
bbs_ACW$SS<-bbs_ACW$lctn_nm
bbs_ACW$lctn_nm<-NULL
wt_ACW<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distactivewell.csv", header=TRUE)
wt_ACW$SS<-wt_ACW$location
wt_ACW$location<-NULL

allpts_ACW<-rbind(bam_ACW,bbs_ACW,wt_ACW)
hist(allpts_ACW$NEAR_DIST)

#Industrial Facilities
bam_IF<-read.csv("0_data/0_raw/distance to footprint/bampoints_distindustrial.csv", header=TRUE)
bam_IF$SS<-bam_IF$SS_V4
bam_IF$SS_V4<-NULL
bbs_IF<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distindustrial.csv", header=TRUE)
bbs_IF$SS<-bbs_IF$lctn_nm
bbs_IF$lctn_nm<-NULL
wt_IF<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distindustrial.csv", header=TRUE)
wt_IF$SS<-wt_IF$location
wt_IF$location<-NULL

allpts_IF<-rbind(bam_IF,bbs_IF,wt_IF)
hist(allpts_IF$NEAR_DIST)

#Mines
bam_MN<-read.csv("0_data/0_raw/distance to footprint/bampoints_distmines.csv", header=TRUE)
bam_MN$SS<-bam_MN$SS_V4
bam_MN$SS_V4<-NULL
bbs_MN<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distmines.csv", header=TRUE)
bbs_MN$SS<-bbs_MN$lctn_nm
bbs_MN$lctn_nm<-NULL
wt_MN<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distmines.csv", header=TRUE)
wt_MN$SS<-wt_MN$location
wt_MN$location<-NULL

allpts_MN<-rbind(bam_MN,bbs_MN,wt_MN)
hist(allpts_MN$NEAR_DIST)

#Any Road
bam_RD<-read.csv("0_data/0_raw/distance to footprint/bampoints_distanyroad.csv", header=TRUE)
bam_RD$SS<-bam_RD$SS_V4
bam_RD$SS_V4<-NULL
bbs_RD<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distanyroad.csv", header=TRUE)
bbs_RD$SS<-bbs_RD$lctn_nm
bbs_RD$lctn_nm<-NULL
wt_RD<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distanyroad.csv", header=TRUE)
wt_RD$SS<-wt_RD$location
wt_RD$location<-NULL

allpts_RD<-rbind(bam_RD,bbs_RD,wt_RD)
hist(allpts_RD$NEAR_DIST)

#Paved Roads
bam_PVRD<-read.csv("0_data/0_raw/distance to footprint/bampoints_distpavedroad.csv", header=TRUE)
bam_PVRD$SS<-bam_PVRD$SS_V4
bam_PVRD$SS_V4<-NULL
bbs_PVRD<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distpavedroad.csv", header=TRUE)
bbs_PVRD$SS<-bbs_PVRD$lctn_nm
bbs_PVRD$lctn_nm<-NULL
wt_PVRD<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distpavedroad.csv", header=TRUE)
wt_PVRD$SS<-wt_PVRD$location
wt_PVRD$location<-NULL

allpts_PVRD<-rbind(bam_PVRD,bbs_PVRD,wt_PVRD)
hist(allpts_PVRD$NEAR_DIST)

#Gravel Roads
bam_GRRD<-read.csv("0_data/0_raw/distance to footprint/bampoints_distgravelroads.csv", header=TRUE)
bam_GRRD$SS<-bam_GRRD$SS_V4
bam_GRRD$SS_V4<-NULL
bbs_GRRD<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distgravelroads.csv", header=TRUE)
bbs_GRRD$SS<-bbs_GRRD$lctn_nm
bbs_GRRD$lctn_nm<-NULL
wt_GRRD<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distgravelroads.csv", header=TRUE)
wt_GRRD$SS<-wt_GRRD$location
wt_GRRD$location<-NULL

allpts_GRRD<-rbind(bam_GRRD,bbs_GRRD,wt_GRRD)
hist(allpts_GRRD$NEAR_DIST)

#Unimproved Roads
bam_UIRD<-read.csv("0_data/0_raw/distance to footprint/bampoints_distunimprovedroad.csv", header=TRUE)
bam_UIRD$SS<-bam_UIRD$SS_V4
bam_UIRD$SS_V4<-NULL
bbs_UIRD<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distunimprovedroad.csv", header=TRUE)
bbs_UIRD$SS<-bbs_UIRD$lctn_nm
bbs_UIRD$lctn_nm<-NULL
wt_UIRD<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distunimprovedroad.csv", header=TRUE)
wt_UIRD$SS<-wt_UIRD$location
wt_UIRD$location<-NULL

allpts_UIRD<-rbind(bam_UIRD,bbs_UIRD,wt_UIRD)
hist(allpts_UIRD$NEAR_DIST)

#Truck Trails
bam_TT<-read.csv("0_data/0_raw/distance to footprint/bampoints_disttrucktrails.csv", header=TRUE)
bam_TT$SS<-bam_TT$SS_V4
bam_TT$SS_V4<-NULL
bbs_TT<-read.csv("0_data/0_raw/distance to footprint/bbspoints_disttrucktrails.csv", header=TRUE)
bbs_TT$SS<-bbs_TT$lctn_nm
bbs_TT$lctn_nm<-NULL
wt_TT<-read.csv("0_data/0_raw/distance to footprint/wtpoints_disttrucktrails.csv", header=TRUE)
wt_TT$SS<-wt_TT$location
wt_TT$location<-NULL

allpts_TT<-rbind(bam_TT,bbs_TT,wt_TT)
hist(allpts_TT$NEAR_DIST)

#Cropland
bam_CR<-read.csv("0_data/0_raw/distance to footprint/bampoints_distcropland.csv", header=TRUE)
bam_CR$SS<-bam_CR$SS_V4
bam_CR$SS_V4<-NULL
bbs_CR<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distcropland.csv", header=TRUE)
bbs_CR$SS<-bbs_CR$lctn_nm
bbs_CR$lctn_nm<-NULL
wt_CR<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distcropland.csv", header=TRUE)
wt_CR$SS<-wt_CR$location
wt_CR$location<-NULL

allpts_CR<-rbind(bam_CR,bbs_CR,wt_CR)
hist(allpts_CR$NEAR_DIST)

#Harvest
bam_HV<-read.csv("0_data/0_raw/distance to footprint/bampoints_distharvest.csv", header=TRUE)
bam_HV$SS<-bam_HV$SS_V4
bam_HV$SS_V4<-NULL
bbs_HV<-read.csv("0_data/0_raw/distance to footprint/bbspoints_distharvest.csv", header=TRUE)
bbs_HV$SS<-bbs_HV$lctn_nm
bbs_HV$lctn_nm<-NULL
wt_HV<-read.csv("0_data/0_raw/distance to footprint/wtpoints_distharvest.csv", header=TRUE)
wt_HV$SS<-wt_HV$location
wt_HV$location<-NULL

allpts_HV<-rbind(bam_HV,bbs_HV,wt_HV)
hist(allpts_HV$NEAR_DIST)

allpts_ABW$Footprint<-"Abandoned Wells"
allpts_ACW$Footprint<-"Active Wells"
allpts_CR$Footprint<-"Cultivated"
allpts_CS$Footprint<-"Conventional Seismic"
allpts_GRRD$Footprint<-"Gravel Roads"
allpts_HV$Footprint<-"Harvest"
allpts_IF$Footprint<-"Industrial Facilities"
allpts_LIS$Footprint<-"Low Impact Seismic"
allpts_MN$Footprint<-"Mines"
allpts_PL$Footprint<-"Pipelines"
allpts_PVRD$Footprint<-"Paved Roads"
allpts_RD$Footprint<-"Any Road"
allpts_TL$Footprint<-"Transmission Lines"
allpts_TT$Footprint<-"Truck Trails"
allpts_UIRD$Footprint<-"Unimproved Roads"

allpts_allfootprint<-bind_rows(allpts_ABW,
                           allpts_ACW,
                           allpts_CR,
                           allpts_CS,
                           allpts_GRRD,
                           allpts_HV,
                           allpts_IF,
                           allpts_LIS,
                           allpts_MN,
                           allpts_PL,
                           allpts_PVRD,
                           allpts_RD,
                           allpts_TL,
                           allpts_TT,
                           allpts_UIRD)

allpts_allfootprint$Footprint<-as.factor(allpts_allfootprint$Footprint)

p1<-ggplot(allpts_allfootprint, aes(x=NEAR_DIST, fill=Footprint)) + 
  geom_histogram(alpha = 0.5) + 
  theme_bw() + 
  scale_fill_viridis_d() + 
  facet_wrap(~Footprint, scales = "free_x") +
  ylab("Number of Points in Study Area (Upland, Lowland)") +
  xlab("Nearest Distance to Different Footprint")
ggsave(p1, file="0_data/0_raw/distance to footprint/NearestFootprintDistances_Histogram.png", units="in", width=12, height=5)

allpts_allfootprint$FID<-NULL
allpts_allfootprint$NEAR_FID<-NULL
allpts_allfootprint.red<-allpts_allfootprint[,c("SS","Footprint","NEAR_DIST")]
allpts_allfootprint.wide<-allpts_allfootprint.red%>%
  pivot_wider(names_from = Footprint,
              values_from = NEAR_DIST)
str(allpts_allfootprint.wide)#18877x16

harvestdata<-allpts_HV[,c("SS","YEAR")]
harvestdata$HarvestYear<-harvestdata$YEAR
harvestdata$YEAR<-NULL

allpts_allfootprint.wideB<-merge(allpts_allfootprint.wide, harvestdata, by=c("SS"))

save(allpts_allfootprint.wideB, file="0_data/1_processed/point counts/allpoints_footprintdistance.RData")

#add footprint to habitat data then filter point counts
load("0_data/1_processed/point counts/allpoints_footprintdistance.RData")
load("0_data/1_processed/point counts/bampoints.studyarea.habitatdata.RData")
load("0_data/1_processed/point counts/bbspoints.studyarea.habitatdata.RData")
load("0_data/1_processed/point counts/wtpoints.studyarea.habitatdata.RData")

names(bam.data)

allpoints.habitat<-rbind(bam.data,bbs.data,wt.data)
names(allpoints.habitat)
str(allpoints.habitat)#90927 obs. of  283 variables

allpoints.habitat.footprint<-merge(allpoints.habitat,
                                   allpts_allfootprint.wideB,
                                   by=c("SS"))
str(allpoints.habitat.footprint)#90927 obs. of  299 variables


#Filter to upland only
nrow(allpoints.habitat.footprint[(allpoints.habitat.footprint$Species_Abie_Ama_v1+
                                    allpoints.habitat.footprint$Species_Abie_Bal_v1+
                                    allpoints.habitat.footprint$Species_Abie_Las_v1+
                                    allpoints.habitat.footprint$Species_Popu_Bal_v1+
                                    allpoints.habitat.footprint$Species_Popu_Tre_v1+
                                    allpoints.habitat.footprint$Species_Pinu_Ban_v1+
                                    allpoints.habitat.footprint$Species_Pinu_Con_v1+
                                    allpoints.habitat.footprint$Species_Pice_Gla_v1)>50,])#58985 visits

upland<-allpoints.habitat.footprint[(allpoints.habitat.footprint$Species_Abie_Ama_v1+
                                       allpoints.habitat.footprint$Species_Abie_Bal_v1+
                                       allpoints.habitat.footprint$Species_Abie_Las_v1+
                                       allpoints.habitat.footprint$Species_Popu_Bal_v1+
                                       allpoints.habitat.footprint$Species_Popu_Tre_v1+
                                       allpoints.habitat.footprint$Species_Pinu_Ban_v1+
                                       allpoints.habitat.footprint$Species_Pinu_Con_v1+
                                       allpoints.habitat.footprint$Species_Pice_Gla_v1)>50,]

#Filter to deciduous and deciduous-leaning upland forest
decid.lean<-upland[upland$SpeciesGroups_Broadleaf_Spp_v1>50,]
nrow(decid.lean)#46719 visits

#Filter to coniferous and coniferous-leaning upland forest
conif.lean<-upland[upland$SpeciesGroups_Needleleaf_Spp_v1>50,]
nrow(conif.lean)#11145

#several hundred upland points have forest cover >50% but 
#deciduous and coniferous cover are each <50%
conif.long<-conif.lean[,c("Abandoned Wells",
                          "Active Wells",
                          "Conventional Seismic",
                          "Cultivated",#check that you're not replacing one of the raster predictors
                          "Gravel Roads",
                          "Harvest",
                          "Industrial Facilities",
                          "Low Impact Seismic",
                          "Mines",
                          "Pipelines",
                          "Paved Roads",
                          "Any Road",
                          "Transmission Lines",
                          "Truck Trails",
                          "Unimproved Roads")]%>%
  pivot_longer(cols=`Abandoned Wells`:`Unimproved Roads`,
               names_to="Footprint", 
               values_to="NEAR_DIST")

nrow(conif.long)#167175
p2<-ggplot(conif.long, aes(x=NEAR_DIST, fill=Footprint)) + 
  geom_histogram(alpha = 0.5) + 
  theme_bw() + 
  scale_fill_viridis_d() + 
  facet_wrap(~Footprint, scales = "free_x") +
  ylab("Number of Coniferous Upland Points") +
  xlab("Nearest Distance to Different Footprint")
ggsave(p2, file="0_data/0_raw/distance to footprint/NearestFootprintDistances_Histogram_ConiferousUplands.png", units="in", width=12, height=5)

decid.long<-decid.lean[,c("Abandoned Wells",
                          "Active Wells",
                          "Conventional Seismic",
                          "Cultivated",
                          "Gravel Roads",
                          "Harvest",
                          "Industrial Facilities",
                          "Low Impact Seismic",
                          "Mines",
                          "Pipelines",
                          "Paved Roads",
                          "Any Road",
                          "Transmission Lines",
                          "Truck Trails",
                          "Unimproved Roads")]%>%
  pivot_longer(cols=`Abandoned Wells`:`Unimproved Roads`,
               names_to="Footprint", 
               values_to="NEAR_DIST")

nrow(decid.long)
p3<-ggplot(decid.long, aes(x=NEAR_DIST, fill=Footprint)) + 
  geom_histogram(alpha = 0.5) + 
  theme_bw() + 
  scale_fill_viridis_d() + 
  facet_wrap(~Footprint, scales = "free_x") +
  ylab("Number of Deciduous Upland Points") +
  xlab("Nearest Distance to Different Footprint")
ggsave(p3, file="0_data/0_raw/distance to footprint/NearestFootprintDistances_Histogram_DeciduousUplands.png", units="in", width=12, height=5)

#histograms of footprint distance from coniferous
#and deciduous forest point counts look similar. Coniferous
#point counts skewed to be closer on average to some kinds of
#footprint.