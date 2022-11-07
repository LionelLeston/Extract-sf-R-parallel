# library(data.table)
# library(dismo)
# library(ggplot2)
# library(ggspatial)
# library(raster)
# library(rpart)
# library(maptools)
# library(rgdal)
# theme_set(theme_bw())
library(sf)
library(multidplyr)
library(dplyr)
library(tidyr)
library(parallel)
#I think units, then spatstat is loaded with sf or another package
#when loaded after units, spatstat's own units function supersedes the
#creation of "units" objects, preventing weighted.mean function from working
pkey.sa<-read.csv("0_data/1_processed/point counts/allpkey.studyarea-noABMIgrids.csv", header=TRUE)
str(pkey.sa)#83439
pkey.sa$X.1<-NULL
pkey.sa$X.2<-NULL
#create sf features in "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs" projection
pkey.sf <- st_as_sf(x = pkey.sa,                         
               coords = c("x.alb10tm", "y.alb10tm"),
               crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#human footprint shapefiles extracted from ArcGIS
#hf2018<-st_read("0_data/0_raw/human footprint 2018/hf2018extent.shp")
#hf2018<-st_transform(hf2018, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#conventional seismic lines
c.seismic<-st_read("0_data/0_raw/human footprint 2018/ConventionalSeismic.shp")
c.seismic<-st_transform(c.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# pkey.sf.summary<-list()
# pkey.sf.demo<-pkey.sf[1:100,]
# for (i in 1:nrow(pkey.sf.demo)){
#   pkey.sf.i<-pkey.sf.demo[i,]
#   PKEY.i<-pkey.sf.i$PKEY
#   SS.i<-pkey.sf.i$SS
#   YEAR.i<-pkey.sf.i$YEAR
#   #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
#   c.seismic.y<-c.seismic%>%
#     filter(YEAR<YEAR.i)
#   #1b. also after time-filtering, create 150-m and 565-m buffers
# p.150<-st_intersection(polyfile.y, b150)
# p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
# area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Recalc_Area))
# 
# b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
# p.565<-st_intersection(polyfile.y, b565)
# p.565$Recalc_Area<-st_area(p.565)#recalculate areas of clipped polygons
# area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Recalc_Area))

#   # pkey.sf.area[[i]]<-data.frame(PKEY=PKEY.i,
#   #                             SS=SS.i,
#   #                             YEAR=YEAR.i,
#   #                             AREA150=area.150,
#   #                             PROP150=area.150/(3.14*150*150),
#   #                             AREA565=area.565,
#   #                             PROP565=area.565/(3.14*565*565))
#   
#   #to filter amount of footprint type within 150 and 565 m
#   #2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
#   maxdist<-1000
#   b.maxdist<-st_buffer(pkey.sf.i,maxdist)
#   p.maxdist<-st_intersection(c.seismic.y, b.maxdist)
#   #3. estimate minimum distance among the time-and-buffer-filtered polygons
#   #4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.
#   nearest<-ifelse((nrow(p.maxdist)==0), maxdist, 
#                   min(st_distance(pkey.sf.i,p.maxdist)))
#   pkey.sf.summary[[i]]<-data.frame(PKEY=PKEY.i,
#                               SS=SS.i,
#                               YEAR=YEAR.i,
#                               AREA150=area.150,
#                               PROP150=area.150/(3.14*150*150),
#                               AREA565=area.565,
#                               PROP565=area.565/(3.14*565*565),
#                               NEAR.DIST=nearest)
# }
# summaries.old<-do.call(rbind,pkey.sf.summary)

#wrap in a function
areadist<-function(pointfile, polyfile, maxdist){
  pkey.sf.summary<-list()
  #pkey.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    pkey.sf.i<-pointfile[i,]#pkey.sf.demo[i,]
    PKEY.i<-pkey.sf.i$PKEY
    SS.i<-pkey.sf.i$SS
    YEAR.i<-pkey.sf.i$YEAR
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(pkey.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection

    #area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Recalc_Area))
    
    b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
    p.565<-st_intersection(polyfile.y, b565)
    p.565$Recalc_Area<-st_area(p.565)#recalculate areas of clipped polygons
    area.565 <- ifelse((nrow(p.565)==0), 0, (p.565 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    #area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Recalc_Area))

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

#Modified areadist function for harvest blocks, where age of harvests within 150 m and 565 m is calculated
#Might also be applied to other footprint if age of those footprints is of interest
#Note weighted.mean function not used because it didn't work with sf package
areadist.age<-function(pointfile, polyfile, maxdist, oldestyear){
  pkey.sf.summary<-list()
  #pkey.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    pkey.sf.i<-pointfile[i,]#pkey.sf.demo[i,]
    PKEY.i<-pkey.sf.i$PKEY
    SS.i<-pkey.sf.i$SS
    YEAR.i<-pkey.sf.i$YEAR
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)#should this filtering occur on the polygons after isolating those within 150 and 565 m of each point count?
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(pkey.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    p.150$YEAR<-ifelse(p.150$YEAR==0, oldestyear, p.150$YEAR)
    p.150$Age<-YEAR.i-p.150$YEAR #ifelse((nrow(p.150))==0,NA,(YEAR.i-p.150$YEAR))   #will be null if nrow(p.150)==0
    p.150$AgeXArea<-p.150$Age*p.150$Recalc_Area#ifelse((nrow(p.150))==0,NA,(p.150$Age*p.150$Recalc_Area))#will be null if nrow(p.150)==0
    #calculate area amount after dissolving polygons
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    sum.recalc.area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Recalc_Area))
    
    SumAgeXArea150<-sum(p.150$AgeXArea)#will be null if nrow(p.150)==0
    MeanWtAge.150<-ifelse((nrow(p.150)==0), NA, SumAgeXArea150/sum.recalc.area.150)#area.150
    #MeanWtHarvestAge will be used to modify WtForestAge if not NA
    b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
    p.565<-st_intersection(polyfile.y, b565)
    p.565$Recalc_Area<-st_area(p.565)#recalculate areas of clipped polygons
    #calculate area amount after dissolving polygons
    area.565 <- ifelse((nrow(p.565)==0), 0, (p.565 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    sum.recalc.area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Recalc_Area))
    p.565$YEAR<-ifelse(p.565$YEAR==0, oldestyear, p.565$YEAR)
    p.565$Age<-YEAR.i-p.565$YEAR   #will be null if nrow(p.565)==0
    p.565$AgeXArea<-p.565$Age*p.565$Recalc_Area #will be null if nrow(p.565)==0
    SumAgeXArea565<-sum(p.565$AgeXArea)#will be null if nrow(p.565)==0
    MeanWtAge.565<-ifelse((nrow(p.565)==0), NA, SumAgeXArea565/sum.recalc.area.565)#area.565
    #MeanWtHarvestAge will be used to modify WtForestAge if not NA
    
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
                                     MEANAGE.150=MeanWtAge.150,
                                     AREA565=area.565,
                                     PROP565=area.565/(3.14*565*565),
                                     MEANAGE.565=MeanWtAge.565,
                                     NEAR.DIST=nearest)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(rbind,pkey.sf.summary)
  return(summaries)
}

#SEISMIC<-areadist(pkey.sf, c.seismic, 1000)#takes ~5 minutes
 
#Try splitting the polygon features shapefile into smaller chunks.
str(c.seismic)
#Bounding box:  xmin: 317461.7 ymin: 5936908 xmax: 830294.2 ymax: 6431509
#CRS:           +proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs

#https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r
#create a global object with same extent as human footprint layer
#(actually used conventional seismic line since its extent is a bit bigger)
hf2018_bb <- matrix(c(317461.7, 6431509,
                      830294.2, 6431509,
                      830294.2, 5936908,
                      317461.7, 5936908,
                      317461.7, 6431509), byrow = TRUE, ncol = 2) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(., crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#then make a grid whose cells (~100x100 km) can be used to intersect shapefile and make tiles
hf2018_grid_10x10 <- st_make_grid(hf2018_bb, n = c(100, 100), 
                                 crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs", 
                                 what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

#now go through each cell in the grid and use it to intersect part of the
#the point counts. For each group of point counts, generate a buffer 
#equal to maxdist and use this multiple-point-counts-buffer to intersect
#footprint layer; then use this final intersection and point-counts-group
#in the areadist function
#pointcountgroup<-list()
maxdist<-1000
oldestyear<-min(c.seismic$YEAR[c.seismic$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(c.seismic, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/conventional seismic areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Low Impact Seismic Lines
li.seismic<-st_read("0_data/0_raw/human footprint 2018/LowImpactSeismic.shp")
li.seismic<-st_transform(li.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-1000
oldestyear<-min(li.seismic$YEAR[li.seismic$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(li.seismic, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/low impact seismic areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Paved Roads
paved.roads<-st_read("0_data/0_raw/human footprint 2018/PavedRoads.shp")
paved.roads<-st_transform(paved.roads, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-1000
oldestyear<-min(paved.roads$YEAR[paved.roads$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(paved.roads, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/paved road areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Gravel Roads
gravel.roads<-st_read("0_data/0_raw/human footprint 2018/GravelRoads.shp")
gravel.roads<-st_transform(gravel.roads, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(gravel.roads[st_is_valid(gravel.roads)==FALSE,])#19
gravel.roads<-gravel.roads[st_is_valid(gravel.roads)==TRUE,]

maxdist<-1000
oldestyear<-min(gravel.roads$YEAR[gravel.roads$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(gravel.roads, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/gravel road areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Unimproved Roads
unimproved.roads<-st_read("0_data/0_raw/human footprint 2018/UnimprovedRoads.shp")
unimproved.roads<-st_transform(unimproved.roads, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(unimproved.roads[st_is_valid(unimproved.roads)==FALSE,])#5
unimproved.roads<-unimproved.roads[st_is_valid(unimproved.roads)==TRUE,]

maxdist<-1000
oldestyear<-min(unimproved.roads$YEAR[unimproved.roads$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(unimproved.roads, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/unimproved road areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#TruckTrails
truck.trails<-st_read("0_data/0_raw/human footprint 2018/TruckTrails.shp")
truck.trails<-st_transform(truck.trails, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(truck.trails[st_is_valid(truck.trails)==FALSE,])#0
truck.trails<-truck.trails[st_is_valid(truck.trails)==TRUE,]

maxdist<-1000
oldestyear<-min(truck.trails$YEAR[truck.trails$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(truck.trails, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/truck trail areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Pipelines
pipelines<-st_read("0_data/0_raw/human footprint 2018/pipelines2018layer.shp")
pipelines<-st_transform(pipelines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#check for corrupt geometries
nrow(pipelines[st_is_valid(pipelines)==FALSE,])#49
pipelines<-pipelines[st_is_valid(pipelines)==TRUE,]

maxdist<-1000
oldestyear<-min(pipelines$YEAR[pipelines$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(pipelines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/pipeline areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Transmission Lines
transmissionlines<-st_read("0_data/0_raw/human footprint 2018/transmissionlines2018layer.shp")
transmissionlines<-st_transform(transmissionlines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(transmissionlines[st_is_valid(transmissionlines)==FALSE,])#0
transmissionlines<-transmissionlines[st_is_valid(transmissionlines)==TRUE,]

maxdist<-1000
oldestyear<-min(transmissionlines$YEAR[transmissionlines$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(transmissionlines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/transmission line areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Active Wells
activewells<-st_read("0_data/0_raw/human footprint 2018/activewells2018layer.shp")
activewells<-st_transform(activewells, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(activewells[st_is_valid(activewells)==FALSE,])#0
activewells<-activewells[st_is_valid(activewells)==TRUE,]

maxdist<-1000
oldestyear<-min(activewells$YEAR[activewells$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(activewells, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/active well areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Abandoned Wells
abandonedwells<-st_read("0_data/0_raw/human footprint 2018/abandonedwells2018layer.shp")
abandonedwells<-st_transform(abandonedwells, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(abandonedwells[st_is_valid(abandonedwells)==FALSE,])#0
abandonedwells<-abandonedwells[st_is_valid(abandonedwells)==TRUE,]

maxdist<-1000
oldestyear<-min(abandonedwells$YEAR[abandonedwells$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(abandonedwells, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/abandoned well areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Harvest
harvest<-st_read("0_data/0_raw/human footprint 2018/harvestareas2018layer.shp")
harvest<-st_transform(harvest, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(harvest[st_is_valid(harvest)==FALSE,])#61
harvest<-harvest[st_is_valid(harvest)==TRUE,]

maxdist<-1000
oldestyear<-min(harvest$YEAR[harvest$YEAR>0])
#Uses areadist.age function instead of areadist function
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(harvest, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/harvest areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Industrial
industrial<-st_read("0_data/0_raw/human footprint 2018/industry2018layer.shp")
industrial<-st_transform(industrial, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(industrial[st_is_valid(industrial)==FALSE,])#15
industrial<-industrial[st_is_valid(industrial)==TRUE,]

maxdist<-1000
oldestyear<-min(industrial$YEAR[industrial$YEAR>0])
#Uses areadist.age function instead of areadist function
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(industrial, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/industrial areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Residential
residential<-st_read("0_data/0_raw/human footprint 2018/residential2018layer.shp")
residential<-st_transform(residential, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(residential[st_is_valid(residential)==FALSE,])#8
residential<-residential[st_is_valid(residential)==TRUE,]

maxdist<-1000
oldestyear<-min(residential$YEAR[residential$YEAR>0])
#Uses areadist.age function instead of areadist function
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(residential, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/residential areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Mines
mines<-st_read("0_data/0_raw/human footprint 2018/mines2018layer.shp")
mines<-st_transform(mines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(mines[st_is_valid(mines)==FALSE,])#12
mines<-mines[st_is_valid(mines)==TRUE,]

maxdist<-1000
oldestyear<-min(mines$YEAR[mines$YEAR>0])#1950
#Uses areadist.age function instead of areadist function
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(mines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/mine areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Cultivated Lands (in extreme southeast and west there is a lot)
cultivation<-st_read("0_data/0_raw/human footprint 2018/cultivation2018layer.shp")
cultivation<-st_transform(cultivation, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#corrupt geometries detected the first time loops were run for this file
nrow(cultivation[st_is_valid(cultivation)==FALSE,])#91
cultivation<-cultivation[st_is_valid(cultivation)==TRUE,]

maxdist<-1000
oldestyear<-min(cultivation$YEAR[cultivation$YEAR>0])#1937
#Uses areadist.age function instead of areadist function
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(pkey.sf, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(cultivation, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/cultivation areadist.age/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Now for each footprint, combine the smaller files into one large file
#for each footprint type's area, weighted mean age, and nearest distance
#to each point count survey. There should be 83439 observations in each 
#footprint file.

#conventional seismic
names1<-list.files(path = "0_data/1_processed/human footprint/conventional seismic areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/conventional seismic areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/conventional-seismic-areadist.age-allsurveys.csv")

#low impact seismic
names1<-list.files(path = "0_data/1_processed/human footprint/low impact seismic areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/low impact seismic areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/low-impact-seismic-areadist.age-allsurveys.csv")

#pipelines
names1<-list.files(path = "0_data/1_processed/human footprint/pipeline areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/pipeline areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/pipeline-areadist.age-allsurveys.csv")

#transmission lines
names1<-list.files(path = "0_data/1_processed/human footprint/transmission line areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/transmission line areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/transmission-line-areadist.age-allsurveys.csv")

#paved roads
names1<-list.files(path = "0_data/1_processed/human footprint/paved road areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/paved road areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/paved-road-areadist.age-allsurveys.csv")

#gravel roads
names1<-list.files(path = "0_data/1_processed/human footprint/gravel road areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/gravel road areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/gravel-road-areadist.age-allsurveys.csv")

#unimproved roads
names1<-list.files(path = "0_data/1_processed/human footprint/unimproved road areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/unimproved road areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/unimproved-road-areadist.age-allsurveys.csv")

#truck trails
names1<-list.files(path = "0_data/1_processed/human footprint/truck trail areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/truck trail areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/truck-trail-areadist.age-allsurveys.csv")

#active wells
names1<-list.files(path = "0_data/1_processed/human footprint/active well areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/active well areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/active-well-areadist.age-allsurveys.csv")

#abandoned wells
names1<-list.files(path = "0_data/1_processed/human footprint/abandoned well areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/abandoned well areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/abandoned-well-areadist.age-allsurveys.csv")

#industrial facilities
names1<-list.files(path = "0_data/1_processed/human footprint/industrial areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/industrial areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/industrial-areadist.age-allsurveys.csv")

#mines
names1<-list.files(path = "0_data/1_processed/human footprint/mine areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/mine areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/mine-areadist.age-allsurveys.csv")

#residential
names1<-list.files(path = "0_data/1_processed/human footprint/residential areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/residential areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/residential-areadist.age-allsurveys.csv")

#cultivation
names1<-list.files(path = "0_data/1_processed/human footprint/cultivation areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/cultivation areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/cultivation-areadist.age-allsurveys.csv")

#harvest
names1<-list.files(path = "0_data/1_processed/human footprint/harvest areadist.age/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/human footprint/harvest areadist.age/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/human footprint/harvest-areadist.age-allsurveys.csv")


#merge files
conventional.seismic<-read.csv("0_data/1_processed/human footprint/conventional-seismic-areadist.age-allsurveys.csv", header=TRUE)
str(conventional.seismic)
conventional.seismic$AREA150.conventional.seismic<-conventional.seismic$AREA150
conventional.seismic$AREA150<-NULL
conventional.seismic$PROP150.conventional.seismic<-conventional.seismic$PROP150
conventional.seismic$PROP150<-NULL
conventional.seismic$MEANAGE.150.conventional.seismic<-conventional.seismic$MEANAGE.150
conventional.seismic$MEANAGE.150<-NULL
conventional.seismic$AREA565.conventional.seismic<-conventional.seismic$AREA565
conventional.seismic$AREA565<-NULL
conventional.seismic$PROP565.conventional.seismic<-conventional.seismic$PROP565
conventional.seismic$PROP565<-NULL
conventional.seismic$MEANAGE.565.conventional.seismic<-conventional.seismic$MEANAGE.565
conventional.seismic$MEANAGE.565<-NULL
conventional.seismic$NEAR.DIST.conventional.seismic<-conventional.seismic$NEAR.DIST
conventional.seismic$NEAR.DIST<-NULL
str(conventional.seismic)

low.impact.seismic<-read.csv("0_data/1_processed/human footprint/low-impact-seismic-areadist.age-allsurveys.csv", header=TRUE)
str(low.impact.seismic)
low.impact.seismic$AREA150.low.impact.seismic<-low.impact.seismic$AREA150
low.impact.seismic$AREA150<-NULL
low.impact.seismic$PROP150.low.impact.seismic<-low.impact.seismic$PROP150
low.impact.seismic$PROP150<-NULL
low.impact.seismic$MEANAGE.150.low.impact.seismic<-low.impact.seismic$MEANAGE.150
low.impact.seismic$MEANAGE.150<-NULL
low.impact.seismic$AREA565.low.impact.seismic<-low.impact.seismic$AREA565
low.impact.seismic$AREA565<-NULL
low.impact.seismic$PROP565.low.impact.seismic<-low.impact.seismic$PROP565
low.impact.seismic$PROP565<-NULL
low.impact.seismic$MEANAGE.565.low.impact.seismic<-low.impact.seismic$MEANAGE.565
low.impact.seismic$MEANAGE.565<-NULL
low.impact.seismic$NEAR.DIST.low.impact.seismic<-low.impact.seismic$NEAR.DIST
low.impact.seismic$NEAR.DIST<-NULL
str(low.impact.seismic)

pipeline<-read.csv("0_data/1_processed/human footprint/pipeline-areadist.age-allsurveys.csv", header=TRUE)
str(pipeline)
pipeline$AREA150.pipeline<-pipeline$AREA150
pipeline$AREA150<-NULL
pipeline$PROP150.pipeline<-pipeline$PROP150
pipeline$PROP150<-NULL
pipeline$MEANAGE.150.pipeline<-pipeline$MEANAGE.150
pipeline$MEANAGE.150<-NULL
pipeline$AREA565.pipeline<-pipeline$AREA565
pipeline$AREA565<-NULL
pipeline$PROP565.pipeline<-pipeline$PROP565
pipeline$PROP565<-NULL
pipeline$MEANAGE.565.pipeline<-pipeline$MEANAGE.565
pipeline$MEANAGE.565<-NULL
pipeline$NEAR.DIST.pipeline<-pipeline$NEAR.DIST
pipeline$NEAR.DIST<-NULL
str(pipeline)

transmission.line<-read.csv("0_data/1_processed/human footprint/transmission-line-areadist.age-allsurveys.csv", header=TRUE)
str(transmission.line)
transmission.line$AREA150.transmission.line<-transmission.line$AREA150
transmission.line$AREA150<-NULL
transmission.line$PROP150.transmission.line<-transmission.line$PROP150
transmission.line$PROP150<-NULL
transmission.line$MEANAGE.150.transmission.line<-transmission.line$MEANAGE.150
transmission.line$MEANAGE.150<-NULL
transmission.line$AREA565.transmission.line<-transmission.line$AREA565
transmission.line$AREA565<-NULL
transmission.line$PROP565.transmission.line<-transmission.line$PROP565
transmission.line$PROP565<-NULL
transmission.line$MEANAGE.565.transmission.line<-transmission.line$MEANAGE.565
transmission.line$MEANAGE.565<-NULL
transmission.line$NEAR.DIST.transmission.line<-transmission.line$NEAR.DIST
transmission.line$NEAR.DIST<-NULL
str(transmission.line)

paved.road<-read.csv("0_data/1_processed/human footprint/paved-road-areadist.age-allsurveys.csv", header=TRUE)
str(paved.road)
paved.road$AREA150.paved.road<-paved.road$AREA150
paved.road$AREA150<-NULL
paved.road$PROP150.paved.road<-paved.road$PROP150
paved.road$PROP150<-NULL
paved.road$MEANAGE.150.paved.road<-paved.road$MEANAGE.150
paved.road$MEANAGE.150<-NULL
paved.road$AREA565.paved.road<-paved.road$AREA565
paved.road$AREA565<-NULL
paved.road$PROP565.paved.road<-paved.road$PROP565
paved.road$PROP565<-NULL
paved.road$MEANAGE.565.paved.road<-paved.road$MEANAGE.565
paved.road$MEANAGE.565<-NULL
paved.road$NEAR.DIST.paved.road<-paved.road$NEAR.DIST
paved.road$NEAR.DIST<-NULL
str(paved.road)

gravel.road<-read.csv("0_data/1_processed/human footprint/gravel-road-areadist.age-allsurveys.csv", header=TRUE)
str(gravel.road)
gravel.road$AREA150.gravel.road<-gravel.road$AREA150
gravel.road$AREA150<-NULL
gravel.road$PROP150.gravel.road<-gravel.road$PROP150
gravel.road$PROP150<-NULL
gravel.road$MEANAGE.150.gravel.road<-gravel.road$MEANAGE.150
gravel.road$MEANAGE.150<-NULL
gravel.road$AREA565.gravel.road<-gravel.road$AREA565
gravel.road$AREA565<-NULL
gravel.road$PROP565.gravel.road<-gravel.road$PROP565
gravel.road$PROP565<-NULL
gravel.road$MEANAGE.565.gravel.road<-gravel.road$MEANAGE.565
gravel.road$MEANAGE.565<-NULL
gravel.road$NEAR.DIST.gravel.road<-gravel.road$NEAR.DIST
gravel.road$NEAR.DIST<-NULL
str(gravel.road)

unimproved.road<-read.csv("0_data/1_processed/human footprint/unimproved-road-areadist.age-allsurveys.csv", header=TRUE)
str(unimproved.road)
unimproved.road$AREA150.unimproved.road<-unimproved.road$AREA150
unimproved.road$AREA150<-NULL
unimproved.road$PROP150.unimproved.road<-unimproved.road$PROP150
unimproved.road$PROP150<-NULL
unimproved.road$MEANAGE.150.unimproved.road<-unimproved.road$MEANAGE.150
unimproved.road$MEANAGE.150<-NULL
unimproved.road$AREA565.unimproved.road<-unimproved.road$AREA565
unimproved.road$AREA565<-NULL
unimproved.road$PROP565.unimproved.road<-unimproved.road$PROP565
unimproved.road$PROP565<-NULL
unimproved.road$MEANAGE.565.unimproved.road<-unimproved.road$MEANAGE.565
unimproved.road$MEANAGE.565<-NULL
unimproved.road$NEAR.DIST.unimproved.road<-unimproved.road$NEAR.DIST
unimproved.road$NEAR.DIST<-NULL
str(unimproved.road)

truck.trail<-read.csv("0_data/1_processed/human footprint/truck-trail-areadist.age-allsurveys.csv", header=TRUE)
str(truck.trail)
truck.trail$AREA150.truck.trail<-truck.trail$AREA150
truck.trail$AREA150<-NULL
truck.trail$PROP150.truck.trail<-truck.trail$PROP150
truck.trail$PROP150<-NULL
truck.trail$MEANAGE.150.truck.trail<-truck.trail$MEANAGE.150
truck.trail$MEANAGE.150<-NULL
truck.trail$AREA565.truck.trail<-truck.trail$AREA565
truck.trail$AREA565<-NULL
truck.trail$PROP565.truck.trail<-truck.trail$PROP565
truck.trail$PROP565<-NULL
truck.trail$MEANAGE.565.truck.trail<-truck.trail$MEANAGE.565
truck.trail$MEANAGE.565<-NULL
truck.trail$NEAR.DIST.truck.trail<-truck.trail$NEAR.DIST
truck.trail$NEAR.DIST<-NULL
str(truck.trail)

active.well<-read.csv("0_data/1_processed/human footprint/active-well-areadist.age-allsurveys.csv", header=TRUE)
str(active.well)
active.well$AREA150.active.well<-active.well$AREA150
active.well$AREA150<-NULL
active.well$PROP150.active.well<-active.well$PROP150
active.well$PROP150<-NULL
active.well$MEANAGE.150.active.well<-active.well$MEANAGE.150
active.well$MEANAGE.150<-NULL
active.well$AREA565.active.well<-active.well$AREA565
active.well$AREA565<-NULL
active.well$PROP565.active.well<-active.well$PROP565
active.well$PROP565<-NULL
active.well$MEANAGE.565.active.well<-active.well$MEANAGE.565
active.well$MEANAGE.565<-NULL
active.well$NEAR.DIST.active.well<-active.well$NEAR.DIST
active.well$NEAR.DIST<-NULL
str(active.well)

abandoned.well<-read.csv("0_data/1_processed/human footprint/abandoned-well-areadist.age-allsurveys.csv", header=TRUE)
str(abandoned.well)
abandoned.well$AREA150.abandoned.well<-abandoned.well$AREA150
abandoned.well$AREA150<-NULL
abandoned.well$PROP150.abandoned.well<-abandoned.well$PROP150
abandoned.well$PROP150<-NULL
abandoned.well$MEANAGE.150.abandoned.well<-abandoned.well$MEANAGE.150
abandoned.well$MEANAGE.150<-NULL
abandoned.well$AREA565.abandoned.well<-abandoned.well$AREA565
abandoned.well$AREA565<-NULL
abandoned.well$PROP565.abandoned.well<-abandoned.well$PROP565
abandoned.well$PROP565<-NULL
abandoned.well$MEANAGE.565.abandoned.well<-abandoned.well$MEANAGE.565
abandoned.well$MEANAGE.565<-NULL
abandoned.well$NEAR.DIST.abandoned.well<-abandoned.well$NEAR.DIST
abandoned.well$NEAR.DIST<-NULL
str(abandoned.well)

industrial<-read.csv("0_data/1_processed/human footprint/industrial-areadist.age-allsurveys.csv", header=TRUE)
str(industrial)
industrial$AREA150.industrial<-industrial$AREA150
industrial$AREA150<-NULL
industrial$PROP150.industrial<-industrial$PROP150
industrial$PROP150<-NULL
industrial$MEANAGE.150.industrial<-industrial$MEANAGE.150
industrial$MEANAGE.150<-NULL
industrial$AREA565.industrial<-industrial$AREA565
industrial$AREA565<-NULL
industrial$PROP565.industrial<-industrial$PROP565
industrial$PROP565<-NULL
industrial$MEANAGE.565.industrial<-industrial$MEANAGE.565
industrial$MEANAGE.565<-NULL
industrial$NEAR.DIST.industrial<-industrial$NEAR.DIST
industrial$NEAR.DIST<-NULL
str(industrial)

mine<-read.csv("0_data/1_processed/human footprint/mine-areadist.age-allsurveys.csv", header=TRUE)
str(mine)
mine$AREA150.mine<-mine$AREA150
mine$AREA150<-NULL
mine$PROP150.mine<-mine$PROP150
mine$PROP150<-NULL
mine$MEANAGE.150.mine<-mine$MEANAGE.150
mine$MEANAGE.150<-NULL
mine$AREA565.mine<-mine$AREA565
mine$AREA565<-NULL
mine$PROP565.mine<-mine$PROP565
mine$PROP565<-NULL
mine$MEANAGE.565.mine<-mine$MEANAGE.565
mine$MEANAGE.565<-NULL
mine$NEAR.DIST.mine<-mine$NEAR.DIST
mine$NEAR.DIST<-NULL
str(mine)

residential<-read.csv("0_data/1_processed/human footprint/residential-areadist.age-allsurveys.csv", header=TRUE)
str(residential)
residential$AREA150.residential<-residential$AREA150
residential$AREA150<-NULL
residential$PROP150.residential<-residential$PROP150
residential$PROP150<-NULL
residential$MEANAGE.150.residential<-residential$MEANAGE.150
residential$MEANAGE.150<-NULL
residential$AREA565.residential<-residential$AREA565
residential$AREA565<-NULL
residential$PROP565.residential<-residential$PROP565
residential$PROP565<-NULL
residential$MEANAGE.565.residential<-residential$MEANAGE.565
residential$MEANAGE.565<-NULL
residential$NEAR.DIST.residential<-residential$NEAR.DIST
residential$NEAR.DIST<-NULL
str(residential)

cultivation<-read.csv("0_data/1_processed/human footprint/cultivation-areadist.age-allsurveys.csv", header=TRUE)
str(cultivation)
cultivation$AREA150.cultivation<-cultivation$AREA150
cultivation$AREA150<-NULL
cultivation$PROP150.cultivation<-cultivation$PROP150
cultivation$PROP150<-NULL
cultivation$MEANAGE.150.cultivation<-cultivation$MEANAGE.150
cultivation$MEANAGE.150<-NULL
cultivation$AREA565.cultivation<-cultivation$AREA565
cultivation$AREA565<-NULL
cultivation$PROP565.cultivation<-cultivation$PROP565
cultivation$PROP565<-NULL
cultivation$MEANAGE.565.cultivation<-cultivation$MEANAGE.565
cultivation$MEANAGE.565<-NULL
cultivation$NEAR.DIST.cultivation<-cultivation$NEAR.DIST
cultivation$NEAR.DIST<-NULL
str(cultivation)


harvest<-read.csv("0_data/1_processed/human footprint/harvest-areadist.age-allsurveys.csv", header=TRUE)
str(harvest)
harvest$AREA150.harvest<-harvest$AREA150
harvest$AREA150<-NULL
harvest$PROP150.harvest<-harvest$PROP150
harvest$PROP150<-NULL
harvest$MEANAGE.150.harvest<-harvest$MEANAGE.150
harvest$MEANAGE.150<-NULL
harvest$AREA565.harvest<-harvest$AREA565
harvest$AREA565<-NULL
harvest$PROP565.harvest<-harvest$PROP565
harvest$PROP565<-NULL
harvest$MEANAGE.565.harvest<-harvest$MEANAGE.565
harvest$MEANAGE.565<-NULL
harvest$NEAR.DIST.harvest<-harvest$NEAR.DIST
harvest$NEAR.DIST<-NULL
str(harvest)

m1<-merge(conventional.seismic, low.impact.seismic, by=c("PKEY","SS","YEAR"))
m2<-merge(m1, pipeline, by=c("PKEY","SS","YEAR"))
m3<-merge(m2, transmission.line, by=c("PKEY","SS","YEAR"))
m4<-merge(m3, paved.road, by=c("PKEY","SS","YEAR"))
m5<-merge(m4, gravel.road, by=c("PKEY","SS","YEAR"))
m6<-merge(m5, unimproved.road, by=c("PKEY","SS","YEAR"))
m7<-merge(m6, truck.trail, by=c("PKEY","SS","YEAR"))
m8<-merge(m7, active.well, by=c("PKEY","SS","YEAR"))
m9<-merge(m8, abandoned.well, by=c("PKEY","SS","YEAR"))
m10<-merge(m9, industrial, by=c("PKEY","SS","YEAR"))
m11<-merge(m10, mine, by=c("PKEY","SS","YEAR"))
m12<-merge(m11, residential, by=c("PKEY","SS","YEAR"))
m13<-merge(m12, cultivation, by=c("PKEY","SS","YEAR"))
m14<-merge(m13, harvest, by=c("PKEY","SS","YEAR"))
write.csv(m14, file="0_data/1_processed/human footprint/all-footprints-areadist.age-allsurveys.csv")

