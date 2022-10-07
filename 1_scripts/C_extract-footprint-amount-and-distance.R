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
 
#https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r
#create a global object with same extent as human footprint layer 
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
pointcountgroup<-list()
maxdist<-1000
for (j in 888:nlevels(as.factor(hf2018_grid_10x10$ID))){
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
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/human footprint/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

