# Extract-sf-R-parallel
A repository for extracting simple features (spatial analysis in R: distances to and amounts of footprint) with or without parallel processing.
I'll be attempting to extract both distance to and amount of different footprint types (within 150 m and 565 m perhaps, as in the ABMI models). I can do those in ArcGIS but I think I will try to do them in R, both for reproducibility and so that I can use the code on short notice as new data becomes available.

For distance estimation, I will need to find a way to reduce the memory used by the analysis since it involves creating a pairwise distance matrix between each point and polygon feature, then measuring the minimum distance from each polygon to each point. Since I have information about when each point count occurred and when each footprint type was built, I think I can loop through the points and for each point and footprint type: 
1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
3. estimate minimum distance among the time-and-buffer-filtered polygons
4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.

I expect that buffering to get amounts of each footprint type within 150 and 565 m of each point will use less memory if I loop through the points and recycle the buffer for each point. I'd still do step 1 above, use the year of the survey then filter to only the footprint polygons that are older than that point count survey.

So I can use simple features in R to do that, and I've created a function to do so, but it still is very slow. It takes 5 minutes to do 100 points and that's with me limiting the estimation of nearest distance to those features within 1000 m of a point, to save time and memory. 

So we come back to the reason for this GitHub repository: to try to find a way to run sf operations in parallel on a laptop with a Windows operating system. I've tried looking at a few different methods for running this analysis via parallel processing but it looks like the dplyr package's functions do not work when I try to run them on individual points located on separate threads.

The multidplyr package is supposed to enable dplyr functions to work with parallel processing, but apparently doesn't work when those functions are applied to simple features objects rather than ordinary data frames.

Someone on Stack Exchange developed a custom function they called st_par for parallel processing of sf operations but it looks like it only works for Macs (mclapply function called in it).

I've tried using parLapply as well, and alternatively the boot function (selecting rows numbered from 1:N rather than doing actual bootstrapping).

Any ideas? Worst case scenario is that I still do my area and distance estimates within ArcGIS.
