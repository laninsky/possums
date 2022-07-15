## Geospatial analyses of mitogenome clades

### Converting coordinates to decimal long/lat
We received the geographic locations for some samples in New Zealand Transverse Mercator Projection. For displaying on maps with R packages, these were converted to decimal long/lat (New Zealand Geodetic Datum 2000 (version 20180701)), using https://www.geodesy.linz.govt.nz/concord/index.cgi

For some samples, a verbal description of location was given, and an approximate long/lat was used based on these descriptions. We also added a little "jitter" for these samples so that different haplotypes could be seen on the map.

Following this, the excel spreadsheet was saved as a tab-delimited file before reading into R and creating some figures using [geospatial_code.R](geospatial_code.R). The necessary columns for this code are:
- Location # Sample locality in text format
- Longitude # Longitude without jitter added
- Jitter_Long # Jitter added for samples with same Lat/Long
- Jitter_Lat # Jitter added for samples with same Lat/Long
- RAxML run 1 # The mitochondrial phylogeny clades

e.g.
```
# A tibble: 295 × 5
   Location    Longitude Jitter_Long Jitter_Lat `RAxML run 1`  
   <chr>           <dbl>       <dbl>      <dbl> <chr>          
 1 Buffer Zone      171.        171.      -45.9 Blue (NSW)     
 2 Buffer Zone      171.        171.      -45.9 Green (Takashi)
 3 Buffer Zone      171.        171.      -45.9 Green (Takashi)
 4 Buffer Zone      171.        171.      -45.9 Red (mixed)    
 5 Buffer Zone      171.        171.      -45.9 Red (mixed)    
 6 Buffer Zone      171.        171.      -45.9 Blue (NSW)     
 7 Buffer Zone      171.        171.      -45.9 Blue (NSW)     
 8 Buffer Zone      171.        171.      -45.9 Blue (NSW)     
 9 Buffer Zone      171.        171.      -45.9 Blue (NSW)     
10 Buffer Zone      171.        171.      -45.9 Red (mixed)    
# … with 285 more rows
```
