## Google Earth Engine Code Share ##

### This repository provides code for performing the Google Earth Engine analysis detailed in the below manuscript, which is accepted for publication in Global Change Biology. This code is shared in the spirit of open science. ### 

Anderson, K., Fawcett, D., Cugulliere, A., Benford, S. and Jones, D. (2019) Vegetation expansion in the subnival Hindu Kush Himalaya, Global Change Biology (in press)


## Part 1: Has the spatial extent of subnival vegetation changed and, if so, at what rate and where? ##

### P140/R40-41 region ##

This code performs the analysis to quantify time-series change in two tiles of Landsat data (path 140, rows 40 and 41; see manuscript Figure 2). The region referred to as P140/R40-41 covered an area on the Nepal/Tibet border centred on Mount Everest.

First, declare the imports as follows:

```
var pathrowextent = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[87.51708984375, 26.504988828743404],
          [88.3026123046875, 29.54000879252545],
          [86.407470703125, 29.81205076752506],
          [85.660400390625, 26.770135082241445]]]),
    SRTM = ee.Image("USGS/SRTMGL1_003"),
    OverallSnowMask = ee.Image("users/dfawcett/HimalayaSnowMask"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    heightlowres = ee.FeatureCollection("users/dfawcett/heightlowres"),
    nepalborder = ee.FeatureCollection("users/dfawcett/NepalBorder"),
    LS8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    LS7SR = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    LS5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR");
  ```
