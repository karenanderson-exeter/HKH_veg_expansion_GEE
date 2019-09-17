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

Next, run the analysis:


```
//Extracting vegetated pixel fractions over Landsat path 140 rows 40-41 region

var startYear=1993
var endYear=2018

//list of years
var years = ee.List.sequence(startYear, endYear);

//adjust elevation thresholds
var elemin=4150
var elemax=4500
//var elemin=4500
//var elemax=5000
//var elemin=5000
//var elemax=5500
//var elemin=5500
//var elemax=6000

//region of interest
var roi=pathrowextent;

//create elevation and aspect mask based on SRTM 90 m resolution
var elemask =SRTM90.clip(roi).gt(elemin).and(SRTM90.clip(roi).lt(elemax));
var aspect = ee.Terrain.aspect(SRTM.clip(roi));
var aspectmask = aspect.gte(45).and(aspect.lte(315));

//clip SRTM data to region
var SRTMclip = SRTM.clip(roi);

//display bounds of path-row region
var empty = ee.Image().byte();
var outlines = empty.paint({
  featureCollection: roi,
  color: 'FFFFFF',
  width: 2
});

Map.addLayer(outlines);

//Select Landsat datasets within region, for October and November

var LS5collROI = LS5SR
.filter(ee.Filter.calendarRange(10,11,'month'))
.filterBounds(roi)
.select(['B5','B4', 'B3', 'B2','pixel_qa'], ['SWIR','NIR', 'RED','GREEN','pixel_qa']);

var LS7collROI = LS7SR 
.filter(ee.Filter.calendarRange(10,11,'month'))
.filterBounds(roi)
.select(['B5','B4', 'B3', 'B2','pixel_qa'], ['SWIR','NIR', 'RED','GREEN','pixel_qa']);

var LS8collROI = LS8SR
.filter(ee.Filter.calendarRange(10,11,'month'))
.filterBounds(roi)
.select(['B6','B5', 'B4', 'B3','pixel_qa'], ['SWIR','NIR', 'RED','GREEN','pixel_qa']);//different band designations than previous missions

//masking clouds, haze, snow and shadow

//extract quality bits information from QA band

var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    return image.select([0], [newName])
                  .bitwise_and(pattern)
                  .right_shift(start);
};

// Function to mask out undesired pixels
var maskClouds = function(image) {
  
  // Select the QA band
  var QA = image.select('pixel_qa');

  // Get the internal_cloud_algorithm_flag bit.
  var internal_cloud_algorithm_flag = getQABits(QA, 1, 1, 'internal_cloud_algorithm_flag');
  var internal_snow_algorithm_flag = getQABits(QA,4,4,'internal_snow_algorithm_flag')
    var internal_shadow_algorithm_flag = getQABits(QA, 3, 3, 'internal_shadow_algorithm_flag');
  var internal_cloud_algorithm_flag_AOT = getQABits(QA, 6, 7, 'internal_cloud_algorithm_flag_AOT');
  // Create a mask for the image.
  var mask = internal_cloud_algorithm_flag;
  var aotmask = internal_cloud_algorithm_flag_AOT.lt(2);
  var shadmask = internal_shadow_algorithm_flag.not();
  var snowmask = internal_snow_algorithm_flag.not();
  // Return an image masking out undesired areas
  return image.mask(mask.and(aotmask).and(shadmask).and(snowmask));
};

//apply QA flag masking
var LS5and7collROI=ee.ImageCollection(LS5collROI.merge(LS7collROI))
var LS5and7collROIcmasked=LS5and7collROI.map(maskClouds)
var LS8collROIcmasked=LS8collROI.map(maskClouds)

//make an empty image per year (needed so script doesn't fail for years with missing data)
function makeImgsWithDates(year){
  var newimg=ee.Image(0).set("system:time_start", ee.Date(ee.Number(year).format('%d').cat('-01-01')).millis());
  return newimg.rename('NDVI').updateMask(0);
}

var emptyimgs = years.map(makeImgsWithDates);

//function to compute NDVI
var calcNDVI = function(image){
  var ndvi= image.normalizedDifference(['NIR','RED']).rename('NDVI');
  return image.addBands(ndvi)
}

//function to compute NDVI and correct to ETM+ (for LS8) using Roy et al, (2016)
var calccorrNDVI = function(image){
  var ndvi= image.normalizedDifference(['NIR','RED']).rename('NDVI').multiply(0.9589).add(0.0029).rename('NDVI');
  return image.addBands(ndvi)
}

//compute NDVI
var LS5and7NDVI = ee.ImageCollection(LS5and7collROIcmasked
.map(calcNDVI).merge(emptyimgs));//NDVI for LS7 

//compute LS8 NDVI and apply correction
var LS8NDVI = ee.ImageCollection(LS8collROIcmasked
.map(calccorrNDVI).merge(emptyimgs));


var LSNDVI=LS5and7NDVI.merge(LS8NDVI)


//function to calculate median NDVI per year
var calculateAnnualNDVI = function(year){
  var currentNDVI=LSNDVI.select('NDVI').filter(ee.Filter.calendarRange({start:year, end:year,field:'year'}));
  var medianNDVI=currentNDVI.select('NDVI').median().set('system:time_start',ee.Date(ee.Number(year).format('%d').cat('-01-01')).millis());//set a date per image (1. Jan every year)
  var NDVImask= medianNDVI.select('NDVI').gte(0.1);
  var NDVIthresh=medianNDVI.updateMask(NDVImask).rename('NDVIthresh');
  return medianNDVI.addBands(NDVIthresh);
};

//apply function to list of years
var LSNDVIcollfinal= ee.ImageCollection(years.map(calculateAnnualNDVI));


//generate fraction of vegetated unmasked pixels
var getVegPixelCount = function(img2red){
  
  //mask areas outside the elevation band
  img2red=img2red.updateMask(SRTMclip.gt(elemin).and(SRTMclip.lt(elemax)));//.updateMask(OverallSnowMask);//);
  
  //mask areas on north-facing slopes
  img2red=img2red.updateMask(aspectmask); 
  
  //get count of total unmasked pixels and green unmasked pixels
  var totalPixelCount =img2red.select(['NDVI','NDVIthresh']).reduceRegions({
    reducer: ee.Reducer.count(),
    collection:roi,
    scale: 90
  });
 

return totalPixelCount
}

//apply calculation of counts to all images (years)
var totalList = LSNDVIcollfinal.map(getVegPixelCount).toList(26);

//Conversions to unpack collection and calculate vegetated pixel fractions (ideally optimise in a way which does not require for loop)
var resultList = ee.List([]);
for (var i=0; i<26; i++){
  var colreducer = ee.Reducer.toList();
  var totalpix = ee.Array(ee.FeatureCollection(totalList.get(i)).reduceColumns(colreducer, ['NDVI']).get('list'));
  var greenpix = ee.Array(ee.FeatureCollection(totalList.get(i)).reduceColumns(colreducer, ['NDVIthresh']).get('list'));
  var greenfrac = greenpix.divide(totalpix);
 resultList=resultList.add(greenfrac); 
}

//Chart the resulting array, export from displayed chart
var resultarray=ee.FeatureCollection(ee.Array.cat(resultList,1));
//Export.table.toDrive(resultarray) //currently not functional
var arraychart =ui.Chart.array.values(resultarray,1)
print(arraychart)
```
