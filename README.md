## Google Earth Engine Code Share ##

### This repository provides code for performing the Google Earth Engine analysis detailed in the below manuscript, which is currently under consideration for publication in Global Change Biology. This code is shared in the spirit of open science. ### 

*Anderson, K., Fawcett, D., Cugulliere, A., Benford, S. and Jones, D. (2019) Vegetation expansion in the subnival Hindu Kush Himalaya, Global Change Biology (in review).*

#### Note, within these scripts we implement the Landsat-7 to Landsat-8 intercalibration, documented in Roy et al (2016). If you are interested in this, you can find their original paper here: ####

*Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S. and Egorov, A. (2016) Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity. Remote sensing of Environment, 185, pp.57-70. [https://www.sciencedirect.com/science/article/pii/S0034425715302455]*


## Part 1: what is the extent of the subnival zone? ##

We used the Moderate Resolution Imaging Spectrometer (MODIS) fractional snow product to answer this question at two scales: 

1. the "P140/R40-41" region, which is the area defined by two tiles of Landsat data (path 140, rows 40 and 41; see manuscript Figure 2). The region referred to as P140/R40-41 covered an area on the Nepal/Tibet border centred on Mount Everest.

2. the national extent of Nepal. 

We selected recent years (2013-2017) to generate a product describing the median snow covered area in late Summer (August and September), when the snow cover is at a minimum (determined using: http://geoapps.icimod.org/HKHSnowCover/). This output was used to represent permanent snow-covered areas for the ROI. The permanent snow-covered area was combined with the Randolph Glacier Inventory (RGI) to calculate the spatial extent of permanent snow and ice cover. This was compared to the extent of the entire subnival zone, represented by the total area above 4150 m.a.s.l using the Shuttle Radar Topography Mission (SRTM) 30 m gridded dataset as a measurement of height above mean sea level. 

### P140/R40-41 extent ###

First, declare the imports as follows:

```javascript
var SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    geometry = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[87.51708984375, 26.504988828743404],
          [88.3026123046875, 29.54000879252545],
          [86.407470703125, 29.81205076752506],
          [85.660400390625, 26.770135082241445]]]),
    modisWater = ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24"),
    subGlaciers = ee.FeatureCollection("users/dfawcett/allAsiaGlaciers_pathrowclip");
```


Next, run the script:
```javascript
//Landsat path 140 rows 40-41 snow cover: compute snow cover fractions and ice cover
//last modified 16/09/2019

var hkh = ee.FeatureCollection("ft:1Q3InAXAA3LAa_K_VLSbafnDiofJfhpbv1k8wqZMw");

//glacier shapefiles from the Randolph Glacier Inventory
var glaciersShp= subGlaciers
var glaciersShpclip = glaciersShp.geometry().intersection({right:geometry,maxError:100});

//adjust elevation thresholds
var elemin=4150
var elemax=10000

//daily MODIS Terra snow cover dataset, reference: http://nsidc.org/data/MOD10A1
var snowCover=ee.ImageCollection("MODIS/006/MOD10A1").filterBounds(hkh).filterDate('2013-01-01','2018-01-01');

//create elevation and aspect mask based on SRTM 90 m resolution
var elemask =SRTM90.clip(hkh).gt(elemin).and(SRTM90.clip(hkh).lt(elemax));

//create vector region at 1 km resolution
var fullmask = elemask//.and(aspectmask);
fullmask= fullmask.reproject({crs:'EPSG:32645',scale:1000})

var vectorregion = fullmask.updateMask(fullmask).reduceToVectors({maxPixels: 2032924972,geometry:hkh});
Map.addLayer(vectorregion)

//mask pixels in each observation which have bad data or detected as cloud / ocean
var maskBadCover = function(img){
  return img.select("NDSI_Snow_Cover")
  .updateMask(img.select("NDSI_Snow_Cover").lte(100)) //pixels with unmasked NDSI value
  .updateMask(img.select("NDSI_Snow_Cover_Basic_QA").lte(0)) //only pixels QA flagged as "best" 
}

//median snow cover per year (autumn)
var yearlymedAutumn= ee.List([]);
for (var i=0; i<5; i++){
 var year = 2013+i;

var autumnSnowCover = snowCover.filter(ee.Filter.calendarRange({start:8, end:9,field:'month'}));
var autumnSnowCoverMasked = autumnSnowCover.map(maskBadCover);
var medAutumnSnow=autumnSnowCoverMasked.median();
//mask water regions
var medAutumnSnowNoWater=medAutumnSnow.updateMask(modisWater.select('water_mask').not())
yearlymedAutumn=yearlymedAutumn.add(medAutumnSnowNoWater)
}

//median over all years
var medianAutumnSnowNoWater=ee.ImageCollection(yearlymedAutumn).median()
//snow pixels classified as greater than 30% fractional cover, only for display
var snowpixelsAutumn=medianAutumnSnowNoWater.gte(30)
snowpixelsAutumn=snowpixelsAutumn.selfMask()

//region above 4150 m intersected with roi
var subregion=vectorregion.geometry().intersection({right:geometry,maxError:100});
Map.addLayer(subregion)
//snow cover fraction over whole area
var snowfrac = medianAutumnSnowNoWater.reduceRegion({reducer:ee.Reducer.mean(),geometry:subregion,scale:1000,tileScale:4,bestEffort:true,maxPixels:1000000000})
//snow cover fraction on glacier area
var glaciersnowfrac=medianAutumnSnowNoWater.reduceRegion({reducer:ee.Reducer.mean(),geometry:glaciersShpclip,scale:1000,tileScale:4,bestEffort:true,maxPixels:1000000000})
//total area
var allarea = subregion.area({maxError:100})
//snow covered area
var snowarea = allarea.multiply(ee.Number(snowfrac.get('NDSI_Snow_Cover')).divide(100))
//glacier area
var glacierarea=glaciersShpclip.area()
//snow covered glacier area
var snowglacierarea=glacierarea.multiply(ee.Number(glaciersnowfrac.get('NDSI_Snow_Cover')).divide(100))

print(glaciersnowfrac)
Map.addLayer(subregion)
print(allarea)
print(snowfrac)
print(snowarea)
print(glacierarea)
print(snowglacierarea)

Map.addLayer(elemask.selfMask().clip(geometry),{palette:'orange'},"subnival zone")

Map.addLayer(glaciersShp,{color: '0000FF'},'glaciers')

Map.addLayer(glaciersShpclip,{color: '0000FF'},'glaciers')

Map.addLayer(snowpixelsAutumn.clip(geometry),{palette:['white'],min:0,max:1},"snow ext autumn");

var empty = ee.Image().byte();

var outlines = empty.paint({
  featureCollection: geometry,
  color: 'FFFFFF',
  width: 2
}); 

Map.addLayer(outlines);
```
#### Example output: ####

![P140-R40/41 subnival zone analysis example](Github_P140R4041_snow.JPG?raw=true "P140-R40/41 subnival zone analysis example")



### Nepal extent ###

First, declare the imports as follows:

```javascript
var SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    modisWater = ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24"),
    nepal = ee.FeatureCollection("users/dfawcett/NepalBorder"),
    nepalGlaciers = ee.FeatureCollection("users/dfawcett/Nepal_Glaciers");
```

Then, run the script:

```javascript
//Nepalese snow cover: compute snow cover fractions and ice cover
//last modified 16/09/2019

//glacier shapefiles from the Randolph Glacier Inventory
var glaciersShp= ee.FeatureCollection("users/dfawcett/Nepal_Glaciers"); //dominic to share the file?

//adjust elevation thresholds
var elemin=4150
var elemax=10000

//daily MODIS Terra snow cover dataset, reference: http://nsidc.org/data/MOD10A1
var snowCover=ee.ImageCollection("MODIS/006/MOD10A1").filterBounds(nepal.geometry()).filterDate('2013-01-01','2018-01-01');

//create elevation and aspect mask based on SRTM 90 m resolution
var elemask =SRTM90.clip(nepal.geometry()).gt(elemin).and(SRTM90.clip(nepal.geometry()).lt(elemax));

//create vector region at 1 km resolution
var fullmask = elemask//.and(aspectmask);
fullmask= fullmask.reproject({crs:'EPSG:32645',scale:1000})
var vectorregion = fullmask.updateMask(fullmask).reduceToVectors({maxPixels: 2032924972,geometry:nepal.geometry()});

//mask pixels in each observation which have bad data or detected as cloud / ocean
var maskBadCover = function(img){
  return img.select("NDSI_Snow_Cover")
  .updateMask(img.select("NDSI_Snow_Cover").lte(100)) //pixels with unmasked NDSI value
  .updateMask(img.select("NDSI_Snow_Cover_Basic_QA").lte(0)) //only pixels QA flagged as "best" 
}

//median snow cover per year (autumn)
var yearlymedAutumn= ee.List([]);
for (var i=0; i<5; i++){
 var year = 2013+i;

var autumnSnowCover = snowCover.filter(ee.Filter.calendarRange({start:8, end:9,field:'month'}));
var autumnSnowCoverMasked = autumnSnowCover.map(maskBadCover);
var medAutumnSnow=autumnSnowCoverMasked.median();
//mask water regions
var medAutumnSnowNoWater=medAutumnSnow.updateMask(modisWater.select('water_mask').not())
yearlymedAutumn=yearlymedAutumn.add(medAutumnSnowNoWater)
}

//median over all years
var medianAutumnSnowNoWater=ee.ImageCollection(yearlymedAutumn).median()
//snow pixels classified as greater than 30% fractional cover, only for display
var snowpixelsAutumn=medianAutumnSnowNoWater.gte(30)
snowpixelsAutumn=snowpixelsAutumn.selfMask()

//region above 4150 m intersected with roi
var subregion=vectorregion.geometry().intersection({right:nepal.geometry(),maxError:100});

//snow cover fraction over whole area
var snowfrac = medianAutumnSnowNoWater.reduceRegion({reducer:ee.Reducer.mean(),geometry:subregion,scale:1000,tileScale:4,bestEffort:true,maxPixels:1000000000})
//snow cover fraction on glacier area
var glaciersnowfrac=medianAutumnSnowNoWater.reduceRegion({reducer:ee.Reducer.mean(),geometry:nepalGlaciers.geometry(),scale:1000,tileScale:4,bestEffort:true,maxPixels:1000000000})
//total area
var allarea = subregion.area({maxError:100})
//snow covered area
var snowarea = allarea.multiply(ee.Number(snowfrac.get('NDSI_Snow_Cover')).divide(100))
//glacier area
var glacierarea=nepalGlaciers.geometry().area()
//snow covered glacier area
var snowglacierarea=glacierarea.multiply(ee.Number(glaciersnowfrac.get('NDSI_Snow_Cover')).divide(100))


Map.addLayer(subregion)
print(allarea)
print(snowfrac)
print(snowarea)
print(glacierarea)
print(snowglacierarea)

Map.addLayer(elemask.selfMask(),{palette:'orange'},"subnival zone")

Map.addLayer(glaciersShp,{color: '0000FF'},'glaciers')

Map.addLayer(snowpixelsAutumn.clip(nepal.geometry()),{palette:['white'],min:0,max:1},"snow ext autumn");

var empty = ee.Image().byte();

var outlines = empty.paint({
  featureCollection: nepal,
  color: 'FFFFFF',
  width: 2
}); 

Map.addLayer(outlines);
```

#### Example output: ####

![Nepal subnival zone analysis example](Github_Nepalsnowresult.JPG?raw=true "Nepal subnival zone analysis example")



## Part 2: Has the spatial extent of subnival vegetation changed and, if so, at what rate and where? ##

### P140/R40-41 region ###

This code performs the analysis to quantify time-series change in two tiles of Landsat data (path 140, rows 40 and 41; see manuscript Figure 2). The region referred to as P140/R40-41 covered an area on the Nepal/Tibet border centred on Mount Everest.

First, declare the imports as follows:

```javascript
var pathrowextent = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[87.51708984375, 26.504988828743404],
          [88.3026123046875, 29.54000879252545],
          [86.407470703125, 29.81205076752506],
          [85.660400390625, 26.770135082241445]]]),
    SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    nepalborder = ee.FeatureCollection("users/dfawcett/NepalBorder"),
    LS8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    LS7SR = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    LS5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR");
  ```

Next, run the analysis:


```javascript
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

//Chart the resulting array
var resultarray=ee.FeatureCollection(ee.Array.cat(resultList,1));
var arraychart =ui.Chart.array.values(resultarray,1)
print(arraychart)
//Export result array as table
Export.table.toDrive(ee.FeatureCollection(ee.Feature(null, { "data": resultarray })))
```

### Nepal region ###

This code performs the analysis to quantify time-series change across the areal extent of the country of Nepal. For this, a KMZ shapefile delineating the boundary of Nepal was uploaded into GEE and used to constrain the analysis. 

First, declare the imports as follows:

```javascript
var SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    nepalborder = ee.FeatureCollection("users/dfawcett/NepalBorder"),
    LS8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    LS7SR = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    LS5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR");
```

Next, run the analysis:


```javascript
//Extracting vegetated pixel fractions over Landsat path 140 rows 40-41 region
//last modified 16/09/2019

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
var roi=nepalborder.geometry();

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

//Chart the resulting array
var resultarray=ee.FeatureCollection(ee.Array.cat(resultList,1));
var arraychart =ui.Chart.array.values(resultarray,1)
print(arraychart)
//Export array as table
Export.table.toDrive(ee.FeatureCollection(ee.Feature(null, { "data": resultarray })))

```

### HKH extent analysis ###

This code performs the analysis to quantify time-series change across the entirety of the Hindu Kush Himalayan region. A different method was followed here because GEE allocates users a fixed processing capacity, so to be able to measure change over the entire HKH, a random sampling method using regions of interest (ROIs) was necessary. We defined 100 circular ROIs with a 5 km radius and randomly deployed these within each of the four height bands previously described. The total area covered by the ROIs that were used to sample the satellite data record equalled 31,416 km2 (overlap of ROIs and height bands not taken into account). Figure 3 in the manuscript shows the spatial distribution of the different height bands across the HKH sampled using the circular ROIs. 

First, declare the imports as follows:

```javascript
var SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    LS5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
    LS7SR = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    LS8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");
```

Next, perform the analysis:

```javascript
//Extracting vegetated pixel fractions over the Hindu-Kush Himalaya (HKH)
//last modified: 16/09/2019

//SR version, NDVI threshold: 0.1, clouds, shadows, aerosols and snow masked

//load Hindu-Kush Himalaya boundaries shapefile
var hkh = ee.FeatureCollection("ft:1Q3InAXAA3LAa_K_VLSbafnDiofJfhpbv1k8wqZMw"); // dominic to share the file?

//set start and end of analysis (data too sparse before 1993)
var startYear=1993
var endYear=2018

//list of years
var years = ee.List.sequence(startYear, endYear);
print(years)

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
var roi = hkh

//create elevation and aspect mask based on SRTM 90 m resolution
var elemask =SRTM90.clip(hkh).gt(elemin).and(SRTM90.clip(hkh).lt(elemax));
var aspect = ee.Terrain.aspect(SRTM.clip(hkh));
var aspectmask = aspect.gte(45).and(aspect.lte(315)); //mask north-facing slopes

//clip SRTM data to region 
var SRTMclip = SRTM.clip(hkh);

//mask for ROI selection, only elevation based
var fullmask = elemask
fullmask= fullmask.reproject({crs:'EPSG:32645',scale:1000}) //resample to 1km for vectorisation

//convert height band area to vector for ROI placement
var vectorregion = fullmask.updateMask(fullmask).reduceToVectors({maxPixels: 2032924972,geometry:hkh});

//display extent of height-band vector region on the map
Map.addLayer(vectorregion,{color: "FF0000"});

//create random points within the height band (random seed here: 2222), then create a 5 km buffer around them
//for testing and if encountering memory issues, reduce the number of random points below (currently 100)
var randpts= ee.FeatureCollection.randomPoints(vectorregion.geometry(),100,2222)
var randptsbuff =randpts.map(function(geom){
  return geom.buffer(5000)});
  
//define random buffered points as new ROI extent
roi=randptsbuff;

var roiout=randptsbuff;

//display location of random points on map
Map.addLayer(randptsbuff.geometry(),{color:"00FF00"});


//////////////////Landsat data preprocessing

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
  // Select the QA band.
  var QA = image.select('pixel_qa');

  // Get the bit flags
  var internal_cloud_algorithm_flag = getQABits(QA, 1, 1, 'internal_cloud_algorithm_flag');
  var internal_snow_algorithm_flag = getQABits(QA,4,4,'internal_snow_algorithm_flag')
  var internal_shadow_algorithm_flag = getQABits(QA, 3, 3, 'internal_shadow_algorithm_flag');
  var internal_cloud_algorithm_flag_AOT = getQABits(QA, 6, 7, 'internal_cloud_algorithm_flag_AOT');
  // Create a mask for the image
  var mask = internal_cloud_algorithm_flag;
  var aotmask = internal_cloud_algorithm_flag_AOT.lt(2);
  var shadmask = internal_shadow_algorithm_flag.not();
  var snowmask = internal_snow_algorithm_flag.not();
  // Return an image masking out undesired areas
  return image.mask(mask.and(aotmask).and(shadmask).and(snowmask));
};

//function to compute NDVI
var calcNDVI = function(image){
  var ndvi= image.normalizedDifference(['NIR','RED']).rename('NDVI');
  return image.addBands(ndvi)
}

//function to compute NDVI and correct to ETM+ (for LS8) using Roy et al, (2016)
var calccorrNDVI = function(image){
  var ndvi= image.normalizedDifference(['NIR','RED']).multiply(0.9589).add(0.0029).rename('NDVI');
  return image.addBands(ndvi)
}

//make an empty image per year (needed so script doesn't fail for years with missing data)
function makeImgsWithDates(year){
  var newimg=ee.Image(0).set("system:time_start", ee.Date(ee.Number(year).format('%d').cat('-01-01')).millis());
  return newimg.rename('NDVI').updateMask(0);
}

var LSNDVI=ee.ImageCollection(ee.Image(0))

//function to calculate median NDVI per year for desired months
var calculateAnnualNDVI = function(year){
  var currentNDVI=LSNDVI.select('NDVI').filter(ee.Filter.calendarRange({start:year, end:year,field:'year'}));
  var medianNDVI=currentNDVI.select('NDVI').median().set('system:time_start',ee.Date(ee.Number(year).format('%d').cat('-01-01')).millis());//set a date per image (1. Jan every year)
  var NDVImask= medianNDVI.select('NDVI').gte(0.1);//apply NDVI threshold
  var NDVIthresh=medianNDVI.updateMask(NDVImask).rename('NDVIthresh');
  return medianNDVI.addBands(NDVIthresh)
};

//count the number of vegetated pixels (NDVI>0.1)
var getVegPixelCount = function(img2red){
  
  //mask areas outside the elevation band
  img2red=img2red.updateMask(SRTMclip.gt(elemin).and(SRTMclip.lt(elemax)));
  
  //mask areas on north-facing slopes
  img2red=img2red.updateMask(aspectmask); 
  
  //get count of total unmasked pixels and green unmasked pixels
  var totalPixelCount =img2red.select(['NDVI','NDVIthresh']).reduceRegion({
    geometry: roiout.geometry(),
    reducer: ee.Reducer.count(),
    scale: 90
  });
 return roiout.set('NDVI',totalPixelCount.get('NDVI'),'NDVIthresh',totalPixelCount.get('NDVIthresh'),'date',img2red.get('system:time_start'))

return totalPixelCount
}

var resultArrays = ee.List([]);
  
////////////////main function, ROIs as input
var mainFunc = function(roicoll){

roiout=ee.Feature(roicoll)

//Select Landsat datasets within region, for October and November

roi=roiout.geometry()
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
.select(['B6','B5', 'B4', 'B3','pixel_qa'], ['SWIR','NIR', 'RED','GREEN','pixel_qa']);

//apply QA flag masking
var LS5and7collROI=ee.ImageCollection(LS5collROI.merge(LS7collROI))
var LS5and7collROIcmasked=ee.ImageCollection(LS5and7collROI.map(maskClouds))
var LS8collROIcmasked= ee.ImageCollection(LS8collROI.map(maskClouds))

//create an empty image per year
var emptyimgs = years.map(makeImgsWithDates);

//compute LS5 and 7 NDVI
var LS5and7NDVI = ee.ImageCollection(LS5and7collROIcmasked
.map(calcNDVI).merge(emptyimgs));//NDVI for LS7

//compute LS8 NDVI and apply correction
var LS8NDVI = ee.ImageCollection(LS8collROIcmasked
.map(calccorrNDVI).merge(emptyimgs));

LSNDVI=LS5and7NDVI.merge(LS8NDVI)

//compute yearly composites of NDVI
var LSNDVIcollfinal= ee.ImageCollection(years.map(calculateAnnualNDVI));

//generate fraction of vegetated unmasked pixels per ROI

//apply calculation of counts per ROIs to all images
var totalList = LSNDVIcollfinal.map(getVegPixelCount).toList(26);

//Conversions to unpack collection and calculate vegetated pixel fractions (ideally optimise in a way which does not require for loop)
var resultListGreenFrac = ee.List([]);
var resultListTot = ee.List([]);
var yearlistout = ee.List([]);
for (var i=0; i<26; i++){
  var colreducer = ee.Reducer.toList();
  var totalpix = ee.Number(ee.Feature(totalList.get(i)).get('NDVI'));
  var greenpix = ee.Number(ee.Feature(totalList.get(i)).get('NDVIthresh'));
  var yearout=ee.Number(ee.Feature(totalList.get(i)).get('date'))
  var greenfrac = greenpix.divide(totalpix);
  
 resultListGreenFrac=resultListGreenFrac.add(greenfrac); 
 resultListTot=resultListTot.add(totalpix)
 yearlistout=yearlistout.add(yearout)
}

roiout=roiout.set('totalpix',resultListTot)
roiout=roiout.set('greenfrac',resultListGreenFrac)
roiout=roiout.set('years',yearlistout)

return(roiout)

}
var iterpoints = randptsbuff.toList(100)
var resultCollection= iterpoints.map(mainFunc)

//export results (charts of results currently not possible due to memory limit)

//export points for geographic evaluation
Export.table.toDrive(randpts,'randpts');

//export tables of total unmasked pixels and fraction of vegetated pixels
var greenfractionvalues = ee.FeatureCollection(resultCollection).select(["greenfrac"], null, false)
var totalpixvalues = ee.FeatureCollection(resultCollection).select(["totalpix"], null, false)
Export.table.toDrive(greenfractionvalues, "HKH100_90m_".concat(elemin.toString(),"_",elemax.toString(),"_NDVI01_seed2222_Roycorr"))
Export.table.toDrive(totalpixvalues,"HKH100_90m_".concat(elemin.toString(),"_",elemax.toString(),"_totalpix_seed2222_Roycorr"))
```
#### Example output: ####

![HKH-extent vegetation fraction](Github_HKH_extent_vegfraction.JPG?raw=true "HKH-extent vegetation fraction")

## Part 3: validation ##

Here, we provide the code for linking the satellite derived vegetation product with observations of vegetation presence/absence derived from photographic data. 

First, declare the imports as follows:

```javascript
var geometry = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[87.51708984375, 26.504988828743404],
          [88.3026123046875, 29.54000879252545],
          [86.407470703125, 29.81205076752506],
          [85.660400390625, 26.770135082241445]]]),
    SRTM = ee.Image("USGS/SRTMGL1_003"),
    SRTM90 = ee.Image("CGIAR/SRTM90_V4"),
    LS8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    LS7SR = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    LS5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
    valpoints = ee.FeatureCollection("users/dfawcett/GEEvalpointsfinalext");
```

Next, perform the analysis:

```javascript
//Validation of NDVI threshold based vegetation classification using in-situ imagery over Landsat path 140 rows 40-41
//last modified 17/09/2019

//range of years
var startYear=1993
var endYear=2017
var years = ee.List.sequence(startYear, endYear);

//adjust elevation thresholds
var elemin=4150
var elemax=10000

//region of interest
var roi = geometry

//create elevation and aspect mask based on SRTM 90 m resolution
var elemask =SRTM90.clip(roi).gt(elemin).and(SRTM90.clip(roi).lt(elemax));
var aspect = ee.Terrain.aspect(SRTM.clip(roi));
var aspectmask = aspect.gte(45).and(aspect.lte(315));

//clip SRTM data to region
var SRTMclip = SRTM.clip(roi);

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
  
  // Select the QA band.
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
  // Return an image masking out cloudy areas.
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


//function to return binary vegetation mask
var calculateAnnualNDVI = function(year){
  var currentNDVI=LSNDVI.select('NDVI').filter(ee.Filter.calendarRange({start:year, end:year,field:'year'}));
  var medianNDVI=currentNDVI.select('NDVI').median().set('system:time_start',ee.Date(ee.Number(year).format('%d').cat('-01-01')).millis());//set a date per image (1. Jan every year)
  var NDVImask= medianNDVI.select('NDVI').gte(0.1).rename('NDVIthresh');
  return NDVImask
};

//apply function to list of years
var LSNDVIcollfinal= ee.ImageCollection(years.map(calculateAnnualNDVI));

//get vegetation extent for 2017
var img2017=ee.Image(LSNDVIcollfinal.toList(25).get(24))

var imgToValidate=img2017.select('NDVIthresh')
  .updateMask(SRTMclip.gt(elemin).and(SRTMclip.lt(elemax)))
  .updateMask(aspectmask)

//get vegetation extent for 1993
var img1993=ee.Image(LSNDVIcollfinal.toList(25).get(0))

var SRTMhillshade = ee.Terrain.hillshade(SRTMclip)

//display hillshade (aesthetic only)
Map.addLayer(SRTMhillshade,{},'SRTM')

//display 2017 and 1993 vegetation extents
Map.addLayer(img2017.select('NDVIthresh')
  .updateMask(SRTMclip.gt(elemin).and(SRTMclip.lt(elemax)))
  .updateMask(aspectmask).updateMask(img2017.select('NDVIthresh').gt(0))
  .clip(roi),{palette:['FF0000']},'vegetation 2017')
 Map.addLayer(img1993.select('NDVIthresh')
   .updateMask(SRTMclip.gt(elemin).and(SRTMclip.lt(elemax)))
   .updateMask(aspectmask).updateMask(img1993.select('NDVIthresh').gt(0))
   .clip(roi),{palette:['006CFF']},'vegetation 1993')

//sample values (presence/absence) at validation point locations
var valpointssampled = imgToValidate.sampleRegions({collection: valpoints, scale: 30})

//display validation points
print(valpointssampled)
Map.addLayer(valpoints,{},'validation points')
```
#### Example output: ####

![Validation points](Github_validation.JPG?raw=true "Validation points")
