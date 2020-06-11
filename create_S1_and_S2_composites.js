///////////////////////////////
// Create the Sentinel-1 and Sentinel-2 composites used for the classification of oil palm plantations
//
// The code generates the Sentinel-1 and Sentinel-2 composites and exports the output to Google Drive
//
// Adrià Descals - a.descals@creaf.uab.cat
// CREAF - Centre de Recerca Ecològica i Aplicacions Forestals

var textSections = '' +
  "\n SECTION 1 - SETUP " +
  "\n SECTION 2 - Create Sentinel-1 Composite " +
  "\n SECTION 3 - Create Sentinel-2 Composite " +
  "\n SECTION 4 - Export Composite to Drive " 
print(textSections)

Map.setOptions('satellite')


//_______________________________________________________________________________________________________________________
// SECTION 1 - SETUP

var endDate = '2020-01-01' // last day of the composite
var spanS2 = 6 // window size (in months) for the Sentinel-2 composite
var spanS1 = 6 // window size (in months) for the Sentinel-1 composite

// define the export area
var exportArea = ee.Geometry.Point([-73.37027048247097, 3.886242442012931]).buffer(10000).bounds();



//_______________________________________________________________________________________________________________________
// SECTION 2 - S1 Composite
{

// CODE BUILT ON THE CONTRIBUTIONS IN THE GEE FORUM. LINKS TO THE GEE FORUM THREADS:

// https://groups.google.com/forum/#!msg/google-earth-engine-developers/FuccCgkDUXQ/HuRkhp2YEQAJ
// https://groups.google.com/forum/#!msg/google-earth-engine-developers/cO0o2yoGdr0/Gp02pZhoGgAJ

var S1comp = function (endDate, span, mode) {


//FUNCTIONS ===================================================================
var thin_ASC = function (image) {

    var edgeSize = 5000;
    
    // y = a * x + b; b = 0; a = -0.2247366657 
    var a = 0.2247366657;

    // band names for displacement image algorithm
    var bandNames = ['dx', 'dy'];

    var mask = ee.Image(1).clip(image.geometry().buffer(-50));
    
    var dx1 = edgeSize;
    var dy1 = a * dx1;

    var disp1 = ee.Image([dx1, dy1]).rename(bandNames);
    var disp1Mask = mask.displace(disp1);

    var dx2 = -1 * edgeSize;
    var dy2 = a * dx2;

    var disp2 = ee.Image([dx2, dy2]).rename(bandNames);
    var disp2Mask = mask.displace(disp2);

    var maskOut = disp1Mask.multiply(disp2Mask);

    return image.mask(maskOut);
};

var thin_DSC = function (image) {

    var edgeSize = 5000;
    
    // y = a * x + b; b = 0; a = -0.2247366657 
    var a = -0.2247366657;

    // band names for displacement image algorithm
    var bandNames = ['dx', 'dy'];

    var mask = ee.Image(1).clip(image.geometry()
    .buffer(-50)
    );
    
    var dx1 = edgeSize;
    var dy1 = a * dx1;

    var disp1 = ee.Image([dx1, dy1]).rename(bandNames);
    var disp1Mask = mask.displace(disp1);

    var dx2 = -1 * edgeSize;
    var dy2 = a * dx2;

    var disp2 = ee.Image([dx2, dy2]).rename(bandNames);
    var disp2Mask = mask.displace(disp2);

    var maskOut = disp1Mask.multiply(disp2Mask);

    return image.mask(maskOut);
};

// LIA CORRECTION
function calcLIA(image) {
  // We can use the gradient of the "angle" band of the S1 image to derive the S1 azimuth angle.
  var s1_inc = image.select('angle');
  var s1_azimuth = ee.Terrain.aspect(s1_inc)
                           .reduceRegion(ee.Reducer.mean(), s1_inc.get('system:footprint'), 1000)
                           .get('aspect');
  // Here we derive the terrain slope and aspect
  var srtm = ee.Image("CGIAR/SRTM90_V4")
  var srtm_slope = ee.Terrain.slope(srtm).select('slope').rename('srtm_slope');
  var srtm_aspect = ee.Terrain.aspect(srtm).select('aspect').rename('srtm_aspect');
  // And then the projection of the slope
  var slope_projected = srtm_slope.multiply(ee.Image.constant(s1_azimuth).subtract(srtm_aspect).multiply(Math.PI/180).cos()).rename('slope_projected');
  // And finally the local incidence angle
  var lia = s1_inc.subtract(ee.Image.constant(90).subtract(ee.Image.constant(90).subtract(slope_projected))).abs().rename('lia');
  return image.addBands([lia]); //srtm_slope,srtm_aspect,slope_projected
}

function toGamma0(image) {
  var gamma0 = image.select('VH').subtract(image.select('lia').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VH_LIACorr_Gamma0');
  var gamma0_VV = image.select('VV').subtract(image.select('lia').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VV_LIACorr_Gamma0');
  return image.addBands([gamma0]).addBands(gamma0_VV);
}

function sqrCosine(image) {
  var sqrCos = image.select('VH_Natural').multiply(ee.Number(30.0).multiply(Math.PI/180.0).cos().pow(2)).divide(image.select('lia').multiply(Math.PI/180.0).cos().pow(2)).rename('VH_LIACorr_SqrCosine_Natural');
  var sqrCos_VV = image.select('VV_Natural').multiply(ee.Number(30.0).multiply(Math.PI/180.0).cos().pow(2)).divide(image.select('lia').multiply(Math.PI/180.0).cos().pow(2)).rename('VV_LIACorr_SqrCosine_Natural');
  return image.addBands([sqrCos]).addBands([sqrCos_VV]);
}

// // Functions to convert from/to dB
function toNatural(image) {
  var im1 = ee.Image(10.0).pow(image.select('VH').divide(10.0)).rename('VH_Natural');
  var im2 = ee.Image(10.0).pow(image.select('VV').divide(10.0)).rename('VV_Natural');
  return image.addBands(im1).addBands(im2); //.addBands(image.select(['angle','lia']));
}

function toDB(image) {
  var im1 = ee.Image(image.select('VH_LIACorr_SqrCosine_Natural')).log10().multiply(10.0).rename('VH_LIACorr_SqrCosine_dB')
  var im2 = ee.Image(image.select('VV_LIACorr_SqrCosine_Natural')).log10().multiply(10.0).rename('VV_LIACorr_SqrCosine_dB')
  return image.addBands(im1).addBands(im2);
}

// Transform the VV and HV for a better land cover visualization
var make_s1comp = function (GRD, bandVV, bandHV) {
  
  var GRD_vv = GRD.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                .select(bandVV);

  var GRD_vh = GRD.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                .select(bandHV);
  
  var percThresh = 50 // do the median of all images
  var converter = ee.Image(-1);
  var five = ee.Image(5.5);
  var GRD_t0_vh = GRD_vh
    //.reduce(ee.Reducer.sampleStdDev())
    .reduce(ee.Reducer.percentile([percThresh]))//.median()
    .multiply(converter)
    .subtract(five);
  var GRD_t0_vv = GRD_vv
    //.reduce(ee.Reducer.sampleStdDev())
    .reduce(ee.Reducer.percentile([percThresh]))//.median()
    .multiply(converter);

  var img = GRD_t0_vh.addBands(GRD_t0_vv).rename(['t0_vh','t0_vv'])
return img;

};

// ===================================================================

/////////////////////

if (mode=='DSC'){
    
  var GRD_DSC = ee.ImageCollection('COPERNICUS/S1_GRD')
            .filter(ee.Filter.eq('instrumentMode', 'IW'))
            .filter(ee.Filter.eq('resolution', 'H'))
            .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
            .filter(ee.Filter.date(ee.Date(endDate).advance(-span, 'month'), ee.Date(endDate)))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))            
            .map(calcLIA)
            .map(toGamma0) //.map(toNatural).map(sqrCosine).map(toDB)
            .map(thin_DSC)
            //.filterBounds(geometry)
  var GRD = GRD_DSC;
            
}else if(mode == 'ASC'){          
  
  var GRD_ASC = ee.ImageCollection('COPERNICUS/S1_GRD')
            .filter(ee.Filter.eq('instrumentMode', 'IW'))
            .filter(ee.Filter.eq('resolution', 'H'))
            .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
            .filter(ee.Filter.date(ee.Date(endDate).advance(-span, 'month'), ee.Date(endDate)))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))            
            .map(calcLIA)
            .map(toGamma0) //.map(toNatural).map(sqrCosine).map(toDB)
            .map(thin_ASC)
            //.filterBounds(geometry)
  var GRD = GRD_ASC
  }

////////////////////////////////////////////
//var GRD = GRD_DSC.merge(GRD_ASC)

////////////////////////////////////////////

//var img1 = make_s1comp(GRD, 'VV', 'VH')
var img2 = make_s1comp(GRD, 'VV_LIACorr_Gamma0', 'VH_LIACorr_Gamma0')

// Map.addLayer(img1, {min:4, max: 11,"bands":["t0_vv","t0_vh","t0_vv"]}, 'sentinel 1', true);
// Map.addLayer(img2, {min:6, max: 10,"bands":["t0_vv","t0_vh","t0_vv"]}, 'sentinel 1 corr', true);

return img2.set('system_time_start',ee.Date(endDate).millis());

};

 
// call the functions for ASC and DSC
var img1 = S1comp(endDate , spanS1, 'ASC')
var img2 = S1comp(endDate , spanS1, 'DSC')


var s1 = img2.rename(['t0_vh_DSC','t0_vv_DSC']).addBands(img1.rename(['t0_vh_ASC','t0_vv_ASC']));
//var count = img1.gt(-1).multiply(10).unmask(0).add(img2.gt(-1).multiply(1).unmask(0)).rename(['count_vh','count_vv'])

// Make the average between ASC and DSC
var s1_all = s1
var s1ASC = s1_all.select('t0_vh_ASC').addBands(s1_all.select('t0_vv_ASC'))
var s1DSC = s1_all.select('t0_vh_DSC').addBands(s1_all.select('t0_vv_DSC'))
var s1ASC = s1ASC.unmask(-999);
var s1ASC = s1ASC.where(s1ASC.eq(-999),s1DSC)
var s1DSC = s1DSC.unmask(-999);
var s1DSC = s1DSC.where(s1DSC.eq(-999),s1ASC)
var s1_vv = s1ASC.select('t0_vv_ASC').add(s1DSC.select('t0_vv_DSC')).divide(2).rename("t0_vv")
var s1_vh = s1ASC.select('t0_vh_ASC').add(s1DSC.select('t0_vh_DSC')).divide(2).rename("t0_vh")
var s1 = (s1_vv.addBands(s1_vh))
var s1 = s1.updateMask(s1.select('t0_vh').neq(-999))

// scale and convert to uint8
var S1_export = s1.multiply(12).uint8()


}



//_______________________________________________________________________________________________________________________
// SECTION 3 - Create Sentinel-2 Composite

{
// call S2 collection
var col = ee.ImageCollection('COPERNICUS/S2_SR')
           .filterDate(ee.Date(endDate).advance(-spanS2, 'month'), ee.Date(endDate))
                  
var S2 = col.map(function(im){
  
  // Mask with Cloud Probability Map 
  var cloudProbMask = im.select('MSK_CLDPRB').lt(60);

  // Mask with the Scene Classification Map 
  var SCL = im.select('SCL')
  var SCL_mask = (SCL.eq(3).or(SCL.eq(1)).or(SCL.eq(8)).or(SCL.eq(9)).or(SCL.eq(10)).or(SCL.eq(11))).not()
  
  // Mask with an additional cloud shadow mask
  var shadow_mask = im.select(['B11']).lt(1200)
    .or(im.select(['B4']).lt(200))
  
  // Compute ndvi
  var ndvi = im.normalizedDifference(['B8', 'B4']).rename('NDVI');
  
  return im
          .addBands(ndvi)
          .updateMask(cloudProbMask)
          .updateMask(shadow_mask.not())
          .updateMask(SCL_mask)
          .copyProperties(im);
})


// Create a quality mosaic with NDVI
var s2_qmosaic1 = S2.qualityMosaic('NDVI')

// var imageVisParam2 = {"opacity":1,"bands":["B4","B3","B2"],"min":180,"max":1320,"gamma":1}
// Map.addLayer(s2_qmosaic1,imageVisParam2,'s2_qmosaic1',false)

// Scale and convert to uint8
var S2_export = s2_qmosaic1.select('B.*').divide(2000).multiply(255).uint8()//.addBands(S2_count.uint8())


}




//_______________________________________________________________________________________________________________________
// SECTION - Export Composite to Drive
{

// display export area
Map.addLayer(exportArea,{},'export area')
Map.centerObject(exportArea,11)

// generate tasks 
Export.image.toDrive({
  image: S1_export,
  description: 'Sentinel-1_composite_VV_VH',
  scale: 10,
  region: exportArea,
  maxPixels:527496940,
});

Export.image.toDrive({
  image: S2_export.select('B4'),
  description: 'Sentinel-2_composite_B4',
  scale: 10,
  region: exportArea,
  maxPixels:527496940,
});


}
