# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:21:12 2020

@author: rulrich1
"""

#Import system modules...
import sys, os, os.path, glob, re
import rasterio
import matplotlib
from rasterio.mask import mask
from rasterio.plot import show
import fiona
from fiona.crs import from_epsg
import pycrs
import geopandas as gpd
import shapely
from shapely.geometry import box, mapping
import pandas as pd
import numpy as np
import shapefile

#Read input variables as integers
Begin_year = '2001'
End_year = '2002'

#store output files
output_path = "F:/rulrich/CASA_output_global"
#read in the smallest file to use as geospatial reference
geokwargs = rasterio.open("F:/rulrich/CASA_inputs_global/Elevation/Dem8KM.tif")
#read first band only so that this can be used in raster calculations, do NOT specify the band because all
#of the input files have only one band and when you specify a single band this converts the array into a 2D array and we need a 3D array to write the output correctly
geokwargs1 = geokwargs.read()
geokwargs1.shape
geokwargs.meta
#save the metadata as kwargs so that it can applied later
kwargs = geokwargs.meta

#create a bounding box from this file to use as a shapefile for clipping
#get the bounds using rasterio
bounds = geokwargs.bounds
#convert the bounds to shapely geometry
geom = box(*bounds)
print(geom.wkt)
#create a schema with no properties
schema = {'geometry': 'Polygon', 'properties': {}}
#create a shapefile with fiona

with fiona.open('geom.shape', 'w', driver='ESRI Shapefile', crs=geokwargs.crs.to_dict(), schema=schema) as clipshapeGYE:
    clipshapeGYE.write({'geometry': mapping(geom), 'properties': {}})
    
with fiona.open('geom.shape', 'r') as shapefile:
    Mask_file = [feature['geometry'] for feature in shapefile]
    
topt = rasterio.open("F:/rulrich/CASA_inputs_global/temperature/meantemp/topt.tif", 'r+')
topt1, topt_transf = rasterio.mask.mask(dataset=topt, shapes=Mask_file, crop=True)
topt1.shape    

#input elevation - single TIF file
elev = rasterio.open("F:/rulrich/CASA_inputs_global/Elevation/Dem8Km.tif", 'r+')
elev1, elev_transf = rasterio.mask.mask(dataset=elev, shapes=Mask_file, crop=True)
elev1.shape

#input land cover - single TIF file
umd_veg = rasterio.open("F:/rulrich/CASA_inputs_global/LandCover/LandCover8Km.tif", 'r+')
lcmsk, lcmsk_transf = rasterio.mask.mask(dataset=umd_veg, shapes=Mask_file, crop=True)
lcmsk.shape

#input soil texture map - single TIF file
gz_textm = rasterio.open("F:/rulrich/CASA_inputs_global/Texture/Soil8KmTexture.tif", 'r+')
textmsk, textmsk_transf= rasterio.mask.mask(dataset=gz_textm, shapes=Mask_file, crop=True, pad=False)
textmsk.shape

#input monthly temperature averages for the years under consideration
#User should specify directory and then the matching pattern for years of interest
#files should already be ordered chronologically at this point
Trasters = glob.glob("F:/rulrich/CASA_inputs_global/temperature/meantemp/tiff/meantemp_2001*.tif") + glob.glob("F:/rulrich/CASA_inputs_global/temperature/meantemp/tiff/meantemp_2002*.tif")
#Sort monthly average temperature raster arrays into month order for use in later calculations
T_arr = []
T_arr = [None]*(len(Trasters))
for t in Trasters:
   if "jan" in t:
      if Begin_year in t:
         T_arr[0:1] = [t]*1
      else:
         T_arr[12:13] = [t]*1
   elif "feb" in t:
      if Begin_year in t:
         T_arr[1:2] = [t]*1
      else:
         T_arr[13:14] = [t]*1
   elif "mar" in t:
      if Begin_year in t:
         T_arr[2:3] = [t]*1
      else:
         T_arr[14:15] = [t]*1
   elif "apr" in t:
      if Begin_year in t:
         T_arr[3:4] = [t]*1
      else:
         T_arr[15:16] = [t]*1
   elif "may" in t:
      if Begin_year in t:
         T_arr[4:5] = [t]*1
      else:
         T_arr[16:17] = [t]*1
   elif "jun" in t:
      if Begin_year in t:
         T_arr[5:6] = [t]*1
      else:
         T_arr[17:18] = [t]*1
   elif "jul" in t:
      if Begin_year in t:
         T_arr[6:7] = [t]*1
      else:
         T_arr[18:19] = [t]*1
   elif "aug" in t:
      if Begin_year in t:
         T_arr[7:8] = [t]*1
      else:
         T_arr[19:20] = [t]*1
   elif "sep" in t:
      if Begin_year in t:
         T_arr[8:9] = [t]*1
      else:
         T_arr[20:21] = [t]*1
   elif "oct" in t:
      if Begin_year in t:
         T_arr[9:10] = [t]*1
      else:
         T_arr[21:22] = [t]*1
   elif "nov" in t:
      if Begin_year in t:
         T_arr[10:11] = [t]*1
      else:
         T_arr[22:23] = [t]*1
   elif "dec" in t:
      if Begin_year in t:
         T_arr[11:12] = [t]*1
      else:
         T_arr[23:24] = [t]*1

tempmskAr = []
tempmskAr_transf = []
for t in T_arr:
    file = rasterio.open(t)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    tempmskAr.append(out_file)
    tempmskAr_transf.append(out_transf)
print(tempmskAr)
type(tempmskAr[0])
tempmskAr[0].shape

Tmaxrasters = glob.glob("F:/rulrich/CASA_inputs_global/temperature/maxtemp/tiff/maxtemp_2001*.tif") + glob.glob("F:/rulrich/CASA_inputs_global/temperature/maxtemp/tiff/maxtemp_2002*.tif")
Tmx_arr = []
Tmx_arr = [None]*(len(Tmaxrasters))
for t in Tmaxrasters:
   if "jan" in t:
      if Begin_year in t:
         Tmx_arr[0:1] = [t]*1
      else:
         Tmx_arr[12:13] = [t]*1
   elif "feb" in t:
      if Begin_year in t:
         Tmx_arr[1:2] = [t]*1
      else:
         Tmx_arr[13:14] = [t]*1
   elif "mar" in t:
      if Begin_year in t:
         Tmx_arr[2:3] = [t]*1
      else:
         Tmx_arr[14:15] = [t]*1
   elif "apr" in t:
      if Begin_year in t:
         Tmx_arr[3:4] = [t]*1
      else:
         Tmx_arr[15:16] = [t]*1
   elif "may" in t:
      if Begin_year in t:
         Tmx_arr[4:5] = [t]*1
      else:
         Tmx_arr[16:17] = [t]*1
   elif "jun" in t:
      if Begin_year in t:
         Tmx_arr[5:6] = [t]*1
      else:
         Tmx_arr[17:18] = [t]*1
   elif "jul" in t:
      if Begin_year in t:
         Tmx_arr[6:7] = [t]*1
      else:
         Tmx_arr[18:19] = [t]*1
   elif "aug" in t:
      if Begin_year in t:
         Tmx_arr[7:8] = [t]*1
      else:
         Tmx_arr[19:20] = [t]*1
   elif "sep" in t:
      if Begin_year in t:
         Tmx_arr[8:9] = [t]*1
      else:
         Tmx_arr[20:21] = [t]*1
   elif "oct" in t:
      if Begin_year in t:
         Tmx_arr[9:10] = [t]*1
      else:
         Tmx_arr[21:22] = [t]*1
   elif "nov" in t:
      if Begin_year in t:
         Tmx_arr[10:11] = [t]*1
      else:
         Tmx_arr[22:23] = [t]*1
   elif "dec" in t:
      if Begin_year in t:
         Tmx_arr[11:12] = [t]*1
      else:
         Tmx_arr[23:24] = [t]*1

tmaxmskAr = []
tmaxmskAr_transf = []
for t in Tmaxrasters:
    file = rasterio.open(t)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    tmaxmskAr.append(out_file)
    tmaxmskAr_transf.append(out_transf)
print(tmaxmskAr)
type(tmaxmskAr[0])
tmaxmskAr[0].shape

Tminrasters = glob.glob("F:/rulrich/CASA_inputs_global/temperature/mintemp/tiff/mintemp_2001*.tif") + glob.glob("F:/rulrich/CASA_inputs_global/temperature/mintemp/tiff/mintemp_2002*.tif")
Tmn_arr = []
Tmn_arr = [None]*(len(Tminrasters))
for t in Tminrasters:
   if "jan" in t:
      if Begin_year in t:
         Tmn_arr[0:1] = [t]*1
      else:
         Tmn_arr[12:13] = [t]*1
   elif "feb" in t:
      if Begin_year in t:
         Tmn_arr[1:2] = [t]*1
      else:
         Tmn_arr[13:14] = [t]*1
   elif "mar" in t:
      if Begin_year in t:
         Tmn_arr[2:3] = [t]*1
      else:
         Tmn_arr[14:15] = [t]*1
   elif "apr" in t:
      if Begin_year in t:
         Tmn_arr[3:4] = [t]*1
      else:
         Tmn_arr[15:16] = [t]*1
   elif "may" in t:
      if Begin_year in t:
         Tmn_arr[4:5] = [t]*1
      else:
         Tmn_arr[16:17] = [t]*1
   elif "jun" in t:
      if Begin_year in t:
         Tmn_arr[5:6] = [t]*1
      else:
         Tmn_arr[17:18] = [t]*1
   elif "jul" in t:
      if Begin_year in t:
         Tmn_arr[6:7] = [t]*1
      else:
         Tmn_arr[18:19] = [t]*1
   elif "aug" in t:
      if Begin_year in t:
         Tmn_arr[7:8] = [t]*1
      else:
         Tmn_arr[19:20] = [t]*1
   elif "sep" in t:
      if Begin_year in t:
         Tmn_arr[8:9] = [t]*1
      else:
         Tmn_arr[20:21] = [t]*1
   elif "oct" in t:
      if Begin_year in t:
         Tmn_arr[9:10] = [t]*1
      else:
         Tmn_arr[21:22] = [t]*1
   elif "nov" in t:
      if Begin_year in t:
         Tmn_arr[10:11] = [t]*1
      else:
         Tmn_arr[22:23] = [t]*1
   elif "dec" in t:
      if Begin_year in t:
         Tmn_arr[11:12] = [t]*1
      else:
         Tmn_arr[23:24] = [t]*1

tminmskAr = []
tminmskAr_transf = []
for t in Tminrasters:
    file = rasterio.open(t)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    tminmskAr.append(out_file)
    tminmskAr_transf.append(out_transf)
print(tminmskAr)
type(tminmskAr[0])
tminmskAr[0].shape

Prasters = glob.glob("F:/rulrich/CASA_inputs_global/PPT_tiff/ppt_2001*.tif") + glob.glob("F:/rulrich/CASA_inputs_global/PPT_tiff/ppt_2002*.tif")
#sort files according to digits in the filename
P_arr = []
P_arr = [None]*(len(Prasters))
for p in Prasters:
   if "jan" in p:
      if Begin_year in p:
         P_arr[0:1] = [p]*1
      else:
         P_arr[12:13] = [p]*1
   elif "feb" in p:
      if Begin_year in p:
         P_arr[1:2] = [p]*1
      else:
         P_arr[13:14] = [p]*1
   elif "mar" in p:
      if Begin_year in p:
         P_arr[2:3] = [p]*1
      else:
         P_arr[14:15] = [p]*1
   elif "apr" in p:
      if Begin_year in p:
         P_arr[3:4] = [p]*1
      else:
         P_arr[15:16] = [p]*1
   elif "may" in p:
      if Begin_year in p:
         P_arr[4:5] = [p]*1
      else:
         P_arr[16:17] = [p]*1
   elif "jun" in p:
      if Begin_year in p:
         P_arr[5:6] = [p]*1
      else:
         P_arr[17:18] = [p]*1
   elif "jul" in p:
      if Begin_year in p:
         P_arr[6:7] = [p]*1
      else:
         P_arr[18:19] = [p]*1
   elif "aug" in p:
      if Begin_year in p:
         P_arr[7:8] = [p]*1
      else:
         P_arr[19:20] = [p]*1
   elif "sep" in p:
      if Begin_year in p:
         P_arr[8:9] = [p]*1
      else:
         P_arr[20:21] = [p]*1
   elif "oct" in p:
      if Begin_year in p:
         P_arr[9:10] = [p]*1
      else:
         P_arr[21:22] = [p]*1
   elif "nov" in p:
      if Begin_year in p:
         P_arr[10:11] = [p]*1
      else:
         P_arr[22:23] = [p]*1
   elif "dec" in p:
      if Begin_year in p:
         P_arr[11:12] = [p]*1
      else:
         P_arr[23:24] = [p]*1

pptmskAr = []
pptmskAr_transf = []
for p in P_arr:
    file = rasterio.open(p)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    pptmskAr.append(out_file)
    pptmskAr_transf.append(out_transf)
print(pptmskAr)
type(pptmskAr[0])
pptmskAr[0].shape

    
Srasters = glob.glob("F:/rulrich/CASA_inputs_global/solar/tiff/solar_2001*.tif") + glob.glob("F:/rulrich/CASA_inputs_global/solar/tiff/solar_2002*.tif")
S_arr = []
S_arr = [None]*(len(Srasters))
for s in Srasters:
   if "jan" in s:
      if Begin_year in s:
         S_arr[0:1] = [s]*1
      else:
         S_arr[12:13] = [s]*1
   elif "feb" in s:
      if Begin_year in s:
         S_arr[1:2] = [s]*1
      else:
         S_arr[13:14] = [s]*1
   elif "mar" in s:
      if Begin_year in s:
         S_arr[2:3] = [s]*1
      else:
         S_arr[14:15] = [s]*1
   elif "apr" in s:
      if Begin_year in s:
         S_arr[3:4] = [s]*1
      else:
         S_arr[15:16] = [s]*1
   elif "may" in s:
      if Begin_year in s:
         S_arr[4:5] = [s]*1
      else:
         S_arr[16:17] = [s]*1
   elif "jun" in s:
      if Begin_year in s:
         S_arr[5:6] = [s]*1
      else:
         S_arr[17:18] = [s]*1
   elif "jul" in s:
      if Begin_year in s:
         S_arr[6:7] = [s]*1
      else:
         S_arr[18:19] = [s]*1
   elif "aug" in s:
      if Begin_year in s:
         S_arr[7:8] = [s]*1
      else:
         S_arr[19:20] = [s]*1
   elif "sep" in s:
      if Begin_year in s:
         S_arr[8:9] = [s]*1
      else:
         S_arr[20:21] = [s]*1
   elif "oct" in s:
      if Begin_year in s:
         S_arr[9:10] = [s]*1
      else:
         S_arr[21:22] = [s]*1
   elif "nov" in s:
      if Begin_year in s:
         S_arr[10:11] = [s]*1
      else:
         S_arr[22:23] = [s]*1
   elif "dec" in s:
      if Begin_year in s:
         S_arr[11:12] = [s]*1
      else:
         S_arr[23:24] = [s]*1

solarmskAr = []
solarmskAr_transf = []
for s in Srasters:
    file = rasterio.open(s)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    solarmskAr.append(out_file)
    solarmskAr_transf.append(out_transf)
print(solarmskAr)
type(solarmskAr[0])
solarmskAr[0].shape


#EVI, soil carbon and long-term monthly temperature averages only have 12 inputs 
#Soil carbon and long-term monthly averages only used when calculating nitrogen, which this version does not
Vrasters = glob.glob("C:/CASA-Express/EVI_MODIS/evi_2002*.tif")
V_arr = []
V_arr = [None]*(len(Vrasters))
for v in Vrasters:
    if "jan" in v:
        V_arr[0:1] = [v]*1
    elif "feb" in v:
        V_arr[1:2] = [v]*1
    elif "mar" in v:
        V_arr[2:3] = [v]*1
    elif "apr" in v:
        V_arr[3:4] = [v]*1
    elif "may" in v:
        V_arr[4:5] = [v]*1
    elif "jun" in v:
        V_arr[5:6] = [v]*1
    elif "jul" in v:
        V_arr[6:7] = [v]*1
    elif "aug" in v:
        V_arr[7:8] = [v]*1
    elif "sep" in v:
        V_arr[8:9] = [v]*1
    elif "oct" in v:
        V_arr[9:10] = [v]*1
    elif "nov" in v:
        V_arr[10:11] = [v]*1
    elif "dec" in v:
        V_arr[11:12] = [v]*1
vegmskAr = []
vegmskAr_transf = []
for v in Vrasters:
    file = rasterio.open(v)
    print(file.bounds)
    print(type(file))
    out_file, out_transf = rasterio.mask.mask(dataset=file, shapes=Mask_file, all_touched=False, crop=True)    
    vegmskAr.append(out_file)
    vegmskAr_transf.append(out_transf)
print(vegmskAr)
type(vegmskAr[0])
vegmskAr[0].shape

#Process: T1 Calculation...
#here a warning will come up for T1_square that stops the code when it is issued
#it is an overflow warning, meaning that the squaring operation produces a number that is too large to be 
#stored in a float32 dtype. We must manually specify the dtype ourselves to ensure that overflow is not encountered
T1_square = np.square(topt1, dtype='float64')
T1 = np.where(topt1 < 0, 0, 0.8 + 0.02 * topt1 - 0.0005 * T1_square)
with rasterio.open(output_path + "\\T1.tif", mode='w', **kwargs) as file:
    file.write(T1.astype(rasterio.float32))

#constant values for soil moisture calculations
ETIM = 0.0
lagpfc = 60.0
lambdavalue = 1.75
lamSOM = 0.25
densitySOM = 1300.0
densityssc = 2650.0
#MODIS EVI constant
emax = 0.55

#create rasters outside the loop that are used without change within the loop
#none are saved to a file only held as variables in the coding environment
#6-tiered nested conditional statement for FSC
FSC = np.where(textmsk == 1, 0.478, np.where(textmsk == 2, 0.44, np.where(textmsk == 3, 0.478, np.where(textmsk == 4, 0.48, np.where(textmsk == 5, 0.543, np.where(textmsk == 6, 0.555, 0.44))))))

FWP = np.where(textmsk == 1, 0.13, np.where(textmsk == 2, 0.08, np.where(textmsk == 3, 0.13, np.where(textmsk == 4, 0.17, np.where(textmsk == 5, 0.27, np.where(textmsk == 6, 0.39, 0.13))))))

FFC = np.where(textmsk == 1, 0.23, np.where(textmsk == 2, 0.2, np.where(textmsk == 3, 0.23, np.where(textmsk == 4, 0.30, np.where(textmsk == 5, 0.40, np.where(textmsk == 6, 0.52, 0.27))))))

m_SAX_A = np.where((textmsk == 1) | (textmsk == 2) | (textmsk == 3), 0.002, np.where(textmsk == 4, 0.013, np.where(textmsk == 5, 0.006, np.where(textmsk == 6, 0.001, 0.002))))

m_SAX_B = np.where((textmsk == 1) | (textmsk == 3), -6.54, np.where(textmsk == 6, -13.78, np.where(textmsk == 4, -6.57, np.where(textmsk == 5, -9.47, -5.48))))

LAImax = np.where(lcmsk == 1, 7, np.where((lcmsk == 2) | (lcmsk == 3), 8, np.where((lcmsk == 4) | (lcmsk == 6) | (lcmsk == 7) | (lcmsk == 8) | (lcmsk == 10) | (lcmsk == 13), 5, np.where(lcmsk == 5, 7.5, 0))))

h1 = np.where(lcmsk == 1, 1, np.where((lcmsk == 2) | (lcmsk == 10) | (lcmsk == 13), 3, np.where((lcmsk == 3) | (lcmsk == 4) | (lcmsk == 5) | (lcmsk == 11), 10, np.where((lcmsk == 6) | (lcmsk == 7), 5, 0))))

h2 = np.where(textmsk < 7, 20, 10)

h3 = np.where(textmsk == 7, 0.001, np.where((lcmsk == 4) | (lcmsk == 7) | (lcmsk == 10), 80, np.where(lcmsk == 14, 75, 180)))

#Initialize Layers for soil moisture should be done outside the loop so that it is done one time only
#these layers will be updated at the very end of the code so that the prior month can inform the following month
#function np.zeros_like creates a new array of zeros that is the same shape and data type as the specified array, here 'topt'
PACKt0 = np.zeros_like(topt1)
X1int0 = np.zeros_like(topt1)
X2int0 = np.zeros_like(topt1)
X3int0 = np.zeros_like(topt1) 
X1t0 = np.full_like(topt1, 0.1)
X2t0 = np.full_like(topt1, 0.1)
X3t0 = np.full_like(topt1, 0.1)
RDR2t0 = np.full_like(topt1, 0.1)
RDR3t0 = np.full_like(topt1, 0.1)
TDDavt0 = np.ones_like(topt1)
FDDavt0 = np.ones_like(topt1)
SOILM0t0 = np.zeros_like(topt1)
#at the time of writing these two variables are not updated because they create SOILM0, a flooding variable that is not considered in this model
Minput_0 = np.zeros_like(topt1)
Moutput_0 = np.zeros_like(topt1)
#these are initialized based upon unchanging rasters used in code but will be updated at end of code
SOILM1t0 = FWP * h1 
SOILM2t0 = FWP * h2 
SOILM3t0 = FWP * h3

for i in range(len(tempmskAr)):
    #globally ignore any division by zero - check with Chris
    #globally ignore 'overflow' errors for floating point overflow that is too large to be expressed, too many decimal points for the type of data
    #ignore other invalid values
    with np.errstate (invalid='ignore', divide='ignore', over = 'ignore', under = 'ignore'):
    #set extent and snap raster to ensure that mask file extent perseveres throughout code
    #arcpy.env.extent = tempmskAr[i]
    #arcpy.env.snapRaster = tempmskAr[i]
        if i > 11:     
            j = i - 11
            #X1, X2, X3, SOILM0 calculations are second-tier (orange) variables and will be updated at end of code
            #all masked files are also second-tier variables
            X1 = X1t0 + X1int0 - X1t0
            X2 = X2t0 + X2int0 - X2t0
            X3 = X3t0 + X3int0 - X3t0
            SOILM0 = SOILM0t0 + Minput_0 - Moutput_0
            with rasterio.open(output_path + "\\sm0_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(SOILM0.astype(rasterio.float32))
            
            #third-tier (turquoise variables) begin here: SVPmax, SVPmin, SNOW, days, monthlyLAI, Taiga_BD, Taiga_H, Tundra_BD, Tundra_H
            #create constants to calculate SVPmax
            power1_SVPmax = (0.00738 * tmaxmskAr[i] + 0.8072) ** 8.0
            absolute1_SVPmax = np.abs(1.8 * tmaxmskAr[i] + 48.0)
            SVPmax = 33.8639 * (power1_SVPmax - 0.000019 * absolute1_SVPmax + 0.001316)
            #create constants to calculate SVPmin
            power1_SVPmin = (0.00738 * tminmskAr[i] + 0.8072) ** 8.0
            absolute1_SVPmin = np.abs(1.8 * tminmskAr[i] + 48.0)
            SVPmin = 33.8639 * (power1_SVPmin - 0.000019 * absolute1_SVPmin + 0.001316)
            #SNOW calculation
            SNOW = np.where(tempmskAr[i] < 1, pptmskAr[i], 0)
            with rasterio.open(output_path + "\\snow_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(SNOW.astype(rasterio.float32))
            #Taiga_BD, Taiga_H, Tundra_BD, Tundra_H
            #forest land cover classes return a 1 indicating taiga and 0 indicating tundra
            forest = np.where((lcmsk == 11) | (lcmsk == 1) | (lcmsk == 2) | (lcmsk == 3) |  (lcmsk == 4) | (lcmsk == 5), 1, 0)
            #all non-forest classes for taiga output have -9999 values and all forest classes for tundra output have -9999 values both indicating 'nodata'
            Taiga_BD = np.select(condlist = [(forest == 1) & (j == 1), 
                                             (forest == 1) & (j == 2), 
                                             (forest == 1) & (j == 3), 
                                             (forest == 1) & (j == 4), 
                                             (forest == 1) & (j==5), 
                                             (forest == 1) & (j>=6) & (j<=9),
                                             (forest == 1) & (j==10),
                                             (forest == 1) & (j==11),
                                             (forest == 1) & (j==12)],
                                 choicelist = [181.7, 198.6, 215.4, 232.8, 247.9, 250.0, 231.7, 146.7, 164.1],
                                 default = -9999)
            Taiga_H = np.select(condlist = [(forest == 1) & (j == 1), 
                                             (forest == 1) & (j == 2), 
                                             (forest == 1) & (j == 3), 
                                             (forest == 1) & (j == 4), 
                                             (forest == 1) & (j==5), 
                                             (forest == 1) & (j>=6) & (j<=9),
                                             (forest == 1) & (j==10),
                                             (forest == 1) & (j==11),
                                             (forest == 1) & (j==12)],
                                 choicelist = [8.2, 0.3, 0.6, 1.2, 1.9, 2.0, 1.7, 0.0, 0.1],
                                 default = -9999)
            Tundra_BD = np.select(condlist = [(forest == 0) & (j == 1), 
                                             (forest == 0) & (j == 2), 
                                             (forest == 0) & (j == 3), 
                                             (forest == 0) & (j == 4), 
                                             (forest == 0) & (j==5), 
                                             (forest == 0) & (j>=6) & (j<=9),
                                             (forest == 0) & (j==10),
                                             (forest == 0) & (j==11),
                                             (forest == 0) & (j==12)],
                                 choicelist = [7.7, 9.3, 11.1, 13.3, 15.4, 15.7, 13.9, 5.2, 6.3],
                                 default = -9999)
            Tundra_H = np.select(condlist = [(forest == 0) & (j == 1), 
                                             (forest == 0) & (j == 2), 
                                             (forest == 0) & (j == 3), 
                                             (forest == 0) & (j == 4), 
                                             (forest == 0) & (j==5), 
                                             (forest == 0) & (j>=6) & (j<=9),
                                             (forest == 0) & (j==10),
                                             (forest == 0) & (j==11),
                                             (forest == 0) & (j==12)],
                                 choicelist = [0.2, 0.3, 0.6, 1.2, 1.9, 2.0, 1.7, 0.0, 0.1],
                                 default = -9999)
                
            with rasterio.open(output_path + "\\taiga_bd_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Taiga_BD.astype(rasterio.float32))
            with rasterio.open(output_path + "\\taiga_h_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Taiga_H.astype(rasterio.float32))
            with rasterio.open(output_path + "\\tundra_bd_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Tundra_BD.astype(rasterio.float32))
            with rasterio.open(output_path + "\\tundra_h_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Tundra_H.astype(rasterio.float32))
            #set days constant for second year months - if days only enters solar equation then do NOT need to do for second year months
            if j == 4 or j == 6 or j == 9 or j == 11:
                days = 30
            elif j == 2:
                if (End_year == 2000 or End_year == 2004 or End_year == 2008 or End_year == 2012 or End_year == 2016 or End_year == 2020):
                    days = 29
                else:
                    days = 28
            else:
                days = 31
            #set constant LAI values for all of second year months
            if j == 1 or j == 2 or j == 3 or j == 11 or j == 12:
                monthlyLAI = 0
            elif j == 4:
                monthlyLAI = 0.23
            elif j == 5:
                monthlyLAI = 1.68
            elif j == 6:
                monthlyLAI = 3.85
            elif j == 7:
                monthlyLAI = 4
            elif j == 8:
                monthlyLAI = 2.28
            elif j == 9:
                monthlyLAI = 1.08
            else:
                monthlyLAI = 0.16
            
        #fourth tier variables (black) begin here: T2, Tffmo, solar, PT_a, PT_b, MELT
        #process T2: values between 0 and 1 to correspond to least and most optimum temperatures only done for second year months
            exp1 = np.exp(0.2 * (topt1 - 10 - tempmskAr[i]))
            exp2 = np.exp(0.3 * (tempmskAr[i] - 10 - topt1))
            T2 = np.where((1.1814 / (1 + exp1) / (1 + exp2)) < 0, 0, (1.1814 / (1 + exp1) / (1 + exp2)))
            with rasterio.open(output_path + "\\T2mo_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(T2.astype(rasterio.float32))
                #Tffmo calculation - not saved to a file, create constants that feed Tffmo
            power1_Tffmo = tempmskAr[i] ** 3.0
            #this should yield a constant
            naturallogbasee1_Tffmo = np.log(1 + monthlyLAI)
            #this should yield an array
            naturallogbasee2_Tffmo = np.log(1 + LAImax)
            #had to change this here because the original first argument to test conditions on did not make sense - test with Chris
            Tffmo = np.where(LAImax > monthlyLAI,
                             (tempmskAr[i] + ((-0.11 + 0.96 * tempmskAr[i] - 0.00008 * power1_Tffmo) - tempmskAr[i]) * naturallogbasee1_Tffmo / naturallogbasee2_Tffmo),
                             (tempmskAr[i] + ((-0.11 + 0.96 * tempmskAr[i] - 0.00008 * power1_Tffmo) - tempmskAr[i]) * naturallogbasee2_Tffmo / naturallogbasee2_Tffmo))
            #Process: Solar Calculation only done for second year months
            S = solarmskAr[i] * 0.5 * 0.0864 * float(days)
            with rasterio.open(output_path + "\\Smo_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(S.astype(rasterio.float32))
            #create constants to evaluate PT_a
            power1_PT = (38 - (2*(elev1 / 305)) + 380 / (SVPmax - SVPmin)) ** -1.0
            #PT_a and PT_b calculation
            PT_a = np.where(power1_PT > 0.0001, power1_PT, 0.0001)
            with rasterio.open(output_path + "\\PTa" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(PT_a.astype(rasterio.float32))
            PT_b = 2.5 + 0.14 * (SVPmax - SVPmin) + elev1 / 550
            #MELT calculation
            MELT = np.where((tempmskAr[i] > 0) & (X1t0 > 0), PACKt0, 0)
            with rasterio.open(output_path + "\\Melt_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(MELT.astype(rasterio.float32))
        
            #5th-tier variables begin here (melon-colored): RES1, RES2, RES3, TDDavIN, FDDavIN, PET, PACK
            #RES1 calculation - recall lamSOM is a constant value
            RES1 = np.where(tempmskAr[i] > 0, h1 / 100 / lamSOM, h1 / 100 / 0.5)
            #RES2 calculation - recall lambdavalue is a constant value
            RES2 = (h2 / 100) / lambdavalue
            #RES3 calculation - recall lambdavalue is a constant value
            RES3 = (h3 / 100) / lambdavalue
            #TDDavIN calculation
            TDDavIN = np.where(ETIM == 0, 
                               np.where(Tffmo > 0, Tffmo * float(days), 0), 
                                        np.where(Tffmo > 0, Tffmo, 0))
            #FDDavIN calculation
            FDDavIN = np.where(ETIM == 0, 
                               np.where(Tffmo <= 0, Tffmo *(-1) * float(days), 0), 
                                        np.where(Tffmo <= 0, Tffmo * (-1), 0))
            #Process: Potential Evapotranspiration Calculation, PET-Priestley-Taylor equations uses Tmx and Tmn
            PET = np.where(((PT_a * (tempmskAr[i] + PT_b) * solarmskAr[i] * 36.16704) / (-9.9255 * tempmskAr[i] + 10468.9636) * float(days)) > 0, 
                           (PT_a * (tempmskAr[i] + PT_b) * solarmskAr[i] * 36.16704)/(-9.9255 * tempmskAr[i] + 10468.9636) * float(days), 0)
            with rasterio.open(output_path + "\\pet_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(PET.astype(rasterio.float32))
            #PACK calculation
            PACK = PACKt0 + SNOW - MELT
            with rasterio.open(output_path + "\\pack_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(PACK.astype(rasterio.float32))
        
            #sixth-tier variables (brown) begin here: PPTint, Tundra_D, Taiga_D
            #PPTint calculation - not saved to a file
            PPTint = np.where(PET < pptmskAr[i], 0.05 * monthlyLAI * PET, 0.05 * monthlyLAI * pptmskAr[i])
            with rasterio.open(output_path + "\\PPTint_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(PPTint.astype(rasterio.float32))
            #keep the no data value -9999 as -9999 for the depth output
            Tundra_D = np.where(Tundra_BD == -9999, -9999, PACK / (Tundra_BD/(10 ** 3)))
            with rasterio.open(output_path + "\\tundra_d_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Tundra_D.astype(rasterio.float32))
            #keep the no data value -9999 as -9999 for the depth output
            Taiga_D = np.where(Taiga_BD == -9999, -9999, PACK / (Taiga_BD/(10 ** 3)))
            with rasterio.open(output_path + "\\taiga_d_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(Taiga_D.astype(rasterio.float32))
                
            #7th-tier variables (purple) begin here: Mnet
            #Mnet calculation
            Mnet = np.where((tempmskAr[i] > 0) & (X1t0 > 0), pptmskAr[i] - PPTint + MELT - PET, 0)
        
            #8th-tier variables (green) begin here: Mo1, Min1
            #Min1 calculation
            Min1 = np.where((0.2 * h1 - SOILM1t0) < (Mnet + MELT), 0.2 * h1 - SOILM1t0, Mnet + MELT)
            #Mo1 calculation
            absolute_Mo1 = np.abs(Mnet)
            Mo1 = np.where((SOILM1t0 * (X1 / h1) - (FWP * (X1 / h1) * h1)) < absolute_Mo1, SOILM1t0 * (X1 / h1) - (FWP * (X1 / h1) * h1), absolute_Mo1)
        
            #9th-tier variables (pink) begin here: Moutput_1, Minput_1
            #Minput_1
            Minput_1 = np.where(tempmskAr[i] > 0, 
                                np.where(Mnet > 0,
                                         np.where(Min1 < 0, 0, Min1), 0), 0)
            #Moutput_1
            Moutput_1 = np.where(SOILM0 == 0, Mo1, 0)
    
            #10th-tier variables (darker blue-green): SOILM1, Mo2, Min2
            #SOILM1 - saved to file
            SOILM1 = SOILM1t0 + (Minput_1 - Moutput_1)
            with rasterio.open(output_path + "\\sm1_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(SOILM1.astype(rasterio.float32))
            #Min2 
            Min2 = np.where(((FFC * h2) - SOILM2t0) < ((Mnet + MELT) - Minput_1), FFC * h2 - SOILM2t0, Mnet + MELT - Minput_1)
            #Mo2
            Mo2 = np.where((SOILM2t0 * (X2 / h2) - (FWP * (X2 / h2) * h2)) < ((-1)*Mnet - Moutput_1), SOILM2t0 * (X2 / h2) - (FWP * (X2 / h2) * h2), (-1)*Mnet - Moutput_1 * RDR2t0)
    
            #11th-tier variables (mauve/dark purple): L1, Moutput_2, Minput_2 - none saved to file
            #L1 (recall densitySOM is a constant value)
            L1 = 80 * (SOILM1 + (Minput_1 - Moutput_1)) / (h1 * 1000 / densitySOM * 100) / 100 * densitySOM
            #Minput_2
            Minput_2 = np.where(tempmskAr[i] > 0,
                                np.where(Mnet > 0,
                                         np.where(Min2 < 0, 0, Min2), 0), 0)
            #Moutput_2
            Moutput_2 = np.where((SOILM0 == 0) & (Mnet < 0),
                                 np.where(Mo2 < 0, 0, Mo2), 0)
        
            #12th-tier variables (gray-purple): X1max_FDD, X1min_FDD, X1min_TDD, N1, Mo3, Min3, SOILM2
            X1max_FDD = np.where(((24 * (FDDavIN + FDDavt0) / (L1 * (0.5 * RES1)))*100) < 0, 0, ((24 * (FDDavIN + FDDavt0) / (L1 * (0.5 * RES1)))*100))
        
            X1min_FDD = np.where(h1 < X1max_FDD, h1, X1max_FDD)
            
            X1min_TDD = np.where(((24 * (TDDavIN + TDDavt0) / (L1 * (0.5 * RES1)))*100) >= h1, 0, ((24 * (TDDavIN + TDDavt0) / (L1 * (0.5 * RES1)))*100))
        
            N1 = (h1 / 100) * L1 * (RES1 / 2) / 24
        
            Mo3 = np.where((SOILM3t0 * (X3 / h3) - (FWP * (X3 / h3) * h3)) < (((-1)* Mnet) - (Moutput_1 + Moutput_2)), (SOILM3t0 * (X3 / h3) - (FWP * (X3 / h3) * h3)) * RDR3t0, (((-1)*Mnet - (Moutput_1 + Moutput_2)) * RDR3t0))    
    
            Min3 = np.where((((FFC / 1.5) * h3) - SOILM3t0) < ((Mnet + MELT) - Minput_1 - Minput_2), (((FFC / 1.5) * h3) - SOILM3t0), ((Mnet + MELT) - Minput_1 - Minput_2))
        
            #SOILM2 saved to file
            SOILM2 = SOILM2t0 + Minput_2 - Moutput_2
            with rasterio.open(output_path + "\\sm2_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(SOILM2.astype(rasterio.float32))
            
            #13th-tier variables (greenish-yellow): X1in, L2, Minput_3, Moutput_3, RDR2in, PFC2
            #for X1in we have to use np.select because there are too many positional arguments to use a nested np.where, default is the argument for when all conditions evaluate to false
            X1in = np.select(condlist = [(Tffmo > 0) & (X1 < h1) & (TDDavIN >= N1), (Tffmo > 0) & (X1 < h1) & (TDDavIN < N1), (Tffmo <= 0) & (X1 > 0) & (FDDavIN >= N1), (Tffmo <= 0) & (X1 > 0) & (FDDavIN < N1)], 
                             choicelist = [h1, X1min_TDD, 0, h1 - X1min_FDD], 
                             default = X1)
        
            L2 = 80 * (((SOILM2 + Minput_2 + Moutput_2)/h2 * 1000 / densityssc) * 100 )/ 100 * densityssc
        
            Minput_3 = np.where((tempmskAr[i] > 0) & (textmsk < 7), 
                                np.where(Mnet < 0,
                                         np.where(Min3 < 0, 0, Min3), 0), 0)
        
            Moutput_3 = np.where((SOILM0 == 0) & (Mnet < 0), 
                                 np.where(Mo3 < 0, 0, Mo3), 0)
        
            powerRDR2 = (SOILM2/(h2*FFC)) ** m_SAX_B
            RDR2in = (1 + m_SAX_A) / (1 + (m_SAX_A * powerRDR2))
        
            PFC2 = ((SOILM2 + Minput_2 - Moutput_2) / (FFC * h2)) * 100
        
            #14th-tier variables (dark pink): N2_2, X2in_max, SOILM3, RDR2
            N2_2 = (h2/100) * L2 * ((RES1 + RES2)/2) / 24
        
            X2in_max = np.where((((24 * FDDavt0)/(L2 * (RES1 + 0.5 * RES2))) * 100) < 0, 0, ((24 * FDDavt0)/(L2 * (RES1 + 0.5 * RES2))) * 100)
        
            SOILM3 = SOILM3t0 + Minput_3 - Moutput_3
        
            RDR2 = RDR2t0 + RDR2in - RDR2t0
        
            #15th-tier variables (grey): X2in_min, L3, PFC3, RDR3in, EET
            X2in_min = np.where(h2 < X2in_max, h2, X2in_max)
        
            L3 = 80 * ((((SOILM3 + Minput_3 - Moutput_3) / h3) * 1000 / densityssc) * 100) / 100 * densityssc
        
            PFC3 = ((SOILM3 + Minput_3 - Moutput_3) / (FFC / 1.5 * h3)) * 100
        
            powerRDR3in = (SOILM3 / (h3*(FFC/1.5))) ** m_SAX_B
            RDR3in = (1 + m_SAX_A) / (1 + (m_SAX_A * powerRDR3in))
        
            #EET saved to file
            EET = np.where(pptmskAr[i] < PET, 
                           np.where((pptmskAr[i] + (PET - pptmskAr[i]) * RDR2) < (pptmskAr[i] + ((SOILM1 + SOILM2 + SOILM3) - (FWP * (h1 + h2 + h3)))), (pptmskAr[i] + (PET - pptmskAr[i]) * RDR2), 
                                    (pptmskAr[i] + ((SOILM1 + SOILM2 + SOILM3) - (FWP * (h1 + h2 + h3))))), PET)
            with rasterio.open(output_path + "\\eet_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(EET.astype(rasterio.float32))
        
            #16th-tier variables (maroon): X2in, N3, X3in_min1, X3in_max, RDR3, W
            X2in = np.select(condlist = [(tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 >= N2_2), 
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET > 0) & ((0.5*(EET + SOILM1)/PET) < 0.5), 
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET > 0) & ((0.5*(EET + SOILM1)/PET) >= 0.5),
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET <= 0), 
                                         (tempmskAr[i] > 0) & (X1 == 0) & (X2 > 0) & (FDDavt0 >= N2_2),
                                         (tempmskAr[i] <= 0) & (X1 == 0) & (X2 > 0) & (FDDavt0 < N2_2)],
                             choicelist = [h2, 0.5*((EET + SOILM1)/PET), 0.5, 0, 0, h2 - X2in_min],
                             default = X2)
        
            N3 = (h3/100)*L3*(RES1 + RES2 + RES3 / 2) / 24
        
            X3in_min1 = np.where(h3 < (((24*TDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100), h3, ((24*TDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100)
        
            X3in_max = np.where((((24*FDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100) < 0, 0, ((24*FDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100)
        
            #RDR3 - this one does not make any sense
            RDR3 = RDR3t0 + RDR3in - RDR3t0
            #W saved to file
            W_mx1 = np.where((0.5 + 0.5*SOILM1) > 0.5, 0.5 + 0.5*SOILM1, 0.5)
            W_mn1 = np.where((PET > 0) & ((0.5 * ((EET + SOILM1) / PET)) < 0.5), 0.5 * ((EET + SOILM1) / PET),
                             np.where((PET > 0) & ((0.5 * ((EET + SOILM1) / PET)) >= 0.5), 0.5, 0))
            W_mx2 = np.where((0.5 + W_mn1) > 0.5, 0.5 + W_mn1, 0.5)
            W = np.where((PET > 0) & (EET > 0), W_mx2, W_mx1)
            with rasterio.open(output_path + "\\w_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                file.write(W.astype(rasterio.float32))
        
            #17th-tier variables (very light purple): FDDavOUT, TDDavOUT, X3in_min2, NPP
            FDDavOUT = np.where(TDDavIN > 0, FDDavt0, 
                                np.where((Tffmo <= 0) & (X1>0) & (FDDavt0>=N1), N1,
                                         np.where((Tffmo <= 0) & (X1==0) & (X2>0) & (FDDavt0 >= N2_2), N2_2,
                                                  np.where((Tffmo <= 0) & (X2==0) & (X3>0) & (FDDavt0>=N3), N3, 0))))
        
            TDDavOUT = np.where(FDDavIN > 0, TDDavt0,
                                np.where((Tffmo > 0) & (X1<h1) & (TDDavt0 >= N1), N1,
                                         np.where((Tffmo > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 >= N3), N3, 0)))
        
            X3in_min2 = np.where(h3 < X3in_max, h3, X3in_max)
        
            #NPP saved to file
            NPP = np.where(vegmskAr[j-1] > 0, vegmskAr[j-1]*solarmskAr[i]*T1*T2*emax*W, 0)
            with rasterio.open(output_path + "\\npp_" + str(j) + "_" + str(End_year) + ".tif", mode='w', **kwargs) as file:
                   file.write(NPP.astype(rasterio.float32))
        
            #18th-tier variables (blue-grey): FDDav, TDDav, X3in
            FDDav = FDDavt0 + FDDavIN - FDDavOUT
        
            TDDav = TDDavt0 + TDDavIN - TDDavOUT
        
            X3in = np.where((tempmskAr[i] > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 >= N3), h3,
                            np.where((tempmskAr[i] > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 < N3), X3in_min1,
                                     np.where((tempmskAr[i] <= 0) & (X2==0) & (X3>0) & (FDDavt0 >= N3), 0,
                                              np.where((tempmskAr[i] <= 0) & (X2==0) & (X3 > 0) & (FDDavt0 < N3), h3-X3in_min2, X3))))
        
            #constant raster updates - update initial values sequentially for each month
            #recall Minput_0 and Moutput_0, the flooding variables, are not used in this version and so are not updated here
            PACKt0 = PACK
            X1int0 = X1in
            X2int0 = X2in
            X3int0 = X3in
            X1t0 = X1
            X2t0 = X2
            X3t0 = X3
            RDR2t0 = RDR2
            RDR3t0 = RDR3
            TDDavt0 = TDDav
            FDDavt0 = FDDav
            SOILM0t0 = SOILM0
        
        else:
            #no first-year input will be saved to a file
            j = 0
            #X1, X2, X3, SOILM0 calculations are second-tier (orange) variables and will be updated at end of code
            #all masked files are also second-tier variables
            X1 = X1t0 + X1int0 - X1t0
            X2 = X2t0 + X2int0 - X2t0
            X3 = X3t0 + X3int0 - X3t0
            SOILM0 = SOILM0t0 + Minput_0 - Moutput_0
            
            #third-tier (turquoise variables) begin here: SVPmax, SVPmin, SNOW, days, monthlyLAI
            #NEED TO ADD SNOWBD AND SNOWHARD HERE AS WELL
            #create constants to calculate SVPmax
            power1_SVPmax = (0.00738 * tmaxmskAr[i] + 0.8072) ** 8.0
            absolute1_SVPmax = np.abs(1.8 * tmaxmskAr[i] + 48.0)
            SVPmax = 33.8639 * (power1_SVPmax - 0.000019 * absolute1_SVPmax + 0.001316)
            #create constants to calculate SVPmin
            power1_SVPmin = (0.00738 * tminmskAr[i] + 0.8072) ** 8.0
            absolute1_SVPmin = np.abs(1.8 * tminmskAr[i] + 48.0)
            SVPmin = 33.8639 * (power1_SVPmin - 0.000019 * absolute1_SVPmin + 0.001316)
            #SNOW calculation
            SNOW = np.where(tempmskAr[i] < 1, pptmskAr[i], 0)

            #set days constant for first year months
            if i == 3 or i == 5 or i == 8 or i == 10:
                days = 30
            elif i == 1:
                if (Begin_year == 2000 or Begin_year == 2004 or Begin_year == 2008 or Begin_year == 2012 or Begin_year == 2016 or Begin_year == 2020):
                    days = 29
                else:
                    days = 28
            else:
                days = 31
            #set constant LAI values for all of second year months
            if i == 0 or i == 1 or i == 2 or i == 10 or i == 11:
                monthlyLAI = 0
            elif i == 3:
                monthlyLAI = 0.23
            elif i == 4:
                monthlyLAI = 1.68
            elif i == 5:
                monthlyLAI = 3.85
            elif i == 6:
                monthlyLAI = 4
            elif i == 7:
                monthlyLAI = 2.28
            elif i == 8:
                monthlyLAI = 1.08
            else:
                monthlyLAI = 0.16
            
            #fourth tier variables (black) begin here: Tffmo, solar, PT_a, PT_b, MELT
       
            #Tffmo calculation - not saved to a file, create constants that feed Tffmo
            power1_Tffmo = tempmskAr[i] ** 3.0
            #this should yield a constant
            naturallogbasee1_Tffmo = np.log(1 + monthlyLAI)
            #this should yield an array
            naturallogbasee2_Tffmo = np.log(1 + LAImax)
            #had to change this here because the original first argument to test conditions on did not make sense - test with Chris
            Tffmo = np.where(LAImax > monthlyLAI,
                             (tempmskAr[i] + ((-0.11 + 0.96 * tempmskAr[i] - 0.00008 * power1_Tffmo) - tempmskAr[i]) * naturallogbasee1_Tffmo / naturallogbasee2_Tffmo),
                             (tempmskAr[i] + ((-0.11 + 0.96 * tempmskAr[i] - 0.00008 * power1_Tffmo) - tempmskAr[i]) * naturallogbasee2_Tffmo / naturallogbasee2_Tffmo))
        
            #create constants to evaluate PT_a
            power1_PT = (38 - (2*(elev1 / 305)) + 380 / (SVPmax - SVPmin)) ** -1.0
            #PT_a and PT_b calculation
            PT_a = np.where(power1_PT > 0.0001, power1_PT, 0.0001)
            PT_b = 2.5 + 0.14 * (SVPmax - SVPmin) + elev1 / 550
            #MELT calculation
            MELT = np.where((tempmskAr[i] > 0) & (X1t0 > 0), PACKt0, 0)
        
            #5th-tier variables begin here (melon-colored): RES1, RES2, RES3, TDDavIN, FDDavIN, PET, PACK
            #RES1 calculation - recall lamSOM is a constant value
            RES1 = np.where(tempmskAr[i] > 0, h1 / 100 / lamSOM, h1 / 100 / 0.5)
            #RES2 calculation - recall lambdavalue is a constant value
            RES2 = (h2 / 100) / lambdavalue
            #RES3 calculation - recall lambdavalue is a constant value
            RES3 = (h3 / 100) / lambdavalue
            #TDDavIN calculation
            TDDavIN = np.where(ETIM == 0, 
                               np.where(Tffmo > 0, Tffmo * float(days), 0), 
                               np.where(Tffmo > 0, Tffmo, 0))
            #FDDavIN calculation
            FDDavIN = np.where(ETIM == 0, 
                               np.where(Tffmo <= 0, Tffmo *(-1) * float(days), 0), 
                               np.where(Tffmo <= 0, Tffmo * (-1), 0))
            #Process: Potential Evapotranspiration Calculation, PET-Priestley-Taylor equations uses Tmx and Tmn
            PET = np.where(((PT_a * (tempmskAr[i] + PT_b) * solarmskAr[i] * 36.16704) / (-9.9255 * tempmskAr[i] + 10468.9636) * float(days)) > 0, 
                           (PT_a * (tempmskAr[i] + PT_b) * solarmskAr[i] * 36.16704)/(-9.9255 * tempmskAr[i] + 10468.9636) * float(days), 0)
            #PACK calculation
            PACK = PACKt0 + SNOW - MELT
        
            #sixth-tier variables (brown) begin here: PPTint, DepthSP
            #will need to add DepthSP later when add SnowBd and SnowHard
            #PPTint calculation - not saved to a file
            PPTint = np.where(PET < pptmskAr[i], 0.05 * monthlyLAI * PET, 0.05 * monthlyLAI * pptmskAr[i])

            #7th-tier variables (purple) begin here: Mnet
            #Mnet calculation
            Mnet = np.where((tempmskAr[i] > 0) & (X1t0 > 0), pptmskAr[i] - PPTint + MELT - PET, 0)
        
            #8th-tier variables (green) begin here: Mo1, Min1
            #Min1 calculation
            Min1 = np.where((0.2 * h1 - SOILM1t0) < (Mnet + MELT), 0.2 * h1 - SOILM1t0, Mnet + MELT)
            #Mo1 calculation
            absolute_Mo1 = np.abs(Mnet)
            Mo1 = np.where((SOILM1t0 * (X1 / h1) - (FWP * (X1 / h1) * h1)) < absolute_Mo1, SOILM1t0 * (X1 / h1) - (FWP * (X1 / h1) * h1), absolute_Mo1)
        
            #9th-tier variables (pink) begin here: Moutput_1, Minput_1
            #Minput_1
            Minput_1 = np.where(tempmskAr[i] > 0, 
                                np.where(Mnet > 0,
                                         np.where(Min1 < 0, 0, Min1), 0), 0)
            #Moutput_1
            Moutput_1 = np.where(SOILM0 == 0, Mo1, 0)
    
            #10th-tier variables (darker blue-green): SOILM1, Mo2, Min2
            #SOILM1 - saved to file
            SOILM1 = SOILM1t0 + Minput_1 - Moutput_1
        
            #Min2 
            Min2 = np.where(((FFC * h2) - SOILM2t0) < ((Mnet + MELT) - Minput_1), FFC * h2 - SOILM2t0, Mnet + MELT - Minput_1)
            #Mo2
            Mo2 = np.where((SOILM2t0 * (X2 / h2) - (FWP * (X2 / h2) * h2)) < ((-1)*Mnet - Moutput_1), SOILM2t0 * (X2 / h2) - (FWP * (X2 / h2) * h2), (-1)*Mnet - Moutput_1 * RDR2t0)
    
            #11th-tier variables (mauve/dark purple): L1, Moutput_2, Minput_2 - none saved to file
            #L1 (recall densitySOM is a constant value)
            L1 = 80 * (SOILM1 + (Minput_1 - Moutput_1)) / (h1 * 1000 / densitySOM * 100) / 100 * densitySOM
            #Minput_2
            Minput_2 = np.where(tempmskAr[i] > 0,
                                np.where(Mnet > 0,
                                         np.where(Min2 < 0, 0, Min2), 0), 0)
            #Moutput_2
            Moutput_2 = np.where((SOILM0 == 0) & (Mnet < 0),
                                 np.where(Mo2 < 0, 0, Mo2), 0)
        
            #12th-tier variables (gray-purple): X1max_FDD, X1min_FDD, X1min_TDD, N1, Mo3, Min3, SOILM2
            X1max_FDD = np.where(((24 * (FDDavIN + FDDavt0) / (L1 * (0.5 * RES1)))*100) < 0, 0, ((24 * (FDDavIN + FDDavt0) / (L1 * (0.5 * RES1)))*100))
        
            X1min_FDD = np.where(h1 < X1max_FDD, h1, X1max_FDD)
        
            X1min_TDD = np.where(((24 * (TDDavIN + TDDavt0) / (L1 * (0.5 * RES1)))*100) >= h1, 0, ((24 * (TDDavIN + TDDavt0) / (L1 * (0.5 * RES1)))*100))
        
            N1 = (h1 / 100) * L1 * (RES1 / 2) / 24
        
            Mo3 = np.where((SOILM3t0 * (X3 / h3) - (FWP * (X3 / h3) * h3)) < (((-1)* Mnet) - (Moutput_1 + Moutput_2)), (SOILM3t0 * (X3 / h3) - (FWP * (X3 / h3) * h3)) * RDR3t0, (((-1)*Mnet - (Moutput_1 + Moutput_2)) * RDR3t0))    
    
            Min3 = np.where((((FFC / 1.5) * h3) - SOILM3t0) < ((Mnet + MELT) - Minput_1 - Minput_2), (((FFC / 1.5) * h3) - SOILM3t0), ((Mnet + MELT) - Minput_1 - Minput_2))
        
            #SOILM2 saved to file
            SOILM2 = SOILM2t0 + Minput_2 - Moutput_2
        
            
            #13th-tier variables (greenish-yellow): X1in, L2, Minput_3, Moutput_3, RDR2in, PFC2
            X1in = np.select(condlist = [(Tffmo > 0) & (X1 < h1) & (TDDavIN >= N1), (Tffmo > 0) & (X1 < h1) & (TDDavIN < N1), (Tffmo <= 0) & (X1 > 0) & (FDDavIN >= N1), (Tffmo <= 0) & (X1 > 0) & (FDDavIN < N1)], 
                             choicelist = [h1, X1min_TDD, 0, h1 - X1min_FDD], 
                             default = X1)
        
            L2 = 80 * (((SOILM2 + Minput_2 + Moutput_2)/h2 * 1000 / densityssc) * 100 )/ 100 * densityssc
        
            Minput_3 = np.where((tempmskAr[i] > 0) & (textmsk < 7), 
                                np.where(Mnet < 0,
                                         np.where(Min3 < 0, 0, Min3), 0), 0)
        
            Moutput_3 = np.where((SOILM0 == 0) & (Mnet < 0), 
                                 np.where(Mo3 < 0, 0, Mo3), 0)
        
            powerRDR2 = (SOILM2/(h2*FFC)) ** m_SAX_B
            RDR2in = (1 + m_SAX_A) / (1 + (m_SAX_A * powerRDR2))
        
            PFC2 = ((SOILM2 + Minput_2 - Moutput_2) / (FFC * h2)) * 100
        
            #14th-tier variables (dark pink): N2_2, X2in_max, SOILM3, RDR2
            N2_2 = (h2/100) * L2 * ((RES1 + RES2)/2) / 24
        
            X2in_max = np.where((((24 * FDDavt0)/(L2 * (RES1 + 0.5 * RES2))) * 100) < 0, 0, ((24 * FDDavt0)/(L2 * (RES1 + 0.5 * RES2))) * 100)
        
            SOILM3 = SOILM3t0 + Minput_3 - Moutput_3
        
            RDR2 = RDR2t0 + RDR2in - RDR2t0
        
            #15th-tier variables (grey): X2in_min, L3, PFC3, RDR3in, EET
            X2in_min = np.where(h2 < X2in_max, h2, X2in_max)
        
            L3 = 80 * ((((SOILM3 + Minput_3 - Moutput_3) / h3) * 1000 / densityssc) * 100) / 100 * densityssc
        
            PFC3 = ((SOILM3 + Minput_3 - Moutput_3) / (FFC / 1.5 * h3)) * 100
        
            powerRDR3in = (SOILM3 / (h3*(FFC/1.5))) ** m_SAX_B
            RDR3in = (1 + m_SAX_A) / (1 + (m_SAX_A * powerRDR3in))
        
            #EET
            EET = np.where(pptmskAr[i] < PET, 
                           np.where((pptmskAr[i] + (PET - pptmskAr[i]) * RDR2) < (pptmskAr[i] + ((SOILM1 + SOILM2 + SOILM3) - (FWP * (h1 + h2 + h3)))), (pptmskAr[i] + (PET - pptmskAr[i]) * RDR2), 
                                    (pptmskAr[i] + ((SOILM1 + SOILM2 + SOILM3) - (FWP * (h1 + h2 + h3))))), PET)
        
            #16th-tier variables (maroon): X2in, N3, X3in_min1, X3in_max, RDR3, W
            X2in = np.select(condlist = [(tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 >= N2_2), 
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET > 0) & ((0.5*(EET + SOILM1)/PET) < 0.5), 
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET > 0) & ((0.5*(EET + SOILM1)/PET) >= 0.5),
                                         (tempmskAr[i] > 0) & (X1==h1) & (X2 < h2) & (TDDavt0 < N2_2) & (PET <= 0), 
                                         (tempmskAr[i] > 0) & (X1 == 0) & (X2 > 0) & (FDDavt0 >= N2_2),
                                         (tempmskAr[i] <= 0) & (X1 == 0) & (X2 > 0) & (FDDavt0 < N2_2)],
                             choicelist = [h2, 0.5*((EET + SOILM1)/PET), 0.5, 0, 0, h2 - X2in_min],
                             default = X2)
        
            N3 = (h3/100)*L3*(RES1 + RES2 + RES3 / 2) / 24
        
            X3in_min1 = np.where(h3 < (((24*TDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100), h3, ((24*TDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100)
        
            X3in_max = np.where((((24*FDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100) < 0, 0, ((24*FDDavt0)/(L3 * (RES1 + RES2 + 0.5 * RES3)))* 100)
        
            #RDR3 - this one does not make any sense
            RDR3 = RDR3t0 + RDR3in - RDR3t0
            #W 
            W_mx1 = np.where((0.5 + 0.5*SOILM1) > 0.5, 0.5 + 0.5*SOILM1, 0.5)
            W_mn1 = np.where((PET > 0) & ((0.5 * ((EET + SOILM1) / PET)) < 0.5), 0.5 * ((EET + SOILM1) / PET),
                             np.where((PET > 0) & ((0.5 * ((EET + SOILM1) / PET)) >= 0.5), 0.5, 0))
            W_mx2 = np.where((0.5 + W_mn1) > 0.5, 0.5 + W_mn1, 0.5)
            W = np.where((PET > 0) & (EET > 0), W_mx2, W_mx1)
        
            #17th-tier variables (very light purple): FDDavOUT, TDDavOUT, X3in_min2, NPP
            FDDavOUT = np.where(TDDavIN > 0, FDDavt0, 
                                np.where((Tffmo <= 0) & (X1>0) & (FDDavt0>=N1), N1,
                                         np.where((Tffmo <= 0) & (X1==0) & (X2>0) & (FDDavt0 >= N2_2), N2_2,
                                                  np.where((Tffmo <= 0) & (X2==0) & (X3>0) & (FDDavt0>=N3), N3, 0))))
        
            TDDavOUT = np.where(FDDavIN > 0, TDDavt0,
                                np.where((Tffmo > 0) & (X1<h1) & (TDDavt0 >= N1), N1,
                                         np.where((Tffmo > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 >= N3), N3, 0)))
        
            X3in_min2 = np.where(h3 < X3in_max, h3, X3in_max)
        
            #18th-tier variables (blue-grey): FDDav, TDDav, X3in
            FDDav = FDDavt0 + FDDavIN - FDDavOUT
        
            TDDav = TDDavt0 + TDDavIN - TDDavOUT
        
            X3in = np.where((tempmskAr[i] > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 >= N3), h3,
                            np.where((tempmskAr[i] > 0) & (X2==h2) & (X3 < h3) & (TDDavt0 < N3), X3in_min1,
                                     np.where((tempmskAr[i] <= 0) & (X2==0) & (X3>0) & (FDDavt0 >= N3), 0,
                                              np.where((tempmskAr[i] <= 0) & (X2==0) & (X3 > 0) & (FDDavt0 < N3), h3-X3in_min2, X3))))
        
            #constant raster updates - update initial values sequentially for each month
            #recall Minput_0 and Moutput_0, the flooding variables, are not used in this version and so are not updated here
            PACKt0 = PACK
            X1int0 = X1in
            X2int0 = X2in
            X3int0 = X3in
            X1t0 = X1
            X2t0 = X2
            X3t0 = X3
            RDR2t0 = RDR2
            RDR3t0 = RDR3
            TDDavt0 = TDDav
            FDDavt0 = FDDav
            SOILM0t0 = SOILM0
        


