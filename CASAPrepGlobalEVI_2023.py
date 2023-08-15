"""
Extracts a user specified subdataset from a mod13C2 hdf file, resamples it to 8km or 0.5 degree based on the pixel size and coordinate system of a mask file.
Saves EVI and EVI * 10 (LAI) to an output directory. Use setnull Land water mask.
Written by Cyrus Hiatt March 2012
Updated by Vanessa Genovese May,2017
Updated by Roshni Biswas March, 2019
Updated by Stephanie Pass March, 2023
"""

import sys
import os
import datetime
import arcpy
from arcpy import env
import glob
from arcpy.sa import *

input_directory = 'D:\\CASA_inputs_global\\EVI_MODIS\\EVI_MODIS_hdf'
sub_data_set = 1
mask = "D:\\CASA_inputs_global\\Elevation\\Dem8KM.tif"
outputFolder = 'D:\\CASA_inputs_global\\EVI_MODIS\\EVI_MODIS_testfiles\\'

arcpy.CheckOutExtension("Spatial")
maskXmin = arcpy.GetRasterProperties_management(mask, "LEFT")
maskYmin = arcpy.GetRasterProperties_management(mask, "BOTTOM")
maskXmax = arcpy.GetRasterProperties_management(mask, "RIGHT")
maskYmax = arcpy.GetRasterProperties_management(mask, "TOP")
maskCellSize = arcpy.GetRasterProperties_management(mask, "CELLSIZEY")
arcpy.Extent = str(maskXmin) +" "+ str(maskYmin) +" "+ str(maskXmax) +" "+ str(maskYmax) +" "+ mask #mask serves as our snap raster.
arcpy.env.workspace = input_directory

#create spatial references for reprojecting
out_sr = arcpy.CreateSpatialReference_management(54009)
sr = arcpy.SpatialReference(54009)
insr = arcpy.SpatialReference(54009)

rasters = arcpy.ListRasters("*.hdf")

for raster in rasters:
    arcpy.env.overwriteOutput = True
    rasterNameIndex = raster.split('.', 5)
    newRaster = input_directory + "\\" + rasterNameIndex[0] + "." + rasterNameIndex[1] + "." + rasterNameIndex[2] + "." + rasterNameIndex[3] + "_sd" + ".tif"
    newRaster2 = newRaster[:-7] + ".tif"
    subdataset = arcpy.ExtractSubDataset_management(raster, newRaster, "1")
    #rasterProjectionText = str(arcpy.CreateSpatialReference_management("World_Mollweide", newRaster))
    
    ordinalDay = int(rasterNameIndex[1][5:8])
    year = rasterNameIndex[1][1:5]
    if ordinalDay < 31:
        calendarMonth = "jan"
    elif (ordinalDay > 31) & (ordinalDay < 59):
        calendarMonth = "feb"
    elif (ordinalDay > 59) & (ordinalDay < 90):
        calendarMonth = "mar"
    elif (ordinalDay > 90) & (ordinalDay < 120):
        calendarMonth = "apr"
    elif (ordinalDay > 120) & (ordinalDay < 151):
        calendarMonth = "may"
    elif (ordinalDay > 151) & (ordinalDay < 181):
        calendarMonth = "jun"
    elif (ordinalDay > 181) & (ordinalDay < 212):
        calendarMonth = "jul"
    elif (ordinalDay > 212) & (ordinalDay < 243):
        calendarMonth = "aug"
    elif (ordinalDay > 243) & (ordinalDay < 273):
        calendarMonth = "sep"
    elif (ordinalDay > 273) & (ordinalDay < 304):
        calendarMonth = "oct"
    elif (ordinalDay > 304) & (ordinalDay < 334):
        calendarMonth = "nov"
    else:
        calendarMonth = "dec"
        
    rasterName = year + '.' + calendarMonth + '.tif'   
    reproject = arcpy.ProjectRaster_management(subdataset, newRaster2, sr , "BILINEAR", 8000, "", "", insr)


    evi = Times(reproject, .0001)
    lai = Times(reproject, .001)

    evi.save(outputFolder + "\\evi_" + rasterName)
    lai.save(outputFolder + "\\lai_" + rasterName)
    
    #close out functions
    reproject = None
    evi = None
    lai = None
    subdataset = None
