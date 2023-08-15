#python NCEP_nc_to_tiff_conv.py 

#This script is modified from the original CNEPconv_FLT_Tiff.py that was created by Cyrus Hiatt in March, 2012, Modified by: Vanessa Genovese May, 2017, and Modified by Roshni Biswas Feb, 2019.

#This script reads monthly NCEP NetCDF files from NOAA and converts them to TIFF. Manually input a year and the script calculates what subdatsets are needed for that year. The 12 nc files for one year are processsed into tiff files and saved into a staging folder. 
#Then the Tiffs are opened back up for necessary conversions, reprojected into Mollweide, given an ocean mask, given a 20x20 rectangular focal filter, and saved into a new folder. 

#the solar conversions are currently set for daylight hours in the Nothern Hemisphere. This needs to be changed when looking at other parts of the world.

#Updated and modified by: Davina Dou November, 2022
#Updated and modified by: Stephanie Pass March, 2023
#-------------------------------------------------------------------------------------------------------------------

import arcpy, sys, string
from arcpy.sa import *
from arcpy import env
from arcpy.ia import *
import os
from osgeo import gdal, osr, gdal_array
from scipy.io import netcdf
import numpy as np
from datetime import date, timedelta

#Manually write in the in_file, output_dir, and year desired
#Clear the staging folder before processing a new NC file and year
nc_file = "D:\\NCEP_Monthly_nc_files\\dswrf.sfc.mon.mean.nc"
#use one of these for the nc file
#air.2m.mon.mean.nc
#prate.sfc.mon.mean.nc
#tmin.2m.mon.mean.nc
#tmax.2m.mon.mean.nc
#dswrf.sfc.mon.mean.nc
output_dir = "D:\\CASA_inputs_global\\NCEP_files_tiffs\\"
year = 2016

def get_subdataset_ids(year):
    #This works for monthly NCEP NetCDF files from NOAA that start at January 1979 as the 0th file
    #This def will output the 12 subdatsets to be processed for the inputted year
    start_year = int(1979)
    end_year = int(year - 1979)
    start_subs_needed = int(end_year * 12)
    subdataset_list = list(range(start_subs_needed, start_subs_needed + 12))
    return subdataset_list

    
subdataset_id_list = get_subdataset_ids(int(year))
print(subdataset_id_list)


def process_subdataset(subdataset, output_file_name):
    #This unpacks a single subdataset and writes it to Geotiff format

    # read into numpy array
    band_array = subdataset.ReadAsArray().astype(np.float32)

    # convert no_data values
    band_array[band_array == -28672] = -32768

    # write raster
    out_ds = gdal.GetDriverByName('GTiff').Create(output_file_name,
                                                  subdataset.RasterXSize,
                                                  subdataset.RasterYSize,
                                                  1,  #Number of bands
                                                  gdal.GDT_Float32,
                                                  ['COMPRESS=NONE', 'TILED=YES'])
    #Our version of ArcGIS Pro(2.5) and GDAL take GeoTransform numbers from(0.0, 1.0, 0.0, 0.0, 0.0, 1.0) and projects them correctly starting at the dateline with these geotransform numbers.
    out_ds.SetGeoTransform([0, 1.875, 0, 90, 0, -1.875])
    #Inputs the NetCDF coordinate system so it projects correctly
    out_ds_srs = osr.SpatialReference()
    out_ds_srs.ImportFromEPSG(4326)
    out_ds.SetProjection(out_ds_srs.ExportToWkt())
    out_ds.GetRasterBand(1).WriteArray(band_array)
    out_ds.GetRasterBand(1).SetNoDataValue(-32768)
   

    out_ds = None  #close dataset to write to disc

    return output_file_name

def process_nc_file(nc_file, output_dir, subdataset_id_list, year):
    #this def takes the 12 subdatasets, gets the variable name, renames the files, runs them through the process subdataset function, and saves them in the staging folder.
    nc_ds = gdal.Open(nc_file)
    subdataset_list = nc_ds.GetSubDatasets()
    output_tif_names = []
    # get the variable
    variable = os.path.basename(nc_file.split('.')[0])
    print(variable)
    #get the variable into CASA format
    if variable == 'tmin':
        newvar = 'mintemp'
    elif variable == 'tmax':
        newvar = 'maxtemp'
    elif variable == 'prate':
        newvar = 'ppt'
    elif variable == 'dswrf':
        newvar = 'solar'
    elif variable == 'air':
        newvar = 'meantemp'
    else:
            print('Error, variable is not recognized')
    #this runs the 12 subdatsets in a loop through the process_subdataset function
    for i in range(len(subdataset_id_list)):
        subdataset = gdal.Open(subdataset_list[subdataset_id_list[i]][0],gdal.GA_ReadOnly)
        #confirm correct dates
        confirm_tim_var = subdataset_list[subdataset_id_list[i]][1]
        get_time_var = confirm_tim_var.split()
        time_var = get_time_var[4]
        time_stamp = time_var.split('"')
        hours = time_stamp[1]
        days = int(hours)/24
        start = date(1800,1,1) 
        delta = timedelta(days)
        offset = start + delta
        print(offset.strftime('%Y-%m-%d %H:%M:%S'))
        #renames the files and places them in the staging folder to be saved
        month = ''
        if i + 1 == 1:
            month = 'jan'
        elif i + 1 == 2:
            month = 'feb'
        elif i + 1 == 3:
            month = 'mar'
        elif i + 1 == 4:
            month = 'apr'
        elif i + 1 == 5:
            month = 'may'
        elif i + 1 == 6:
            month = 'jun'
        elif i + 1 == 7:
            month = 'jul'
        elif i + 1 == 8:
            month = 'aug'
        elif i + 1 == 9:
            month = 'sep'
        elif i + 1 == 10:
            month = 'oct'
        elif i + 1 == 11:
            month = 'nov'
        else:
            month = 'dec'
        output_file_name = output_dir + '\stg\\' + newvar + '_' + str(year) + '.' + month  + "_" + str(i+1).zfill(2) + '.tif'
        output_tif_names.append(process_subdataset(subdataset, output_file_name))
    return output_tif_names

process = process_nc_file(nc_file, output_dir, subdataset_id_list, year)

#the next part of the script opens up the tiff files from the staging folder and reprojects them

#check spatial license
arcpy.CheckOutExtension("Spatial")
input_directory = "D:\\CASA_inputs_global\\NCEP_files_tiffs\\stg"

def reproject_raster(raster, casa_variable):
    #The tiffs are opened back up for necessary conversions, reprojected into Mollweide, given an ocean mask, given a 20x20 rectangular focal filter, and saved into a new folder. 
    output_dir = 'D:\\CASA_inputs_global\\NCEP_files_tiffs\\'
    arcpy.env.overwriteOutput = True
    #applythemask
    nulltemplate_raster = "D:\\CASA_inputs_global\\NCEP_files_tiffs\\NCEP_ocean_mask\\SetNull_Dem8KM1.tif"
    arcpy.env.snapRaster = nulltemplate_raster
    arcpy.env.extent = nulltemplate_raster
    arcpy.env.cellSize = nulltemplate_raster
    newraster = raster[:-7] + ".tif"
    #reads the month substring of the file to do solar conversions later
    month = newraster[11:14]
    print(month)
    #set the out raster projection to mollweide
    out_sr = arcpy.CreateSpatialReference_management(54009)
    sr = arcpy.SpatialReference(54009)
    #The input projection is WGS 84 EPSG 4326
    insr = arcpy.SpatialReference(4326)
    #reproject raster to mollweide with 8000 size cells
    projection = arcpy.ProjectRaster_management(raster , newraster, sr, "NEAREST", 8000 , "#", "#", insr)
    # apply the focal statistics to the newly projected mollweide raster
    focal_raster = FocalStatistics(projection, NbrRectangle(20, 20, "CELL") , "MEAN", "DATA")
    #Apply the ocean mask to the raster
    focal_ext = ExtractByMask(focal_raster, nulltemplate_raster)
    #do necessary conversions
    if casa_variable in ('mintemp', 'maxtemp', 'meantemp'):
        arcpy.env.overwriteOutput = True
        from_unit = 'Kelvin'
        to_unit = 'Celsius'
        out_temp_raster = UnitConversion(focal_ext, from_unit, to_unit)
        out_temp_raster.save(output_dir + newraster)
        #close out the functions
        projection = None
        focal_raster = None
        focal_ext = None
        out_temp_raster = None
        return newraster
    elif casa_variable in ('ppt', 'prate'):
        arcpy.env.overwriteOutput = True
        inconstant = 262974.6
        out_prate_raster = Times(focal_ext, inconstant)
        out_prate_raster.save(output_dir + newraster)
        #close out the functions
        projection = None
        focal_raster = None
        focal_ext = None
        out_temp_raster = None
        out_prate_raster = None
        return newraster
    elif casa_variable in ('solar', 'dswrf'):         
        arcpy.env.overwriteOutput = True
        #calculate inconstant for solar for the northern hemisphere
        #change the ave_daylight_min for others parts of the world
        if month == 'jan':
            ave_daylight_min = 580
        elif month == 'feb':
            ave_daylight_min = 622
        elif month == 'mar':
            ave_daylight_min = 683
        elif month == 'apr':
            ave_daylight_min = 757
        elif month == 'may':
            ave_daylight_min = 825
        elif month == 'jun':
            ave_daylight_min = 874
        elif month == 'jul':
            ave_daylight_min = 881
        elif month == 'aug':
            ave_daylight_min = 842
        elif month == 'sep':
            ave_daylight_min = 778
        elif month == 'oct':
            ave_daylight_min = 707
        elif month == 'nov':
            ave_daylight_min = 637
        elif month == 'dec':
            ave_daylight_min = 587
        else:
            print('month file format not recognized')
        #1440 minutes in a day
        daymin = 1440
        #average seconds in a month
        ave_seconds = 2592000
        daylightminutes = ave_daylight_min/daymin
        #calculate MJ per month
        solar_inconstant = daylightminutes * ave_seconds * (10**-6)
        print(solar_inconstant)
        out_solar_raster = Times(focal_ext, solar_inconstant)
        out_solar_raster.save(output_dir + newraster)
        #close out the functions
        projection = None
        focal_raster = None
        focal_ext = None
        out_temp_raster = None
        out_prate_raster = None
        out_solar_raster = None
        return newraster
    else:
        print('casa_variable is not recognized')

def reproject_all_rasters(all_rasters , casa_variable):
    #this runs all 12 rasters through the reproject raster function
    arcpy.env.workspace = input_directory
    arcpy.env.overwriteOutput = True
    rasters = arcpy.ListRasters("*.tif")
    #print(rasters)
    #variable = rasters[0].split('.')[0] # get the variable
    print(casa_variable)
    for raster in rasters:
        reproject_raster(raster , casa_variable)

#this informaiton is to get the variable info for the reprojected raster and then runs the reproject_all_rasters function
arcpy.env.workspace = input_directory
my_rasters = arcpy.ListRasters("*.tif")
print(my_rasters)
casa_variable = my_rasters[0].split('_')[0]

save_off_files = reproject_all_rasters(my_rasters, casa_variable)
