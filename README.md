# NASA-CASA
This repository contains the scripts used to prepare and run CASA (The Carnegie-Ames-Stanford Approach) terrestrial carbon cycle model. 

**Scripts Used**

1.CASAglobal.py

2.CASAPrepGlobalEVI_2023.py

3.NCEP_nc_to_tiff_conv.py

**1. CASAglobal.py**
   
This is the script that takes the prepared files and runs the terrestrial carbon cycle model

**3. CASAPrepGlobalEVI_2023.py**

This script prepares MODIS satellite imagery for the CASA global script. It extracts a user-specified subdataset from a mod13C2 hdf file, resamples it to 8km, and reprojects the file into Mollweide. Then saves EVI and EVI * 10 (LAI) to an output directory.

**3. NCEP_nc_to_tiff_conv.py**

The script is built to open NetCDF nc files from NOAA’s database. These nc files contain
monthly sub datasets for different variables since 1979. The script takes a manually inputted
year and outputs what sub datasets are needed to be processed for that year. The 12 nc files
are processed with GDAL into tiff files and saved into a staging folder. Then the TIFFS are opened back up for necessary conversions, reprojected into Mollweide, given an ocean mask, given a 20x20 rectangular focal filter, and saved into a new folder. A separate text file is provided in the folder listing the specific libraries that were used to run this script. 

Select these 5 NC files from the NOAA website and run each seperately through the script: 

air.2m.mon.mean.nc

tmax.2m.mon.mean.nc

tmin.2m.mon.mean.nc

dswrf.sfc.mon.mean.nc

prate.sfc.mon.mean.nc

These 5 NC files can be found from a list here: https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/

