#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:53:59 2023

@author: x_tilin
"""

import os
import glob
import netCDF4 
import numpy as np
from datetime import datetime

##########################
### load era5 surface data     ###
##########################

era_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project3/pl_mcao1/backup/'
os.chdir(era_path)
sl_file= glob.glob("ERA5-20141215-20141231-sl.nc", recursive=True)
sl_file=''.join(sl_file)
sl_file = netCDF4.Dataset(sl_file)
#sl_file.variables.keys()

###############################
###set up era5 time for loop###
###############################

era_time=sl_file.variables['time'][:]
era_time = "2014-12-15 00:00:00"     ##refer to era5 data
era_time = datetime.strptime(era_time, "%Y-%m-%d %H:%M:%S")
era_time = era_time.timestamp()


start_time = "2014-12-25 00:00:00"     ##refer to era5 data
start_time = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
start_time = start_time.timestamp()

# test_time=era_time+24*60*60*10
# test_time = datetime.fromtimestamp(test_time)
# test_time = test_time.strftime("%Y-%m-%d %H:%M:%S")

end_time = "2014-12-31 23:00:00"     ##refer to era5 data
end_time = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
end_time = end_time.timestamp()

i_start=int((start_time-era_time)/60/60)
i_end=int((end_time-era_time)/60/60)


##########################
###load cmip6 data     ###
##########################

cmip6_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/cmip6/cmip_25km/'
os.chdir(cmip6_path)
cmip6_file= glob.glob("ERA5-20141215-20141231-sl.nc", recursive=True)
cmip6_file=''.join(cmip6_file)
cmip6_file = netCDF4.Dataset(cmip6_file)



for i in range(i_start,i_end):
    
    era_sst=sl_file.variables['sst'][i]
    era_lon=sl_file.variables['lon'][:]
    era_lat=sl_file.variables['lat'][:]














