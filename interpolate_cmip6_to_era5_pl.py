#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 17:10:31 2023

@author: x_tilin
"""


import os
import glob
import netCDF4 
import numpy as np
from datetime import datetime
import pywinter.winter as pyw
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def creat_FILI_for_WPS(lat,lon,u10m,v10m,tp2m,lsm,msl,ci,sst,skt,sd,sp,                
                       temp,uwind,vwind,rh,z,
                       st1,st2,st3,st4,swvl1,swvl2,swvl3,swvl4): 
    ##Geo
    # lat = sl_file.variables['lat'][:] # degrees north
    # lon = sl_file.variables['lon'][:] # degrees east
    dlat = -np.abs(lat[1] - lat[0])
    dlon = np.abs(lon[1] - lon[0])
    winter_geo = pyw.Geo0(lat[0],lon[0],dlat,dlon)
    
    ##2D##
    #read 2D data
    # i=2
    # u10m=sl_file.variables['10u'][i]
    # v10m=sl_file.variables['10v'][i]
    # dtp2m=sl_file.variables['2d'][i]
    # tp2m=sl_file.variables['2t'][i]
    # lsm=sl_file.variables['lsm'][i]
    # msl=sl_file.variables['msl'][i]
    # ci=sl_file.variables['ci'][i]
    # sst=sl_file.variables['sst'][i]
    # skt=sl_file.variables['skt'][i]
    # sd=sl_file.variables['sd'][i]*1000  ## refer to WPS/ungrib/src/SNOW_EC
    # sp=sl_file.variables['sp'][i]
    
    #creat 2D field
    # winter_snow = pyw.V2d('SNOWH',sd,'Physical Snow Depth','m','200100')
    winter_snow = pyw.V2d('SNOW',sd*1000,'Water Equivalent of Accumulated Snow Depth','kg/m^2','200100') ## refer to WPS/ungrib/src/SNOW_EC, so need to*1000
    winter_seaice = pyw.V2d('SEAICE',ci,' Sea-Ice Fraction','fraction','200100')
    winter_skt = pyw.V2d('SKINTEMP',skt,' Sea-Surface Temperature','K','200100')
    winter_pmsl = pyw.V2d('PMSL',msl,'Sea_level Pressure','Pa','200100')
    winter_lsm = pyw.V2d('LANDSEA',lsm,'Land/Sea flag','0/1 Flag','200100')
    winter_sp = pyw.V2d('PSFC',sp,'Surface Pressure','Pa','200100')
    winter_t2m = pyw.V2d('TT',tp2m,'Temperature','K','200100')
    winter_sst = pyw.V2d('SST',skt,'Sea-Surface Temperature','K','200100')  ##have not downloaded cmip6 sst
    winter_u10 = pyw.V2d('UU',u10m)
    winter_v10 = pyw.V2d('VV',v10m)

    
    ###3D###
    #read 3D data
    # temp = pl_file.variables['t'][i,:,:,:]
    # uwind = pl_file.variables['u'][i,:,:,:]
    # vwind = pl_file.variables['v'][i,:,:,:]
    # rh = pl_file.variables['r'][i,:,:,:]
    # z = pl_file.variables['z'][i,:,:,:]
    # plevs=np.array([100000,  92500,  85000,  70000,  60000,  50000,  40000,  30000,  25000,
    #   20000,  15000,  10000,   7000,   5000,   3000,   2000,   1000,    500,
    #     100])
    plevs=np.array([100, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000,
                    50000, 60000, 70000, 85000, 92500, 100000])
           
    # Create winter 3D isobaric fields
    winter_t = pyw.V3dp('TT',temp,plevs)
    winter_u = pyw.V3dp('UU',uwind,plevs)
    winter_v = pyw.V3dp('VV',vwind,plevs)
    winter_rh = pyw.V3dp('RH',rh,plevs) # rh = 1.E2 * (p*q/(q*(1.-eps) + eps))/(svp1*exp(svp2*(t-svpt0)/(T-svp3)))
    winter_z = pyw.V3dp('GHT',z/9.81,plevs)  ##/WPS/ungrib/src/rrpr.F
    
    ###Soil###
    #read soil data
    # st1=sl_file.variables['stl1'][i]
    # st2=sl_file.variables['stl2'][i]
    # st3=sl_file.variables['stl3'][i]
    # st4=sl_file.variables['stl4'][i]
    soilt_lay=np.squeeze(np.array([st1,st2,st3,st4]))
    # swvl1=sl_file.variables['swvl1'][i]
    # swvl2=sl_file.variables['swvl2'][i]
    # swvl3=sl_file.variables['swvl3'][i]
    # swvl4=sl_file.variables['swvl4'][i]
    soilm_lay=np.squeeze(np.array([swvl1,swvl2,swvl3,swvl4]))
    slt_layer = ['000007','007028','028100','100255']
    
    # Create winter 3D soil fields
    winter_soilt_layer = pyw.Vsl('ST',soilt_lay,slt_layer)
    winter_soilm_layer = pyw.Vsl('SM',soilm_lay,slt_layer)
    
    
    # Listing fields
    total_fields = [winter_t2m, winter_u10, winter_v10, winter_sst,
                    winter_lsm, winter_sp, winter_pmsl, winter_skt, winter_seaice, 
                    winter_snow, winter_soilt_layer, winter_soilm_layer,
                    winter_z, winter_t, winter_u, winter_v, winter_rh]
    
    return winter_geo,total_fields

# Out path
# path_out = '/home/x_tilin/snic2021-23-400/users/x_tilin/run/project3/test/wps/pywinter/'
# pyw.cinter('FILE','2014-12-25_02',winter_geo,total_fields,path_out)

###main programma###

##########################
### load era5 surface and plev data     ###
##########################

era_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project3/pl_Erik2002/'
os.chdir(era_path)
pl_era5= glob.glob("ERA5-20021215-20021225-pl.nc", recursive=True)
pl_era5=''.join(pl_era5)
pl_era5 = netCDF4.Dataset(pl_era5)
era_lon=pl_era5.variables['lon'][:]
era_lat=pl_era5.variables['lat'][:]
#pl_era5.variables.keys()

sl_era5= glob.glob("ERA5-20021215-20021225-sl.nc", recursive=True)
sl_era5=''.join(sl_era5)
sl_era5 = netCDF4.Dataset(sl_era5)

###############################
###set up era5 time for loop###
###############################

era_time=pl_era5.variables['time'][:]
era_time = "2002-12-15 00:00:00"     ##refer to era5 data
era_time = datetime.strptime(era_time, "%Y-%m-%d %H:%M:%S")
era_time = era_time.timestamp()

start_time = "2002-12-15 00:00:00"     ##wrf simulation start date + the time you want to add cmip6 anomaly
start_time = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
start_time = start_time.timestamp()


end_time = "2002-12-23 23:00:00"     ##wrf simulation end date
end_time = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
end_time = end_time.timestamp()

i_start=int((start_time-era_time)/60/60)
i_end=int((end_time-era_time)/60/60)


##########################
###load cmip6 data     ###
##########################

cmip6_present_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/cmip6/cmip_25km/so'
os.chdir(cmip6_present_path)
cmip6_present_file= glob.glob("present.nc", recursive=True)
cmip6_present_file=''.join(cmip6_present_file)
cmip6_present_file = netCDF4.Dataset(cmip6_present_file)

cmip6_future_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/cmip6/cmip_25km/so'
os.chdir(cmip6_future_path)
cmip6_future_file= glob.glob("future.nc", recursive=True)
cmip6_future_file=''.join(cmip6_future_file)
cmip6_future_file = netCDF4.Dataset(cmip6_future_file)
#cmip6_future_file.variables.keys()

cmip6_present_ta=np.squeeze(cmip6_present_file.variables['ta'])
cmip6_present_hur=np.squeeze(cmip6_present_file.variables['hur'])
cmip6_present_ts=np.squeeze(cmip6_present_file.variables['ts'])

cmip6_future_ta=np.squeeze(cmip6_future_file.variables['ta'])
cmip6_future_hur=np.squeeze(cmip6_future_file.variables['hur'])
cmip6_future_ts=np.squeeze(cmip6_future_file.variables['ts'])

cmip6_lon=cmip6_future_file.variables['lon'][:]-180
cmip6_lat=cmip6_future_file.variables['lat'][::-1]

cmip6_Anom_ta=(cmip6_future_ta-cmip6_present_ta)[::-1,:,:]
cmip6_Anom_hur=(cmip6_future_ta-cmip6_present_hur)[::-1,:,:]
cmip6_Anom_ts=cmip6_future_ts-cmip6_present_ts

######################################
###interpolate cmip6 into era5 grid###
######################################

cmip6_X, cmip6_Y = np.meshgrid(cmip6_lon, cmip6_lat)
era_X, era_Y = np.meshgrid(era_lon, era_lat)
small_Anom_ts = griddata((cmip6_X.flatten(), cmip6_Y.flatten()), cmip6_Anom_ts.flatten(), (era_X, era_Y), method='linear')

small_Anom_ta=np.empty((cmip6_Anom_ta.shape[0], era_lat.shape[0], era_lon.shape[0]))
small_Anom_hur=np.empty((cmip6_Anom_ta.shape[0],  era_lat.shape[0], era_lon.shape[0]))

for zz in range(0,len(cmip6_Anom_hur)):
    small_Anom_ta[zz,:,:]= griddata((cmip6_X.flatten(), cmip6_Y.flatten()), cmip6_Anom_ta[zz,:,:].flatten(), (era_X, era_Y), method='linear')
    small_Anom_hur[zz,:,:]= griddata((cmip6_X.flatten(), cmip6_Y.flatten()), cmip6_Anom_hur[zz,:,:].flatten(), (era_X, era_Y), method='linear')
    
small_cmip6_ta = np.nan_to_num(small_Anom_ta, nan=0)
small_cmip6_hur = np.nan_to_num(small_Anom_hur, nan=0)
small_cmip6_ts = np.nan_to_num(small_Anom_ts, nan=0)



for i in range(i_start,i_end):
    
    era_lon=pl_era5.variables['lon'][:]
    era_lat=pl_era5.variables['lat'][:]
    
    era_2d_u10m=sl_era5.variables['10u'][i]
    era_2d_v10m=sl_era5.variables['10v'][i]
    era_2d_dpt2m=sl_era5.variables['2d'][i]
    era_2d_tp2m=sl_era5.variables['2t'][i]
    era_2d_lsm=sl_era5.variables['lsm'][i]
    era_2d_msl=sl_era5.variables['msl'][i]
    era_2d_ci=sl_era5.variables['ci'][i]
    era_2d_sst=sl_era5.variables['sst'][i]
    era_2d_skt=sl_era5.variables['skt'][i] #+ small_cmip6_ts #switch with cmip6
    era_2d_sd=sl_era5.variables['sd'][i]
    era_2d_sp=sl_era5.variables['sp'][i]
    
    era_3d_t=pl_era5.variables['t'][i]  #+ small_cmip6_ta #switch with cmip6
    era_3d_u=pl_era5.variables['u'][i]
    era_3d_v=pl_era5.variables['v'][i]
    era_3d_z=pl_era5.variables['z'][i]
    era_3d_q=pl_era5.variables['q'][i] 
    era_3d_r=pl_era5.variables['r'][i] #+ small_cmip6_hur #switch with cmip6
    #  rh = 1.E2 * (p*q/(q*(1.-eps) + eps))/(svp1*exp(svp2*(t-svpt0)/(T-svp3)))
    #era_3d_r= 1.E2 * (era_2d_sp*era_3d_q/(era_3d_q*(1.-0.622) + 0.622))/(611.1*np.exp(17.67*(era_3d_t-273.15)/(era_2d_t-29.65)))

    era_soil_st1=sl_era5.variables['stl1'][i]
    era_soil_st2=sl_era5.variables['stl2'][i]
    era_soil_st3=sl_era5.variables['stl3'][i]
    era_soil_st4=sl_era5.variables['stl4'][i]
    era_soil_swvl1=sl_era5.variables['swvl1'][i]
    era_soil_swvl2=sl_era5.variables['swvl2'][i]
    era_soil_swvl3=sl_era5.variables['swvl3'][i]
    era_soil_swvl4=sl_era5.variables['swvl4'][i]
    
    winter_geo, total_fields = creat_FILI_for_WPS(era_lat,era_lon,
                       era_2d_u10m,era_2d_v10m,era_2d_tp2m,era_2d_lsm,era_2d_msl,era_2d_ci,era_2d_sst,era_2d_skt,era_2d_sd,era_2d_sp,                
                       era_3d_t,era_3d_u,era_3d_v,era_3d_r,era_3d_z,
                       era_soil_st1,era_soil_st2,era_soil_st3,era_soil_st4,era_soil_swvl1,era_soil_swvl2,era_soil_swvl3,era_soil_swvl4)
    
    ###output setup###
    now_time=era_time+i*60*60
    now_time = datetime.fromtimestamp(now_time)
    now_time = now_time.strftime("%Y-%m-%d %H:%M:%S")  
    now_time[0:10]+'_'+now_time[11:13]
    
    path_out = '/home/x_tilin/snic2021-23-400/users/x_tilin/run/project3/pl_Erik2002_ungrib/wps/'
    pyw.cinter('FILE',now_time[0:10]+'_'+now_time[11:13],winter_geo,total_fields,path_out)

    

    
    
    

# nc_file = netCDF4.Dataset('/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project3/pl_mcao1/ERA5-20141215-20141231-pl.nc', 'r+')
# variable_name = 't'
# variable = nc_file.variables[variable_name]
# new_values = variable[:] * 2
# variable[:] = new_values
# nc_file.close()


# levl=range(240,290,5)
# bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
#         rsphere=(6378137.00,6356752.3142),\
#         resolution='l',projection='merc',\
#         lat_0=40.,lon_0=40.,lat_ts=30.)    
    
# fig = plt.figure(figsize=(15,7)) 
# # lon_grid,lat_grid=np.meshgrid(cmip6_lon,cmip6_lat)
# lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

# xx,yy=bm(lon_grid,lat_grid)            
# wspd_contours1 = bm.contourf(xx, yy, era_2d_skt, levels=levl,
#                         extend='both')
# # wspd_contours1 = bm.contourf(xx, yy, cimp6_Anom_ts, levels=levl, 
# #                         extend='both')
# plt.colorbar(wspd_contours1, fraction=0.046,orientation="horizontal", pad=.05)

# bm.drawcoastlines(linewidth=0.25) # Add the geographic boundaries
# bm.drawstates(linewidth=0.25)
# bm.drawcountries(linewidth=0.25)
# bm.drawparallels(np.arange(-80., 80., 20.),labels=[1,0,0,0], fontsize=13)
# bm.drawmeridians(np.arange(-180., 180., 40.),labels=[0,0,0,1],  fontsize=13)
# bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),                                 
#                                   lake_color=np.array([0.59375 ,
#                                                       0.71484375,
#                                                       0.8828125 ]))  
# # plt.title(str( _polarlow_time))
# # plt.savefig(save_dir+'PL_'+str(idx)+"_track_slp{}.png".format(wrf_idx),dpi=600,bbox_inches = 'tight') 
# plt.show()
                  







