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


def specfic2relative_humidity(t,q): ##refer to WRF source code wps/ungrib/src/rrpr.f90
    svp1=611.2
    svp2=17.67
    svp3=29.65
    svpt0=273.15
    eps=0.622
    plevs=np.array([100, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000,
                    50000, 60000, 70000, 85000, 92500, 100000])
    p = np.empty_like(era_3d_t)
    p[:len(plevs), :, :] = np.reshape(plevs, (-1, 1, 1))
    q[q<1.E-10]=1.E-10#q=max(1.E-10,q)
    rh = 1.E2 * (p*q/(q*(1.-eps) + eps))/(svp1*np.exp(svp2*(t-svpt0)/(t-svp3)))
    return rh


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
    winter_sst = pyw.V2d('SST',sst,'Sea-Surface Temperature','K','200100')  ##have not downloaded cmip6 sst
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


end_time = "2002-12-21 05:00:00"     ##wrf simulation end date
end_time = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
end_time = end_time.timestamp()

i_start=int((start_time-era_time)/60/60)
i_end=int((end_time-era_time)/60/60)


##########################
###load cmip6 data     ###
##########################

cmip6_present_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/cmip6/ssp245/ACCESS-CM2'
os.chdir(cmip6_present_path)

cmip6_history_ci= glob.glob('siconc*'+'*historical*', recursive=True)
cmip6_history_ci=''.join(cmip6_history_ci)
cmip6_history_ci = netCDF4.Dataset(cmip6_history_ci)
ci_lon=cmip6_history_ci.variables['longitude'][:]
ci_lat=cmip6_history_ci.variables['latitude'][:]
cmip6_present_ci=np.nanmean(np.squeeze(cmip6_history_ci.variables['siconc']),axis=0)

cmip6_history_sst= glob.glob('tos*'+'*historical*', recursive=True)
cmip6_history_sst=''.join(cmip6_history_sst)
cmip6_history_sst = netCDF4.Dataset(cmip6_history_sst)
sst_lon=cmip6_history_sst.variables['longitude'][:]
sst_lat=cmip6_history_sst.variables['latitude'][:]
cmip6_present_sst=np.nanmean(np.squeeze(cmip6_history_sst.variables['tos']),axis=0)

cmip6_history_ts= glob.glob('ts*'+'*historical*', recursive=True)
cmip6_history_ts=''.join(cmip6_history_ts)
cmip6_history_ts = netCDF4.Dataset(cmip6_history_ts)
cmip6_present_ts=np.nanmean(np.squeeze(cmip6_history_ts.variables['ts']),axis=0)[::-1,:]

cmip6_history_ta= glob.glob('ta*'+'*historical*', recursive=True)
cmip6_history_ta=''.join(cmip6_history_ta)
cmip6_history_ta = netCDF4.Dataset(cmip6_history_ta)
cmip6_present_ta=np.nanmean(np.squeeze(cmip6_history_ta.variables['ta']),axis=0)[::-1,::-1,:]

# cmip6_history_hus= glob.glob('hus*'+'*historical*', recursive=True)
# cmip6_history_hus=''.join(cmip6_history_hus)
# cmip6_history_hus = netCDF4.Dataset(cmip6_history_hus)
# cmip6_present_hus=np.nanmean(np.squeeze(cmip6_history_hus.variables['hus']),axis=0)[::-1,::-1,:]




cmip6_future_path='/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/cmip6/ssp245/ACCESS-CM2'
os.chdir(cmip6_future_path)

cmip6_ssp370_ci= glob.glob('siconc*'+'*ssp245*', recursive=True)
cmip6_ssp370_ci=''.join(cmip6_ssp370_ci)
cmip6_ssp370_ci = netCDF4.Dataset(cmip6_ssp370_ci)
ci_lon=cmip6_ssp370_ci.variables['longitude'][:]
ci_lat=cmip6_ssp370_ci.variables['latitude'][:]
cmip6_future_ci=np.nanmean(np.squeeze(cmip6_ssp370_ci.variables['siconc']),axis=0)

cmip6_ssp370_sst= glob.glob('tos*'+'*ssp245*', recursive=True)
cmip6_ssp370_sst=''.join(cmip6_ssp370_sst)
cmip6_ssp370_sst = netCDF4.Dataset(cmip6_ssp370_sst)
sst_lon=cmip6_ssp370_sst.variables['longitude'][:]
sst_lat=cmip6_ssp370_sst.variables['latitude'][:]
cmip6_future_sst=np.nanmean(np.squeeze(cmip6_ssp370_sst.variables['tos']),axis=0)

cmip6_ssp370_ts= glob.glob('ts*'+'*ssp245*', recursive=True)
cmip6_ssp370_ts=''.join(cmip6_ssp370_ts)
cmip6_ssp370_ts = netCDF4.Dataset(cmip6_ssp370_ts)
cmip6_future_ts=np.nanmean(np.squeeze(cmip6_ssp370_ts.variables['ts']),axis=0)[::-1,:]

cmip6_ssp370_ta= glob.glob('ta*'+'*ssp245*', recursive=True)
cmip6_ssp370_ta=''.join(cmip6_ssp370_ta)
atm_lon=cmip6_ssp370_ts.variables['lon'][:]
atm_lat=cmip6_ssp370_ts.variables['lat'][:][::-1]
cmip6_ssp370_ta = netCDF4.Dataset(cmip6_ssp370_ta)
cmip6_future_ta=np.nanmean(np.squeeze(cmip6_ssp370_ta.variables['ta']),axis=0)[::-1,::-1,:]

# cmip6_ssp370_hus= glob.glob('hus*'+'*ssp370*', recursive=True)
# cmip6_ssp370_hus=''.join(cmip6_ssp370_hus)
# cmip6_ssp370_hus = netCDF4.Dataset(cmip6_ssp370_hus)
# cmip6_future_hus=np.nanmean(np.squeeze(cmip6_ssp370_hus.variables['hus']),axis=0)[::-1,::-1,:]



cmip6_Anom_ta=cmip6_future_ta-cmip6_present_ta
# cmip6_Anom_hus=cmip6_future_hus-cmip6_present_hus
cmip6_Anom_ts=cmip6_future_ts-cmip6_present_ts
cmip6_Anom_sst=cmip6_future_sst-cmip6_present_sst
cmip6_Anom_ci=(cmip6_future_ci-cmip6_present_ci)/100

######################################
###interpolate cmip6 into era5 grid###
######################################
era_X, era_Y = np.meshgrid(era_lon, era_lat)

cmip6_X, cmip6_Y = np.meshgrid(atm_lon, atm_lat)
era_X, era_Y = np.meshgrid(era_lon, era_lat)

cmip6_2_era5_tsAnom = griddata(np.hstack((cmip6_X.flatten()[:,None], cmip6_Y.flatten()[:,None])), cmip6_Anom_ts.flatten(), (era_X, era_Y), method='linear')
cmip6_2_era5_taAnom=np.empty((cmip6_Anom_ta.shape[0], era_lat.shape[0], era_lon.shape[0]))
# cmip6_2_era5_husAnom=np.empty((cmip6_Anom_ta.shape[0],  era_lat.shape[0], era_lon.shape[0]))

for zz in range(0,len(cmip6_Anom_ta)):
    cmip6_2_era5_taAnom[zz,:,:]= griddata(np.hstack((cmip6_X.flatten()[:,None], cmip6_Y.flatten()[:,None])), cmip6_Anom_ta[zz,:,:].flatten(), (era_X, era_Y), method='linear')
    # cmip6_2_era5_husAnom[zz,:,:]= griddata(np.hstack((cmip6_X.flatten()[:,None], cmip6_Y.flatten()[:,None])), cmip6_Anom_hus[zz,:,:].flatten(), (era_X, era_Y), method='linear')
    
anom_cmip6_ta = np.nan_to_num(cmip6_2_era5_taAnom, nan=0)
# anom_cmip6_hus = np.nan_to_num(cmip6_2_era5_husAnom, nan=0)
anom_cmip6_ts = np.nan_to_num(cmip6_2_era5_tsAnom, nan=0)


cmip6_2_era5_sstAnom = griddata(np.hstack((sst_lon.flatten()[:,None], sst_lat.flatten()[:,None])), cmip6_Anom_sst.flatten(), (era_X, era_Y), method='linear')
cmip6_2_era5_ciAnom = griddata(np.hstack((ci_lon.flatten()[:,None], ci_lat.flatten()[:,None])), cmip6_Anom_ci.flatten(), (era_X, era_Y), method='linear')
anom_cmip6_sst= np.nan_to_num(cmip6_2_era5_sstAnom, nan=0)
anom_cmip6_ci = np.nan_to_num(cmip6_2_era5_ciAnom, nan=0)


# cmip6_ci_X, cmip6_ci_Y = np.meshgrid(ci_lon, ci_lat)
# cmip6_sst_X, cmip6_sst_Y = np.meshgrid(sst_lon, sst_lat)
# cmip6_2_era5_sstAnom = griddata(np.hstack((cmip6_sst_X.flatten()[:,None], cmip6_sst_Y.flatten()[:,None])), cmip6_Anom_sst.flatten(), (era_X, era_Y), method='linear')
# cmip6_2_era5_ciAnom = griddata(np.hstack((cmip6_ci_X.flatten()[:,None], cmip6_ci_Y.flatten()[:,None])), cmip6_Anom_ci.flatten(), (era_X, era_Y), method='linear')
# anom_cmip6_sst= np.nan_to_num(cmip6_2_era5_sstAnom, nan=0)
# anom_cmip6_ci = np.nan_to_num(cmip6_2_era5_ciAnom, nan=0)
#####################
### for NorESM2-MM ##
##
# cmip6_2_era5_sstAnom = griddata(np.hstack((sst_lon.flatten()[:,None], sst_lat.flatten()[:,None])), cmip6_Anom_sst.flatten(), (era_X, era_Y), method='linear')
# nan_indices = np.isnan(ci_lon)
# ci_lon_without_nan = ci_lon[~nan_indices]
# nan_indices = np.isnan(ci_lat)
# ci_lat_without_nan = ci_lat[~nan_indices]
# cmip6_Anom_ci_without_nan = cmip6_Anom_ci[~nan_indices]
# cmip6_2_era5_ciAnom = griddata(np.hstack((ci_lon_without_nan.flatten()[:,None], ci_lat_without_nan.flatten()[:,None])), cmip6_Anom_ci_without_nan, (era_X, era_Y), method='linear')
# anom_cmip6_sst= np.nan_to_num(cmip6_2_era5_sstAnom, nan=0)
# anom_cmip6_ci = np.nan_to_num(cmip6_2_era5_ciAnom, nan=0)

# cmip6_present_ci_without_nan = cmip6_present_ci[~nan_indices]
# cmip6_2_era5_ci = griddata(np.hstack((ci_lon_without_nan.flatten()[:,None], ci_lat_without_nan.flatten()[:,None])), cmip6_present_ci_without_nan, (era_X, era_Y), method='linear')

# ##test cmip6 atm,oce,ice
cmip6_2_era5_sst = griddata(np.hstack((sst_lon.flatten()[:,None], sst_lat.flatten()[:,None])), cmip6_present_sst.flatten(), (era_X, era_Y), method='linear')
cmip6_2_era5_ci = griddata(np.hstack((ci_lon.flatten()[:,None], ci_lat.flatten()[:,None])), cmip6_present_ci.flatten(), (era_X, era_Y), method='linear')
cmip6_2_era5_ta = griddata(np.hstack((cmip6_X.flatten()[:,None], cmip6_Y.flatten()[:,None])), cmip6_present_ta[18,:,:].flatten(), (era_X, era_Y), method='linear')
cmip6_2_era5_ts = griddata(np.hstack((cmip6_X.flatten()[:,None], cmip6_Y.flatten()[:,None])), cmip6_present_ts.flatten(), (era_X, era_Y), method='linear')


levl=np.arange(0,1,0.2)
bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc',\
        lat_0=40.,lon_0=40.,lat_ts=30.)    
    
fig = plt.figure(figsize=(15,7)) 
lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

xx,yy=bm(lon_grid,lat_grid)            
wspd_contours1 = bm.contourf(xx, yy, cmip6_2_era5_ci, levels=levl,
                        extend='both')
# wspd_contours1 = bm.contourf(xx, yy, cimp6_Anom_ts, levels=levl, 
#                         extend='both')
plt.colorbar(wspd_contours1, fraction=0.046,orientation="horizontal", pad=.05)

bm.drawcoastlines(linewidth=0.25) # Add the geographic boundaries
bm.drawstates(linewidth=0.25)
bm.drawcountries(linewidth=0.25)
bm.drawparallels(np.arange(-80., 80., 20.),labels=[1,0,0,0], fontsize=13)
bm.drawmeridians(np.arange(-180., 180., 40.),labels=[0,0,0,1],  fontsize=13)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),                                 
                                  lake_color=np.array([0.59375 ,
                                                      0.71484375,
                                                      0.8828125 ]))  
# plt.title(str( _polarlow_time))
# plt.savefig(save_dir+'PL_'+str(idx)+"_track_slp{}.png".format(wrf_idx),dpi=600,bbox_inches = 'tight') 
plt.show()
            
levl=np.arange(260,280,2)
bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc',\
        lat_0=40.,lon_0=40.,lat_ts=30.)    
    
fig = plt.figure(figsize=(15,7)) 
lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

xx,yy=bm(lon_grid,lat_grid)            
wspd_contours1 = bm.contourf(xx, yy, cmip6_2_era5_sst+273.15, levels=levl,
                        extend='both')
# wspd_contours1 = bm.contourf(xx, yy, cimp6_Anom_ts, levels=levl, 
#                         extend='both')
plt.colorbar(wspd_contours1, fraction=0.046,orientation="horizontal", pad=.05)

bm.drawcoastlines(linewidth=0.25) # Add the geographic boundaries
bm.drawstates(linewidth=0.25)
bm.drawcountries(linewidth=0.25)
bm.drawparallels(np.arange(-80., 80., 20.),labels=[1,0,0,0], fontsize=13)
bm.drawmeridians(np.arange(-180., 180., 40.),labels=[0,0,0,1],  fontsize=13)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),                                 
                                  lake_color=np.array([0.59375 ,
                                                      0.71484375,
                                                      0.8828125 ]))  
# plt.title(str( _polarlow_time))
# plt.savefig(save_dir+'PL_'+str(idx)+"_track_slp{}.png".format(wrf_idx),dpi=600,bbox_inches = 'tight') 
plt.show()

levl=np.arange(260,280,2)
bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc',\
        lat_0=40.,lon_0=40.,lat_ts=30.)    
    
fig = plt.figure(figsize=(15,7)) 
lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

xx,yy=bm(lon_grid,lat_grid)            
wspd_contours1 = bm.contourf(xx, yy, cmip6_2_era5_ts, levels=levl,
                        extend='both')
# wspd_contours1 = bm.contourf(xx, yy, cimp6_Anom_ts, levels=levl, 
#                         extend='both')
plt.colorbar(wspd_contours1, fraction=0.046,orientation="horizontal", pad=.05)

bm.drawcoastlines(linewidth=0.25) # Add the geographic boundaries
bm.drawstates(linewidth=0.25)
bm.drawcountries(linewidth=0.25)
bm.drawparallels(np.arange(-80., 80., 20.),labels=[1,0,0,0], fontsize=13)
bm.drawmeridians(np.arange(-180., 180., 40.),labels=[0,0,0,1],  fontsize=13)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),                                 
                                  lake_color=np.array([0.59375 ,
                                                      0.71484375,
                                                      0.8828125 ]))  
# plt.title(str( _polarlow_time))
# plt.savefig(save_dir+'PL_'+str(idx)+"_track_slp{}.png".format(wrf_idx),dpi=600,bbox_inches = 'tight') 
plt.show()

levl=np.arange(260,280,2)
bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc',\
        lat_0=40.,lon_0=40.,lat_ts=30.)    
    
fig = plt.figure(figsize=(15,7)) 
lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

xx,yy=bm(lon_grid,lat_grid)            
wspd_contours1 = bm.contourf(xx, yy, cmip6_2_era5_ta, levels=levl,
                        extend='both')
# wspd_contours1 = bm.contourf(xx, yy, cimp6_Anom_ts, levels=levl, 
#                         extend='both')
plt.colorbar(wspd_contours1, fraction=0.046,orientation="horizontal", pad=.05)

bm.drawcoastlines(linewidth=0.25) # Add the geographic boundaries
bm.drawstates(linewidth=0.25)
bm.drawcountries(linewidth=0.25)
bm.drawparallels(np.arange(-80., 80., 20.),labels=[1,0,0,0], fontsize=13)
bm.drawmeridians(np.arange(-180., 180., 40.),labels=[0,0,0,1],  fontsize=13)
bm.fillcontinents(color=np.array([ 0.9375 , 0.9375 , 0.859375]),                                 
                                  lake_color=np.array([0.59375 ,
                                                      0.71484375,
                                                      0.8828125 ]))  
# plt.title(str( _polarlow_time))
# plt.savefig(save_dir+'PL_'+str(idx)+"_track_slp{}.png".format(wrf_idx),dpi=600,bbox_inches = 'tight') 
plt.show()

for i in range(i_start,i_end):
    
    era_lon=pl_era5.variables['lon'][:]
    era_lat=pl_era5.variables['lat'][:]
    
    era_2d_u10m=sl_era5.variables['10u'][i]
    era_2d_v10m=sl_era5.variables['10v'][i]
    era_2d_dpt2m=sl_era5.variables['2d'][i]
    era_2d_tp2m=sl_era5.variables['2t'][i]
    era_2d_lsm=sl_era5.variables['lsm'][i]
    era_2d_msl=sl_era5.variables['msl'][i]
    era_2d_ci=sl_era5.variables['ci'][i] +anom_cmip6_ci
    era_2d_ci[era_2d_ci < 0] = -1e+30
    era_2d_sst=sl_era5.variables['sst'][i] + anom_cmip6_sst
    era_2d_skt=sl_era5.variables['skt'][i] + anom_cmip6_ts #switch with cmip6
    era_2d_sd=sl_era5.variables['sd'][i]
    era_2d_sp=sl_era5.variables['sp'][i]
    
    era_3d_t=pl_era5.variables['t'][i]  + anom_cmip6_ta #switch with cmip6
    era_3d_u=pl_era5.variables['u'][i]
    era_3d_v=pl_era5.variables['v'][i]
    era_3d_z=pl_era5.variables['z'][i]
    # era_3d_q=pl_era5.variables['q'][i] + anom_cmip6_hus
    # era_3d_r=specfic2relative_humidity(era_3d_t,era_3d_q)
    era_3d_q=pl_era5.variables['q'][i]
    era_3d_r=pl_era5.variables['r'][i] #+ small_cmip6_hur #switch with cmip6


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
    
    path_out = '/home/x_tilin/snic2021-23-400/users/x_tilin/run/project3/pl_Erik2002/ssp245/pl_Erik2002_ACCESS-CM2/wps/'
    pyw.cinter('FILE',now_time[0:10]+'_'+now_time[11:13],winter_geo,total_fields,path_out)

    

    
    
    

# nc_file = netCDF4.Dataset('/home/x_tilin/snic2021-23-400/users/x_tilin/input_data/Boundary/WRF/project3/pl_mcao1/ERA5-20141215-20141231-pl.nc', 'r+')
# variable_name = 't'
# variable = nc_file.variables[variable_name]
# new_values = variable[:] * 2
# variable[:] = new_values
# nc_file.close()
 

#############
##test era5##
#############


# levl=np.arange(0,1,0.2)
# bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
#         rsphere=(6378137.00,6356752.3142),\
#         resolution='l',projection='merc',\
#         lat_0=40.,lon_0=40.,lat_ts=30.)    
    
# fig = plt.figure(figsize=(15,7)) 
# lon_grid,lat_grid=np.meshgrid(era_lon,era_lat)  

# xx,yy=bm(lon_grid,lat_grid)            
# wspd_contours1 = bm.contourf(xx, yy, cmip6_2_era5_ci, levels=levl,
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
                  



##############
##test cmip6##
##############

# levl=np.arange(260,280,2)
# bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
#         rsphere=(6378137.00,6356752.3142),\
#         resolution='l',projection='merc',\
#         lat_0=40.,lon_0=40.,lat_ts=30.)    
    
# fig = plt.figure(figsize=(15,7)) 
# lon_grid,lat_grid=np.meshgrid(atm_lon,atm_lat)


# xx,yy=bm(lon_grid,lat_grid)            
# wspd_contours1 = bm.contourf(xx, yy, cmip6_present_ts, levels=levl,
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
                  

##############
##test ocean##
##############

# levl=np.arange(0,10.1,0.1)
# bm = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=85.,\
#         rsphere=(6378137.00,6356752.3142),\
#         resolution='l',projection='merc',\
#         lat_0=40.,lon_0=40.,lat_ts=30.)    

# # bm = Basemap(projection='spstere', boundinglat=40, lon_0=180, resolution='l')   
# fig = plt.figure(figsize=(15,7)) 

# lon_grid,lat_grid=np.meshgrid(sst_lon,sst_lat)
# xx,yy=bm(lon_grid,lat_grid)             
# wspd_contours1 = bm.contourf(xx, yy, cmip6_present_sst, levels=levl,
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


# plt.contourf(era_3d_q[15])
# plt.colorbar()



