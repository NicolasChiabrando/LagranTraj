#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 18:59:30 2024

@author: vinita
"""

# =============================================================================
#  # import traj_main module 
# =============================================================================
 
#import lagrantraj.trajectories as traj

import lagrantraj.trajectories as traj
# =============================================================================
# # Read the nc file using the read_data
# =============================================================================
LON_nc,LAT_nc,P_nc,data = traj.read_data_ecmwf('TC1279_cont_dec_merged.nc','/home/vinita/VINITA/ECMWF/',['pv','pt'],['u','v','w'])


lon_seeds = ([-80., -79., -78., -77., -76., -75.])
lat_seeds = ([35., 34., 33., 32., 31., 30.])
pres_seeds = [30000]*6
initial_time_step = 83


TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj=traj.compute_trajectories(lon_seeds,lat_seeds,pres_seeds,initial_time_step,
                              LON_nc,LAT_nc,P_nc,data,
                              ['pv'],['u','v','w'],
                              trajectories_duration=36,
                              dt_data=3.,dt_traj=0.5,
                              niter=4,BACKWARD=False)


# =============================================================================
# data_traj = data_traj.isel(n_seeds=upper_indices)
# =============================================================================
P_traj=P_traj
LON_traj=LON_traj
LAT_traj=LAT_traj
n_seeds =    LAT_traj.shape[0]
color=P_traj


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import numpy as np

fig = plt.figure(figsize=(18,12))
ax=plt.subplot(projection=ccrs.NorthPolarStereo())
ax.scatter(LON_traj[:,0],LAT_traj[:,0], c=P_traj[:,0], edgecolors='black',
           cmap='Greens',transform=ccrs.PlateCarree())

#a= plt.contour(geopt.longitude,geopt.latitude[:],(geopt.z[0,26,:,:,1]/100),colors='black',transform=ccrs.PlateCarree())
#plt.clabel(a, inline=1, fontsize=10)
extent = 2500000
ax.set_extent((-extent,extent,-extent,extent),crs=ccrs.NorthPolarStereo())
plt.title(' Trajectories from 300 hpa pressure level', size=26)
ax.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
ax.coastlines()

for i_traj in range(n_seeds):
    points = np.array([LON_traj[i_traj,:], LAT_traj[i_traj,:]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(np.nanmin(color),np.nanmax(color))
    lc = LineCollection(segments, cmap='jet', norm=norm,transform=ccrs.Geodetic())
    lc.set_array(color[i_traj,:])
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
plt.xlim([(np.nanmin(LON_traj))-0.5,(np.nanmax(LON_traj))+0.5])
#print(np.nanmin(LON_traj)
plt.ylim([(np.nanmin(LAT_traj))-0.5,(np.nanmax(LAT_traj))+0.5])
  
#cbar_ax = fig.add_axes([0.92, 0.125, 0.02, 0.755])
colo = fig.colorbar(lc,shrink=0.9)
colo.ax.tick_params(labelsize=23)
colo.set_label(label='Pressure [Hpa]', size=23)
ax.set_extent([-180,180,20,90], ccrs.PlateCarree())
ax.coastlines()
plt.show()
print('ok')

