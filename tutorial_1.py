#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 18:59:30 2024

@author: vinita
"""


#------------------------------------------------------------------------------
# Computing trajectories example
#------------------------------------------------------------------------------
 
#import lagrantraj.trajectories as traj

import lagrantraj.trajectories as traj

#------------------------------------------------------------------------------
# Reading input netcdf file
#------------------------------------------------------------------------------
""" def read_data(File_name,Root_input,list_var,list_var_advec,lat='latitude',lon='longitude',pres='isobaricInhPa')
    list_var are mandatory (u,v,w) wind componenets list_var_advect , additional variables to compute 
    thier values along the trajetories 
    substitude the coordinate variables (e.g lat,lon,press)"""
    
list_var = ['u','v','w']
list_var_advect = ['pv','pt']
filename = 'TC1279_cont_dec_merged.nc'
root_input = '/home/vinita/VINITA/ECMWF/'
root_output = '/home/vinita/VINITA/ECMWF/'
LON_nc,LAT_nc,P_nc,data = traj.read_data(filename,root_input,list_var_advect, list_var,lat='latitude',lon='longitude',pres='isobaricInhPa')


lat_seeds = ([74., 74., 73., 73., 73., 73., 73., 73., 72., 72., 72., 72., 72.,
             72., 72., 72., 71., 71., 71.])
lon_seeds = ([ -55.,  -45.,  -60.,  -56.,  -52.,  -48.,  -44.,  -40.,  -67.,
               -63.,  -59.,  -55.,  -51.,  -47.,  -43.,  -39.,  -71.,  -67.,
               -63.])
pres_seeds = [30000]*19
initial_time_step = 86

""" def compute_trajectories(x0,y0,z0,initial_time_index,
                         LON_nc,LAT_nc,P_nc,data,
                         list_var,list_var_advec,
                         trajectories_duration=None,
                         dt_data=6.,dt_traj=0.5,
                         niter=4,BACKWARD=True):
     Compute Lagrangian trajectories - traj_duration is in hours 
                                     - dt_data is input data temporal resolution in hours 
                                     - dt_traj is output trajectories temporal resolution in hours """
TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj=traj.compute_trajectories(lon_seeds,lat_seeds,pres_seeds,initial_time_step,
                              LON_nc,LAT_nc,P_nc,data,
                              list_var_advect, list_var,
                              trajectories_duration=72,
                              dt_data=3.,dt_traj=0.5,
                              niter=4,BACKWARD=True)



#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------
P_traj=P_traj/100 # in hpa
LON_traj=LON_traj
LAT_traj=LAT_traj
n_seeds =    LAT_traj.shape[0]
color=P_traj


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import numpy as np

fig = plt.figure(figsize=(15,12))
ax=plt.subplot(projection=ccrs.NorthPolarStereo())
ax.scatter(LON_traj[:,0],LAT_traj[:,0], c=P_traj[:,0], edgecolors='black',
           cmap='Greens',transform=ccrs.PlateCarree())

#a= plt.contour(geopt.longitude,geopt.latitude[:],(geopt.z[0,26,:,:,1]/100),colors='black',transform=ccrs.PlateCarree())
#plt.clabel(a, inline=1, fontsize=10)
extent = 2500000
ax.set_extent((-extent,extent,-extent,extent),crs=ccrs.NorthPolarStereo())
plt.title(' Trajectories map (from 300 hpa pressure level) ', size=26)
ax.set_extent([-180, 180,30, 90], ccrs.PlateCarree())
ax.coastlines(linewidth=0.2)

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



#------------------------------------------------------------------------------
# saving Data in NetCDF format
#------------------------------------------------------------------------------
"""save_output_data(Root_output,initial_time_index,
                     list_var,list_var_advec,
                     TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj)"""
    
    
traj.save_output_data(root_output,initial_time_step,
                     list_var_advect,list_var,
                     TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj)    
    
    
    
    
