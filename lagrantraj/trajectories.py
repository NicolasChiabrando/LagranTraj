
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:38:29 2024

@author: vinita
"""


#------------------------------------------------------------------------------
# Packages:
#------------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import interpn, interp1d
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import sys
import time
#------------------------------------------------------------------------------
# Constants:
#------------------------------------------------------------------------------    
Ra=6370000.0 # m : Radius of the Earth
cp=1005.0 # J/kg/K : Specific heat of dry air
Lv=2500000 # J/kg : Latent heat of vaporization at 0◦C
ZR=287.05 #J/K/mol: Gas constant for dry air
ZKP=ZR/cp  # Ra/Cp (pour calcul de la température potentielle)
g=9.81 # N/kg: Acceleration due to gravity at sea level



def read_data(File_name,Root_input,list_var,list_var_advec,lat='latitude',lon='longitude',pres='isobaricInhPa', hPa = True):
    """ Read data """
    print('Read data: ', end='', file = sys.stderr)
    data={}
    
    for i_var in list_var_advec+list_var:
        data[i_var]=   Dataset(Root_input+File_name,format='NETCDF4')
            
            
    LON_nc=data[list_var_advec[0]][lon][:]
    LAT_nc=data[list_var_advec[0]][lat][:]        
    P_nc=data[list_var_advec[0]][pres][:]
    if hPa:
    	P_nc*=100
    
    return LON_nc,LAT_nc,P_nc,data

#------------------------------------------------------------------------------
# Read data multiple files where the filename is stored in a dictionary filename ={}
#------------------------------------------------------------------------------ 
def read_data_multiple(File_name,Root_input,list_var,list_var_advec, hPa = True):
    """ Read data """
    print('Read data: ', end='', file = sys.stderr)
    if isinstance(File_name,dict):
        data={}
        for i_var in list_var_advec+list_var:
            data[i_var]=   Dataset(Root_input+File_name[i_var],format='NETCDF4')
    else:
        data={}
        for i_var in list_var_advec+list_var:
            data[i_var]=   Dataset(Root_input+File_name,format='NETCDF4')
            
    LON_nc=data[list_var_advec[0]].variables['lon'][:]
    LAT_nc=data[list_var_advec[0]].variables['lat'][:]        
    P_nc=data[list_var_advec[0]].variables['level'][:]
    if hPa:
    	P_nc*=100

    return LON_nc,LAT_nc,P_nc,data




#------------------------------------------------------------------------------    
# Trajectories definition
#------------------------------------------------------------------------------        
def compute_trajectories(x0,y0,z0,initial_time_index,
                         LON_nc,LAT_nc,P_nc,data,
                         list_var,list_var_advec,
                         trajectories_duration=None,
                         dt_data=6.,dt_traj=0.5,
                         niter=4,BACKWARD=True):
    """ Compute Lagrangian trajectories - traj_duration is in hours"""
    
    print('Computing backward trajectories...: ', file = sys.stderr)

    # Defining number of seeds
    Number_Seeds=len(x0)

    # Defining trajectories direction
    direction=1
    if BACKWARD:
        direction=-1
        
    # Time step (in seconds): temporal resolution of the trajectories - input are dt_data & dt_traj are in hours
    njour=int(24./dt_data)
    npdt=int(dt_data/dt_traj)
    nb_time_step=int((npdt*njour*(trajectories_duration))/24)
    DT=direction*dt_traj*3600. #DT is in seconds
                
    # Initialization of all variables along trajectories for each seed
    TIME_traj=np.zeros((Number_Seeds,nb_time_step+1))  
    LON_traj=np.zeros((Number_Seeds,nb_time_step+1))
    LAT_traj=np.zeros((Number_Seeds,nb_time_step+1))
    P_traj=np.zeros((Number_Seeds,nb_time_step+1))
    U_traj=np.zeros((Number_Seeds,nb_time_step+1))
    V_traj=np.zeros((Number_Seeds,nb_time_step+1))
    W_traj=np.zeros((Number_Seeds,nb_time_step+1))        
    VAR_traj={}
    for i_var in list_var:
        VAR_traj[i_var]=np.zeros((Number_Seeds,nb_time_step+1))

    # Save the first position (seeding) at time 0
    LON_traj[:,0]=x0    
    LAT_traj[:,0]=y0  
    P_traj[:,0]=z0
    
    # Initialize time variable
    TIME_traj[:,0]=initial_time_index*dt_data

    ipdt=0                                                                              
    for i_ech in range(initial_time_index*npdt,npdt*initial_time_index+direction*int((npdt*njour*trajectories_duration)/24),npdt*direction):
        
        print('\t time step : ',i_ech/npdt, file = sys.stderr)
        
        # Get meteorological field for the current time step and the previous/next one
        U_nc=data[list_var_advec[0]].variables[list_var_advec[0]][int(i_ech/npdt),:,:,:]
        V_nc=data[list_var_advec[1]].variables[list_var_advec[1]][int(i_ech/npdt),:,:,:]
        W_nc=data[list_var_advec[2]].variables[list_var_advec[2]][int(i_ech/npdt),:,:,:]
        VAR_nc={}
        for i_var in list_var:
            VAR_nc[i_var]=data[i_var].variables[i_var][int(i_ech/npdt),:,:,:]
            
        U2_nc=data[list_var_advec[0]].variables[list_var_advec[0]][int(i_ech/npdt)+direction,:,:,:]
        V2_nc=data[list_var_advec[1]].variables[list_var_advec[1]][int(i_ech/npdt)+direction,:,:,:]
        W2_nc=data[list_var_advec[2]].variables[list_var_advec[2]][int(i_ech/npdt)+direction,:,:,:]
        VAR2_nc={}
        for i_var in list_var:
            VAR2_nc[i_var]=data[i_var].variables[i_var][int(i_ech/npdt)+direction,:,:,:]
        
        # Initialize other variables
        if ipdt==0:
            print('yes ipdt = 0', file = sys.stderr)
            start = time.perf_counter()
            
            U_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),U_nc,np.array([z0,y0,x0]).T)
            V_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),V_nc,np.array([z0,y0,x0]).T)
            W_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),W_nc,np.array([z0,y0,x0]).T)
            for i_var in list_var:
                VAR_traj[i_var][:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),VAR_nc[i_var],np.array([z0,y0,x0]).T)

        # Computation
        for ipd in range(0,npdt):    
            #print(f'ipd = {ipd}')            
                                                                                        
            ipdt=ipdt+1
            poi=(ipd+1)/npdt

            for ipoin in range(0,Number_Seeds):
                #print(f'ipoin = {ipoin}')   
                lo=LON_traj[ipoin,ipdt-1]
                la=LAT_traj[ipoin,ipdt-1]
                pre=P_traj[ipoin,ipdt-1]
                #print('lo,la,pre',lo,la,pre)                    
                u_tmp=U_traj[ipoin,ipdt-1]
                v_tmp=V_traj[ipoin,ipdt-1]
                w_tmp=W_traj[ipoin,ipdt-1]
                      
                u1=u_tmp
                v1=v_tmp
                w1=w_tmp
                
                for iter in range(1,niter):
                   # print(f'iter = {iter}')
                    coco=np.cos(la*np.pi/180.0);
                    
                    
                    lo_1=(lo+u_tmp*DT/(Ra*coco)*180.0/np.pi)                                 
                    la_1=la+v_tmp*DT/Ra*180.0/np.pi
                    pre_1=pre+w_tmp*DT
                    #print('lo_1,la_1,pre_1',lo_1,la_1,pre_1)
                    # Ensure coordinates lies within the bounds
                    try:
                        lo_1=lo_1-int(lo_1/180.)*360.
                    except ValueError:
                        lo_1=np.nan
                    #pre_1=np.minimum(pre_1,98000.0)
                    #print('lo_1,la_1,pre_1',lo_1,la_1,pre_1)
                    try:
                        #start = time.perf_counter()
                        u2=interp_4d(U_nc,U2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
                        v2=interp_4d(V_nc,V2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
                        w2=interp_4d(W_nc,W2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
                        #end = time.perf_counter()
                        #print(f"Elapsed time 2.5 : {end - start:.2f} seconds")
                    except ValueError:
                        
                        u2=np.nan
                        v2=np.nan
                        w2=np.nan
                    u_tmp=0.5*(u1 + u2)
                    v_tmp=0.5*(v1 + v2)
                    w_tmp=0.5*(w1 + w2)
                coco=np.cos(la*np.pi/180.0);
                lo=lo+u_tmp*DT/(Ra*coco )*180.0/np.pi
                la=la+v_tmp*DT/Ra*180.0/np.pi
                pre=pre+w_tmp*DT

             #   print(w_tmp)
                # Ensure coordinates lies within the bounds
                try:
                    lo=lo-int(lo/180.)*360.
                except ValueError:
                    lo=np.nan
                #pre=np.minimum(pre,99000.0)
                #start = time.perf_counter()
                
                LAT_traj[ipoin,ipdt]=la
                LON_traj[ipoin,ipdt]=lo
                P_traj[ipoin,ipdt]=pre
                TIME_traj[ipoin,ipdt]=initial_time_index*dt_data+direction*ipdt*dt_traj
                
                try:
                    U_traj[ipoin,ipdt]=interp_4d(U_nc,U2_nc,P_nc,LAT_nc,LON_nc,poi,pre,la,lo)
                    V_traj[ipoin,ipdt]=interp_4d(V_nc,V2_nc,P_nc,LAT_nc,LON_nc,poi,pre,la,lo)
                    W_traj[ipoin,ipdt]=interp_4d(W_nc,W2_nc,P_nc,LAT_nc,LON_nc,poi,pre,la,lo)
                except ValueError:
                    U_traj[ipoin,ipdt]=np.nan
                    V_traj[ipoin,ipdt]=np.nan
                    W_traj[ipoin,ipdt]=np.nan
                for i_var in list_var:
                    try:
                        VAR_traj[i_var][ipoin,ipdt]=interp_4d(VAR_nc[i_var],VAR2_nc[i_var],P_nc,LAT_nc,LON_nc,poi,pre,la,lo)
                    except ValueError:
                        VAR_traj[i_var][ipoin,ipdt]=np.nan
            end = time.perf_counter()
            #print(f"Elapsed time : {end - start:.2f} seconds")
    
   
    
    #print(counter)             
    return TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj

def compute_trajectories_ERA5(x0, y0, p0, initial_time_index,
                         LON_nc, LAT_nc, m_nc, data,
                         list_var, list_var_advec, PS, f,
                         trajectories_duration=None,
                         dt_data=6., dt_traj=0.5,
                         niter=4, BACKWARD=True):
    """ Compute Lagrangian trajectories of ERA5 dataset, 

    with model index in vertical coordinates and lagrangian pressure velocity as W - traj_duration is in hours

    f is the function s.t. p_m = f(m, ps, a, b)
    """
    
    print('Computing backward trajectories...: ', file = sys.stderr)

    # Defining number of seeds
    Number_Seeds=len(x0)

    # Defining trajectories direction
    direction=1
    if BACKWARD:
        direction=-1
        
    # Time step (in seconds): temporal resolution of the trajectories - input are dt_data & dt_traj are in hours
    njour=int(24./dt_data)
    npdt=int(dt_data/dt_traj)
    nb_time_step=int((npdt*njour*(trajectories_duration))/24)
    DT=direction*dt_traj*3600. #DT is in seconds
                
    # Initialization of all variables along trajectories for each seed
    TIME_traj=np.zeros((Number_Seeds,nb_time_step+1))  
    LON_traj=np.zeros((Number_Seeds,nb_time_step+1))
    LAT_traj=np.zeros((Number_Seeds,nb_time_step+1))
    m_traj=np.zeros((Number_Seeds,nb_time_step+1))
    P_traj=np.zeros((Number_Seeds,nb_time_step+1))
    U_traj=np.zeros((Number_Seeds,nb_time_step+1))
    V_traj=np.zeros((Number_Seeds,nb_time_step+1))
    W_traj=np.zeros((Number_Seeds,nb_time_step+1))        
    VAR_traj={}
    for i_var in list_var:
        VAR_traj[i_var]=np.zeros((Number_Seeds,nb_time_step+1))


    
    # Save the first position (seeding) at time 0
    LON_traj[:,0]=x0    
    LAT_traj[:,0]=y0  
    P_traj[:,0]=p0
    m0 = np.zeros(np.size(p0))
    PS0=interpn((LAT_nc,LON_nc),PS[int(initial_time_index), :, :],np.array([y0,x0]).T)
    
    start = time.perf_counter()
    for i in np.arange(np.size(m0)) :
        P0_list = f(PS0[i])
        #print(P0.shape)
        #print(m_nc)
        m0[i] = interp1d(P0_list, m_nc, bounds_error=False)(p0[i])
    m_traj[:,0] = m0
    #end = time.perf_counter()
    #print(f"Elapsed time 1 : {end - start:.2f} seconds")
    # Initialize time variable
    

    TIME_traj[:,0]=initial_time_index*dt_data

    ipdt=0                                                                              
    for i_ech in range(initial_time_index*npdt,npdt*initial_time_index+direction*int((npdt*njour*trajectories_duration)/24),npdt*direction):
        
        print('\t time step : ',i_ech/npdt, file = sys.stderr)
        PS_nc = PS[int(i_ech/npdt), :, :] 
        # Get meteorological field for the current time step and the previous/next one
        U_nc=data[list_var_advec[0]].variables[list_var_advec[0]][int(i_ech/npdt),:,:,:]
        V_nc=data[list_var_advec[1]].variables[list_var_advec[1]][int(i_ech/npdt),:,:,:]
        W_nc=data[list_var_advec[2]].variables[list_var_advec[2]][int(i_ech/npdt),:,:,:]
        VAR_nc={}
        for i_var in list_var:
            VAR_nc[i_var]=data[i_var].variables[i_var][int(i_ech/npdt),:,:,:]
        
        PS2_nc = PS[int(i_ech/npdt) + direction, :, :]

        U2_nc=data[list_var_advec[0]].variables[list_var_advec[0]][int(i_ech/npdt)+direction,:,:,:]
        V2_nc=data[list_var_advec[1]].variables[list_var_advec[1]][int(i_ech/npdt)+direction,:,:,:]
        W2_nc=data[list_var_advec[2]].variables[list_var_advec[2]][int(i_ech/npdt)+direction,:,:,:]
        VAR2_nc={}
        for i_var in list_var:
            VAR2_nc[i_var]=data[i_var].variables[i_var][int(i_ech/npdt)+direction,:,:,:]
        
        # Initialize other variables
        if ipdt==0:
            print('yes ipdt = 0', file = sys.stderr)
            U_traj[:,ipdt]=interpn((m_nc,LAT_nc,LON_nc),U_nc,np.array([m0,y0,x0]).T, bounds_error=False)
            V_traj[:,ipdt]=interpn((m_nc,LAT_nc,LON_nc),V_nc,np.array([m0,y0,x0]).T, bounds_error=False)
            W_traj[:,ipdt]=interpn((m_nc,LAT_nc,LON_nc),W_nc,np.array([m0,y0,x0]).T, bounds_error=False)
            for i_var in list_var:
                VAR_traj[i_var][:,ipdt]=interpn((m_nc,LAT_nc,LON_nc),VAR_nc[i_var],np.array([m0,y0,x0]).T, bounds_error=False)
        print(f'ipdt= {ipdt}')

        # Computation
        for ipd in range(0,npdt):    
                                                      
            #print(f'ipd = {ipd}')                                                                            
            start = time.perf_counter()
            ipdt=ipdt+1
            poi=(ipd+1)/npdt

            for ipoin in range(0,Number_Seeds):
                
                #print(f'ipoin = {ipoin}')
                lo=LON_traj[ipoin,ipdt-1]
                la=LAT_traj[ipoin,ipdt-1]
                pre=P_traj[ipoin,ipdt-1]
                m = m_traj[ipoin,ipdt-1]
                #print('lo,la,pre',lo,la,pre)                    
                u_tmp=U_traj[ipoin,ipdt-1]
                v_tmp=V_traj[ipoin,ipdt-1]
                w_tmp=W_traj[ipoin,ipdt-1]
                      
                u1=u_tmp
                v1=v_tmp
                w1=w_tmp
                #print(w1) 
     
                for iter in range(1,niter):
             #       print(f'iter = {iter}')
                    coco=np.cos(la*np.pi/180.0);
                    
                    #start = time.perf_counter()
                    lo_1=(lo+u_tmp*DT/(Ra*coco)*180.0/np.pi)                                 
                    la_1=la+v_tmp*DT/Ra*180.0/np.pi
                    pre_1=pre+w_tmp*DT
                    #if np.shape(pre_1) == (1,) : 
                    #    pre_1 = pre_1[0]

                    #print(m_1)
                    #end = time.perf_counter()

                    #print(f"Elapsed time 2 : {end - start:.2f} seconds")

                 #   start = time.perf_counter()
                    #print('lo_1,la_1,pre_1',lo_1,la_1,pre_1)
                    
                    # Find the corresponding m_1 
                    try:
                        ps_1=interp_3d(PS_nc, PS2_nc, LAT_nc, LON_nc, poi, la_1, lo_1)
                        p_list_1 = f(ps_1)
                        m_1 = interp1d(p_list_1, m_nc, bounds_error=False)(pre_1)
                    except ValueError:
                        m_1 = np.nan
                    
                    # Ensure coordinates lies within the bounds
                    try:
                        lo_1=lo_1-int(lo_1/180.)*360.
                    except ValueError:
                        lo_1=np.nan
                    #pre_1=np.minimum(pre_1,98000.0)
                    #print('lo_1,la_1,pre_1',lo_1,la_1,pre_1)
                    try:
                        u2=interp_4d(U_nc, U2_nc, m_nc, LAT_nc, LON_nc, poi, m_1, la_1, lo_1)
                        v2=interp_4d(V_nc, V2_nc, m_nc, LAT_nc, LON_nc, poi, m_1, la_1, lo_1)
                        w2=interp_4d(W_nc, W2_nc, m_nc, LAT_nc, LON_nc, poi, m_1, la_1, lo_1)
                    except ValueError:
                        u2 = np.nan
                        v2 = np.nan
                        w2 = np.nan
                    #print(w1,w2)
                    u_tmp = 0.5 * (u1 + u2)
                    v_tmp = 0.5 * (v1 + v2)
                    w_tmp = 0.5 * (w1 + w2)
           #         end = time.perf_counter()

                    #print(f"Elapsed time 2.5 : {end - start:.2f} seconds")
                coco = np.cos(la * np.pi / 180.0);
                lo = lo + u_tmp * DT / (Ra * coco)*180.0 / np.pi
                la = la + v_tmp * DT / Ra * 180.0 / np.pi
                pre = pre + w_tmp * DT
                #print(pre)
            #    start = time.perf_counter()
                #end = time.perf_counter()

                #print(f"Elapsed time 3 : {end - start:.2f} seconds")

                try:
                    ps = interp_3d(PS_nc, PS2_nc, LAT_nc, LON_nc, poi, la, lo)
                    p_list = f(ps)
                    m = interp1d(p_list, m_nc, bounds_error=False)(pre)
                except ValueError:
                    m=np.nan
                # Ensure coordinates lies within the bounds
                try:
                    lo=lo-int(lo/180.)*360.
                except ValueError:
                    lo=np.nan
                #pre=np.minimum(pre,99000.0)
                 
                #start = time.perf_counter()
                LAT_traj[ipoin,ipdt]=la
                LON_traj[ipoin,ipdt]=lo
                P_traj[ipoin,ipdt]=pre
                m_traj[ipoin, ipdt] = m
                TIME_traj[ipoin,ipdt]=initial_time_index*dt_data+direction*ipdt*dt_traj
                try:
                    U_traj[ipoin,ipdt]=interp_4d(U_nc,U2_nc,m_nc,LAT_nc,LON_nc,poi,m,la,lo)
                    V_traj[ipoin,ipdt]=interp_4d(V_nc,V2_nc,m_nc,LAT_nc,LON_nc,poi,m,la,lo)
                    W_traj[ipoin,ipdt]=interp_4d(W_nc,W2_nc,m_nc,LAT_nc,LON_nc,poi,m,la,lo)
                except ValueError:
                    U_traj[ipoin,ipdt]=np.nan
                    V_traj[ipoin,ipdt]=np.nan
                    W_traj[ipoin,ipdt]=np.nan
                for i_var in list_var:
                    try:
                        VAR_traj[i_var][ipoin,ipdt]=interp_4d(VAR_nc[i_var],VAR2_nc[i_var],m_nc,LAT_nc,LON_nc,poi,m,la,lo)
                    except ValueError:
                        VAR_traj[i_var][ipoin,ipdt]=np.nan
            end = time.perf_counter()

            print(f"Elapsed time : {end - start:.2f} seconds")
    
    #print(counter)             
    return TIME_traj, LAT_traj, LON_traj, m_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj




#------------------------------------------------------------------------------    
# 4D (space & time) interpolation
#------------------------------------------------------------------------------        
def interp_4d_old(array_t,array_tpdt,P_nc,LAT_nc,LON_nc,poi,pre,la,lo):
    """ Perform space-time interpolation """
    array_coord=(P_nc,LAT_nc,LON_nc)
    interp_coord=np.array([pre,la,lo]).T
    return poi*interpn(array_coord,array_tpdt,interp_coord, bounds_error=False)+(1-poi)*interpn(array_coord,array_t,interp_coord, bounds_error=False)

def interp_4d(array_t,array_tpdt,P_nc,LAT_nc,LON_nc,poi,pre,la,lo):
    """ Perform space-time interpolation """
    array_t_extended,lon_extended=expand_array(array_t,LON_nc)
    array_tpdt_extended=expand_array(array_tpdt,LON_nc,doLon=False)
    array_coord=(P_nc,LAT_nc,lon_extended)
    interp_coord=np.array([pre,la,lo]).T
    return poi*interpn(array_coord,array_tpdt_extended,interp_coord, bounds_error=False)+(1-poi)*interpn(array_coord,array_t_extended,interp_coord, bounds_error=False)

def expand_array(arrayIn,lon,doLon=True):
    """ Expand array in longitudes to account for periodic BC """
    npre,nlat,nlon=arrayIn.shape
    arrayOut=np.zeros((npre,nlat,nlon+2))
    lonOut=np.zeros(nlon+2)
    arrayOut[:,:,1:nlon+1]=arrayIn
    arrayOut[:,:,0]=arrayIn[:,:,nlon-1]
    arrayOut[:,:,nlon+1]=arrayIn[:,:,0]
    if doLon:
        lonOut[1:nlon+1]=lon
        lonOut[0]=lon[nlon-1]-360.
        lonOut[nlon+1]=lon[0]+360.
        return arrayOut,lonOut
    else:
        return arrayOut

def expand_array_3d(arrayIn,lon,doLon=True):
    """ Expand array in longitudes to account for periodic BC """
    nlat,nlon=arrayIn.shape
    arrayOut=np.zeros((nlat,nlon+2))
    lonOut=np.zeros(nlon+2)
    arrayOut[:,1:nlon+1]=arrayIn
    arrayOut[:,0]=arrayIn[:,nlon-1]
    arrayOut[:,nlon+1]=arrayIn[:,0]
    if doLon:
        lonOut[1:nlon+1]=lon
        lonOut[0]=lon[nlon-1]-360.
        lonOut[nlon+1]=lon[0]+360.
        return arrayOut,lonOut
    else:
        return arrayOut

def interp_3d(array_t,array_tpdt,LAT_nc,LON_nc,poi,la,lo):
    """ Perform 2Dspace-time interpolation """
    array_t_extended,lon_extended=expand_array_3d(array_t,LON_nc)
    array_tpdt_extended=expand_array_3d(array_tpdt,LON_nc,doLon=False)
    array_coord=(LAT_nc,lon_extended)
    interp_coord=np.array([la,lo]).T
    return poi*interpn(array_coord,array_tpdt_extended,interp_coord, bounds_error=False) + (1-poi) * interpn(array_coord,array_t_extended,interp_coord, bounds_error=False)
                        
#------------------------------------------------------------------------------    
# Seeding points coordinates: different values possible
#------------------------------------------------------------------------------    
def generate_seeds(Init,Number,Resolution):
    out=[Init]
    for i in range(1,Number,1):
        out.append( out[i-1] + Resolution)
    return np.asarray(out)







#------------------------------------------------------------------------------
# saving Data in NetCDF format
#------------------------------------------------------------------------------
def save_output_data(Root_output,initial_time_index,dt,niter,Nhor,Np,
                     list_var,list_var_advec,
                     TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj):
    """ Write output data to nc file """
    
    Number_Seeds=LAT_traj.shape[0]
    ipdt=LAT_traj.shape[1]-1
    
    print('Saving the trajectories data: ', end='', file = sys.stderr)
        
    ncdf = Dataset(Root_output+f'Traj_time_step_{initial_time_index}_dt{dt}_niter{niter}_Nhor{Nhor}_Np{Np}.nc','w', format='NETCDF4')
    ncdf.createDimension('n_seeds', Number_Seeds)
    ncdf.createDimension('time_ind', ipdt+1)                        
                  
    for i_var,j_var in zip([TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj],
                           ['time','lat','lon','P']+list_var_advec):
        TMP_out = ncdf.createVariable(j_var, 'f8', ('n_seeds','time_ind'))  
        TMP_out[:]=i_var
                
    for i_var in list_var:
        TMP_out = ncdf.createVariable(i_var, 'f8', ('n_seeds','time_ind'))  
        TMP_out[:]=VAR_traj[i_var]  
            
    #Metadata    
    ncdf.initial_time=initial_time_index
   
    
    ncdf.close()
    print('ok', file = sys.stderr)
        
    return

def save_output_data_ERA5(Root_output,initial_time_index, dt, niter, Nhor, Np,
                     list_var,list_var_advec,
                     TIME_traj, LAT_traj, LON_traj, m_traj, P_traj, U_traj, V_traj, W_traj, VAR_traj):
    """ Write output data to nc file """
    
    Number_Seeds=LAT_traj.shape[0]
    ipdt=LAT_traj.shape[1]-1
    
    print('Saving the trajectories data: ', end='', file = sys.stderr)
        
    ncdf = Dataset(Root_output+f'Traj_time_step_{initial_time_index}_dt{dt}_niter{niter}_Nhor{Nhor}_Np{Np}.nc','w', format='NETCDF4')
    ncdf.createDimension('n_seeds', Number_Seeds)
    ncdf.createDimension('time_ind', ipdt+1)                        
                  
    for i_var,j_var in zip([TIME_traj, LAT_traj, LON_traj, m_traj, P_traj, U_traj, V_traj, W_traj],
                           ['time','lat','lon','m','P']+list_var_advec):
        TMP_out = ncdf.createVariable(j_var, 'f8', ('n_seeds','time_ind'))  
        TMP_out[:]=i_var
                
    for i_var in list_var:
        TMP_out = ncdf.createVariable(i_var, 'f8', ('n_seeds','time_ind'))  
        TMP_out[:]=VAR_traj[i_var]  
            
    #Metadata    
    ncdf.initial_time=initial_time_index
   
    
    ncdf.close()
    print('ok', file = sys.stderr)
    
    return 
