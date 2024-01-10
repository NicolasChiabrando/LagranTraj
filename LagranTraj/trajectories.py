#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Trajectographie pour NetCDF
Ph Arbogast & M Wimmer
February 2018
contact: meryl.wimmer@umr-cnrm.fr
"""

#------------------------------------------------------------------------------
# Packages:
#------------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import interpn
from netCDF4 import Dataset
import xarray as xr
import pandas as pd

#------------------------------------------------------------------------------
# Constants:
#------------------------------------------------------------------------------    
Ra=6370000.0 # m : Radius of the Earth
cp=1005.0 # J/kg/K : Specific heat of dry air
Lv=2500000 # J/kg : Latent heat of vaporization at 0◦C
ZR=287.05 #J/K/mol: Gas constant for dry air
ZKP=ZR/cp  # Ra/Cp (pour calcul de la température potentielle)
g=9.81 # N/kg: Acceleration due to gravity at sea level

#------------------------------------------------------------------------------
# Misc functions
#------------------------------------------------------------------------------    
def loop_over_blocks(start_index,end_index,ndates=1,
                     sourceType="ERA5", trajectories_duration=72,blocktype='NA', n_block=50,block_life='onset'):
    """ Function to loop over the blocks """
    if sourceType=='ERA5':
        file="ERA5_dictionary.csv"
    if sourceType=='model2':
        file="con2_smt_dictionary.csv"
    if sourceType=='model3':
        file="con3_smt_dictionary.csv"
    contour_index,time_index=get_block_indices(sourceType=sourceType,blocktype= blocktype, n_block= n_block, block_life= block_life)
    

    for block_index in range(start_index,end_index+1):
        for idate in range(0,ndates):
            compute_traj(contour_index[block_index],time_index[block_index],date_shift=0
                         ,sourceType=sourceType,trajectories_duration=trajectories_duration)
    return
   
# =============================================================================
# def get_block_number(sourceType='ERA5'):
#     """ Function to loop over the blocks """
#     if sourceType=='ERA5':
#         file="ERA5_dictionary.csv"
#     if sourceType=='model2':
#         file="con2_smt_dictionary.csv"
#     if sourceType=='model3':
#         file="con3_smt_dictionary.csv"    
#     contour_index,time_index=get_block_indices(sourceType=sourceType,blocktype= blocktype, n_block= n_block, block_life= block_life)
#     #contour_index,time_index=get_block_indices_pac(sourceType=sourceType)
#     return np.size(contour_index)
#     
# =============================================================================

def get_block_indices(sourceType='ERA5', blocktype='NA', n_block=50,block_life='onset'):   
    if sourceType == 'ERA5':
        flag_file = pd.read_csv('/media/vinita0503/LaCie/Data/ERA5/block_flag_list_ERA5.csv')
        dic_file = pd.read_csv('/media/vinita0503/LaCie/Data/ERA5/ERA5_dictionary.csv')
        dic_file = dic_file.rename(columns={'contour_index': 'Flag'})
    elif sourceType == 'model3':
        flag_file = pd.read_csv('block_flag_list_con3_smt.csv')
        dic_file = pd.read_csv('con3_smt_dictionary.csv')
        dic_file = dic_file.rename(columns={'contour_index': 'Flag'})
    if blocktype == 'All':
        dic_file = dic_file[1:]
    if blocktype == 'NA':
        dic_file = dic_file[dic_file['true_false']]
    elif blocktype == 'PC':
        df = flag_file.groupby('Flag').first().reset_index()
        flag_f_mask = df[(df['Longitude'] < -130) & (df['Longitude'] > -180) & (df['Latitude'] > 35) & (df['Longitude'] < 70)]
        dic_file = dic_file[dic_file['Flag'].isin(flag_f_mask.Flag)]
        # Extracting block timestep and block_index

    block_idx = np.array(dic_file['Flag'].head(n_block))
    int_timestep_idx = np.array(dic_file['initial_time_step'].head(n_block).tolist())
        
    # Selecting those blocks from flag_file to get lat_c and lon_c
    flag_file = flag_file[flag_file.Flag.isin(block_idx)]
        
        # Only the onset flag values nth == stage of block
    if block_life=='onset':
        timestep_idx = []
        timestep_idx = int_timestep_idx
    
    elif block_life=='maintenance':
        time_add = flag_file.groupby('Flag')['Intensity'].apply(lambda x: x.idxmin() - x.index.min())
        timestep_idx=[int_timestep_idx + time_add for int_timestep_idx, time_add in zip(int_timestep_idx, time_add)]
       
    return block_idx , timestep_idx


def dataSource(sourceType="ERA5"):
    
    """ Return data specifics - type should be "ERA5" or "model" """    
    if sourceType=="ERA5":
        #variable names
        list_var_advec=['u','v','w']
        list_var=['ta','pv','pv_anom']
        # Root for files
        Root_input='/media/vinita0503/LaCie/Data/ERA5/'
        Root_output='/media/vinita0503/LaCie/Data/ERA5/'
        # File name for each variable
        File_name={}
        for i_var in list_var_advec+list_var:
            File_name[i_var]='ERA5_'+i_var+'_DJF.nc'
        flagfile='/media/vinita0503/LaCie/Data/ERA5/block_flag_ERA5.nc'

    if sourceType=="model2":
        list_var_advec=['U','V','OMEGA']
        list_var=['T','pv','pv_anom']
        # Root for files
        Root_input='/gpfswork/rech/ptn/usz67lr/block_trajectories/Data_model2_smt/'
        Root_output='/gpfswork/rech/ptn/usz67lr/block_trajectories/traj_data_model_NA_all/'
        # File name for each variable
        File_name={}
        for i_var in list_var_advec+list_var:
            File_name[i_var]='histhf_'+i_var+'_dynamico_smt.nc'
        flagfile='/gpfswork/rech/ptn/usz67lr/block_trajectories/Data_model2_smt/block_flag_con2_smt.nc'
        
    if sourceType=="model3":
        list_var_advec=['U','V','OMEGA']
        list_var=['T','pv','pv_anom']
        # Root for files
        Root_input='/gpfswork/rech/ptn/usz67lr/block_trajectories/Data_model_teq_58_smt/'
        Root_output='/gpfswork/rech/ptn/usz67lr/block_trajectories/traj_data_model_NA_all_mai/'
        # File name for each variable
        File_name={}
        for i_var in list_var_advec+list_var:
            File_name[i_var]='histhf_'+i_var+'_dynamico3_smt.nc'
        flagfile='/gpfswork/rech/ptn/usz67lr/block_trajectories/Data_model_teq_58_smt/block_flag_con3_smt.nc'
    if sourceType=="ECMWF":
        list_var_advec=['u','v','w']
        list_var=['t','pv','q','pv_anom']
        # Root for files
        Root_input='/media/vinita0503/LaCie/Data/ECMWF_IFS/'
        Root_output='/media/vinita0503/LaCie/Data/ECMWF_IFS/'
        File_name = 'TC639_exp_n.nc'
        
        flagfile='/media/vinita0503/LaCie/Data/ECMWF_IFS/Block_flag_TC639_exp_cf.nc'   
    return list_var_advec,list_var,Root_input,Root_output,File_name,flagfile

#------------------------------------------------------------------------------
# Computing trajectories for a given block ID
#------------------------------------------------------------------------------
def compute_traj(tracking_ID,date_index,date_shift=0,sourceType="ERA5",
                 trajectories_duration=None):
    """ tracking_ID is the index of the flag contour, trajectory duration is in hours"""

    BACKWARD=True     # Backward or Forward trajectories?
    SAVING=True       # Saving output ?
    PLOT_TRAJ=True    # Plot Trajectories ?
    initial_time_index=date_index+date_shift       # Time step !

    # Source specific variables
    list_var_advec,list_var,Root_input,Root_output,File_name,flagfile=dataSource(sourceType=sourceType)
    
    #------------------------------------------------------------------------------
    # Seeding
    #------------------------------------------------------------------------------
    x0,y0,z0=get_seeds(initial_time_index,flagfile,tracking_ID)
    

    #------------------------------------------------------------------------------
    # Read data
    #------------------------------------------------------------------------------
    if sourceType == "ECMWF":
        LON_nc,LAT_nc,P_nc,data=read_data_ecmwf(File_name,Root_input,list_var,list_var_advec)
    else:
        LON_nc,LAT_nc,P_nc,data=read_data(File_name,Root_input,list_var,list_var_advec)

    #------------------------------------------------------------------------------
    # Trajectories calculation
    #------------------------------------------------------------------------------
    TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj=\
    compute_trajectories(x0,y0,z0,initial_time_index,
                             LON_nc,LAT_nc,P_nc,data,
                             list_var,list_var_advec,
                             trajectories_duration=trajectories_duration,
                             dt_data=12.,dt_traj=0.5,
                             niter=4,BACKWARD=True)
    
    
# =============================================================================
#     TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj=\
#         compute_trajectories_old(x0,y0,z0,initial_time_index,
#                              LON_nc,LAT_nc,P_nc,data,
#                              list_var,list_var_advec,
#                              trajectories_duration=72,npdt=12,njour=4,niter=4,BACKWARD=True)
# =============================================================================
    
 
         


    #------------------------------------------------------------------------------
    # saving Data in NetCDF format
    #------------------------------------------------------------------------------
    if SAVING:
        save_output_data(SAVING,Root_output,tracking_ID,initial_time_index,
                         list_var,list_var_advec,
                         TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj)

    print('\t Finish !')
    
    return 

#------------------------------------------------------------------------------
# Read data (two modes possible: xarray or netcdf4)
#------------------------------------------------------------------------------

def read_data_xarray(File_name,Root_input,list_var,list_var_advec):
    """ Read data """
    print('Reading data with xarray...')
    if isinstance(File_name,dict):
        data={}
        for i_var in list_var_advec+list_var:
            ds=xr.open_dataset(Root_input+File_name[i_var])
            data[i_var]=np.array(ds[i_var])
    else:
        data={}
        for i_var in list_var_advec+list_var:
            ds=xr.open_dataset(Root_input+File_name)
            data[i_var]=(ds[i_var])
            
    LON_nc=ds['lon']
    LAT_nc=ds['lat']        
    P_nc=ds['level'][:]*100
    
    return LON_nc,LAT_nc,P_nc,data

def read_data(File_name,Root_input,list_var,list_var_advec):
    """ Read data """
    print('Read data: ', end='')
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
    P_nc=data[list_var_advec[0]].variables['level'][:]*100
    
    return LON_nc,LAT_nc,P_nc,data

def read_data_ecmwf(File_name,Root_input,list_var,list_var_advec):
    """ Read data """
    print('Read data: ', end='')
    data={}
    
    for i_var in list_var_advec+list_var:
        data[i_var]=   Dataset(Root_input+File_name,format='NETCDF4')
            
            
    LON_nc=data[list_var_advec[0]]['longitude'][:]
    LAT_nc=data[list_var_advec[0]]['latitude'][:]        
    P_nc=data[list_var_advec[0]]['isobaricInhPa'][:]*100
    return LON_nc,LAT_nc,P_nc,data
#------------------------------------------------------------------------------
# rename dims and variable name
#------------------------------------------------------------------------------
def rename_coords_var(data,variable=None):
    list_dims=list(data.dims)
    list_var=list(data.data_vars)
    return data.rename({list_dims[0]: 'time',list_dims[1]:'lat' ,list_dims[2]:'lon',list_var[0]:variable})  
    
#------------------------------------------------------------------------------
# saving Data in NetCDF format
#------------------------------------------------------------------------------
def save_output_data(SAVING,Root_output,tracking_ID,initial_time_index,
                     list_var,list_var_advec,
                     TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj):
    """ Write output data to nc file """
    
    Number_Seeds=LAT_traj.shape[0]
    ipdt=LAT_traj.shape[1]-1
    
    print('Saving the trajectories data: ', end='')
        
    ncdf = Dataset(Root_output+'Traj_ID'+str(tracking_ID)+'time_step_'+str(initial_time_index)+'.nc','w', format='NETCDF4')
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
    print('ok')
        
    return

#------------------------------------------------------------------------------
# Seeding
#------------------------------------------------------------------------------
def get_seeds(initial_time_index,flagfile,tracking_ID,sourceType='ECMWF'):
    """ Create seeds """

    print('Creating seeds: ',end='')

    # Initialize variables
    VERIFY_SEEDING=True # -> Create a plot with the seeding
    CROSS_SECTION=False # if False : seeding in a box;
                        # if True : Number_Seeds_Longitude=1 or Number_Seeds_Latitude=1 or Number_Seeds_Longitude=Number_Seeds_Latitude
    Initial_Longitude=-180
    Initial_Latitude= 90.#45                 
    Initial_Pressure=30000 #80000              
    Number_Seeds_Longitude=361
    Number_Seeds_Latitude=91
    Number_Seeds_Pressure=1
    Seeding_Resolution_Longitude=1.0
    Seeding_Resolution_Latitude=-1.0 #0.5
    Seeding_Resolution_Pressure=5000

    # Number of seeds
    Number_Seeds=Number_Seeds_Longitude*Number_Seeds_Latitude*Number_Seeds_Pressure  

    if CROSS_SECTION:
        if Number_Seeds_Longitude==Number_Seeds_Latitude:
            Number_Seeds=Number_Seeds_Pressure*Number_Seeds_Latitude
                
    Latitude_values  =generate_seeds(Initial_Latitude,  Number_Seeds_Latitude,  Seeding_Resolution_Latitude)
    Longitude_values =generate_seeds(Initial_Longitude, Number_Seeds_Longitude, Seeding_Resolution_Longitude)
    Pressure_values  =generate_seeds(Initial_Pressure,  Number_Seeds_Pressure,  Seeding_Resolution_Pressure)
    
    #creating pv_box at initial timestep
    flag_data=xr.open_dataset(flagfile)
    flag_data = rename_coords_var(flag_data,variable='flag')
    block=flag_data.flag
    
        
    lon_slice=np.arange(Initial_Longitude,Number_Seeds_Longitude*Seeding_Resolution_Longitude
                        +Initial_Longitude,Seeding_Resolution_Longitude)
        
    lat_slice=np.arange(Initial_Latitude,Number_Seeds_Latitude*Seeding_Resolution_Latitude
                        +Initial_Latitude,Seeding_Resolution_Latitude)
        
    block=block[initial_time_index,:,:].sel(lat=lat_slice,lon=lon_slice, method='nearest')
        
    index=np.arange(10000,105000,5000)
    if sourceType == 'ECMWF':
        index = index[::-1]
    block_all=xr.concat([block,block,block,block,block,block,block,block,block,block,
                         block,block,block,block,block,block,block,block,block] ,pd.Index(index, name="level"))
        
        
    # applying condition to statisfy PV anom cri
    block_all=block_all.sel(level=Pressure_values)
    
    block_all=block_all.data.reshape(Number_Seeds)
    
    #to eliminate seeding point not satisfying PV cri    
    cri= block_all==tracking_ID
    
        
    #Seeding points coordinates: array with all coordinates
    if CROSS_SECTION:
        if Number_Seeds_Latitude==1:
            x,z=np.meshgrid(Longitude_values,Pressure_values)
            y=Latitude_values[0]
        elif Number_Seeds_Longitude==1:
            y,z=np.meshgrid(Latitude_values,Pressure_values)
            x=Longitude_values[0]
        elif Number_Seeds_Pressure==1:
            x,y=np.meshgrid(Longitude_values,Latitude_values)
            z=Pressure_values[0]
        else:
            x,z=np.meshgrid(Longitude_values,Pressure_values)
            y,z=np.meshgrid(Latitude_values,Pressure_values)        
    else:
        y,z,x=np.meshgrid(Latitude_values,Pressure_values,Longitude_values)

    #box area finish
                
    if isinstance(x,np.ndarray):x=x.reshape(Number_Seeds)
    if isinstance(y,np.ndarray):y=y.reshape(Number_Seeds)
    if isinstance(z,np.ndarray):z=z.reshape(Number_Seeds)
    # criteria to the meshgrid 
    x=x[cri]
    y=y[cri]
    z=z[cri]
        
    #redefining number od seeds 
    Number_Seeds=x.shape[0]
    print('Number of seeds=',Number_Seeds)    
    #  Graphes to verify the seeding
    if VERIFY_SEEDING:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3D(x,y,z,c=z, edgecolors='black', cmap='Greens')
        ax.invert_zaxis()
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        ax.set_zlabel('Pressure [Pa]')
        plt.show()
           
    return x,y,z 

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
    
    print('Computing backward trajectories...: ')

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
        
        print('\t time step : ',i_ech/npdt)
        
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
            print('yes ipdt = 0')
            
            U_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),U_nc,np.array([z0,y0,x0]).T)
            V_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),V_nc,np.array([z0,y0,x0]).T)
            W_traj[:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),W_nc,np.array([z0,y0,x0]).T)
            for i_var in list_var:
                VAR_traj[i_var][:,ipdt]=interpn((P_nc,LAT_nc,LON_nc),VAR_nc[i_var],np.array([z0,y0,x0]).T)

        # Computation
        for ipd in range(0,npdt):    
                                                      
                                                                                        
            ipdt=ipdt+1
            poi=(ipd+1)/npdt

            for ipoin in range(0,Number_Seeds):
                   
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
                        u2=interp_4d(U_nc,U2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
                        v2=interp_4d(V_nc,V2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
                        w2=interp_4d(W_nc,W2_nc,P_nc,LAT_nc,LON_nc,poi,pre_1,la_1,lo_1)
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

                # Ensure coordinates lies within the bounds
                try:
                    lo=lo-int(lo/180.)*360.
                except ValueError:
                    lo=np.nan
                #pre=np.minimum(pre,99000.0)
                
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
    
   
    
    #print(counter)             
    return TIME_traj, LAT_traj, LON_traj, P_traj, U_traj, V_traj, W_traj,VAR_traj

#------------------------------------------------------------------------------    
# 4D (space & time) interpolation
#------------------------------------------------------------------------------        
def interp_4d_old(array_t,array_tpdt,P_nc,LAT_nc,LON_nc,poi,pre,la,lo):
    """ Perform space-time interpolation """
    array_coord=(P_nc,LAT_nc,LON_nc)
    interp_coord=np.array([pre,la,lo]).T
    return poi*interpn(array_coord,array_tpdt,interp_coord)+(1-poi)*interpn(array_coord,array_t,interp_coord)

def interp_4d(array_t,array_tpdt,P_nc,LAT_nc,LON_nc,poi,pre,la,lo):
    """ Perform space-time interpolation """
    array_t_extended,lon_extended=expand_array(array_t,LON_nc)
    array_tpdt_extended=expand_array(array_tpdt,LON_nc,doLon=False)
    array_coord=(P_nc,LAT_nc,lon_extended)
    interp_coord=np.array([pre,la,lo]).T
    return poi*interpn(array_coord,array_tpdt_extended,interp_coord)+(1-poi)*interpn(array_coord,array_t_extended,interp_coord)

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

                        
#------------------------------------------------------------------------------    
# Seeding points coordinates: different values possible
#------------------------------------------------------------------------------    
def generate_seeds(Init,Number,Resolution):
    out=[Init]
    for i in range(1,Number,1):
        out.append( out[i-1] + Resolution)
    return np.asarray(out)


