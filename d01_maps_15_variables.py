#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 14:09:06 2021

@author: kostas
"""

"""
 THE PRESENT CODE OPENS THE WRF FILES OF DOMAIN 01  AND THEN IT 
 CONSTRUCTS MULTIPLE MAPS AND IT SAVES THEM INTO DIFFERENT FILES.

ATTENTION: WHEN SOMEONE NEEDS TO USE THE PRESENT CODE ON HIS MACHINE HE HAS TO
           CHANGE THE !!! LINES WHICH ARE:
               d01_basic_path =......(101 line) 
               HE MAY NEED TO CHANGE THE REST OF THE PATHS THERE, BUT EVERETHING ELSE IS
               ALREADY SET.
               
"""
#%% Import Libraries // Open files and variables // Set Lon-Lat for domain // Set interpolation levels

# -----------------------  Importing Libraries  -------------------------------
import matplotlib
import glob
import os
from matplotlib.cm import get_cmap   
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from wrf import (getvar, to_np, get_cartopy,interplevel, smooth2d, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair,ALL_TIMES,vinterp)
from netCDF4 import Dataset
from metpy.units import units
import metpy.calc as mpcalc
from datetime import datetime
from tqdm import tqdm
#  ---------------------------------------------------------------------------
start_time = datetime.now()
print ("\n (d01 Visualization) Initialization time = ",start_time)
# -----------------------------------------------------------------------------
print("\n  Open files and variables  //  Set Lon-Lat for domain  //  Set interpolation levels \n")

pwd_directory = os.getcwd()  # the current directory where wrfout files exist
# Now we get the files we want from the new directory
filenames = [f for f in os.listdir(pwd_directory) if f.startswith('wrfout_d01')]
# Define the WRF GRID RESOLUTION for the current domain
wrf_grid_resolution =' Res:30km'   # !!! need to change that point

# Open the wrf file and keep the variable we are interested in
wrf_file = Dataset(filenames[0])  

# Lon-Lat processing
lats = wrf_file['XLAT'][:,:,0][0]
lons = wrf_file['XLONG'][:,0,:][0] 
time = wrf_file['XTIME'][:]
lons.min()
lons.max()
lats.max()
lats.min()
space_lat = (lats.max() - lats.min())/len(lats)
space_lon = (lons.max() - lons.min())/len(lons)
lat_made = np.arange(lats.min(), lats.max(), space_lat) 
lon_made = np.arange(lons.min(), lons.max(), space_lon) 
# The dimensions of the arrays must be compatible.

# Mean values for lon, lat for my data:
lon_0 = lon_made.mean()
lat_0 = lat_made.mean()

'''  Extract the wrf variables we want  ''' 
# wrf_Terrain = getvar(wrf_file, "ter", timeidx=ALL_TIMES)   #!!!NEW    Model Terrain Height            
''' Atmospheric Parameters  '''
p              = getvar(wrf_file, "pressure",timeidx=ALL_TIMES)        # Full model pressure [hPa]
z              = getvar(wrf_file, "z", units="dm",timeidx=ALL_TIMES)   # Geopotential height (m)
ua             = getvar(wrf_file, "ua", units="kt",timeidx=ALL_TIMES)  # V-component of Wind on Mass Points
va             = getvar(wrf_file, "va", units="kt",timeidx=ALL_TIMES)   # W-component of Wind on Mass Points
wspd           = getvar(wrf_file, "wspd_wdir", units="kts",timeidx=ALL_TIMES)[0,:] # wind speed
temp           = getvar(wrf_file, "tc", timeidx=ALL_TIMES)              #Temperature in Celsius
slp            = getvar(wrf_file, "slp", timeidx=ALL_TIMES)           # we keep the Sea Level Pressure [hPa] 
cloud_fraction = getvar(wrf_file, "cfrac", timeidx=ALL_TIMES)
# max_reflect    = getvar(wrf_file, "mdbz",timeidx=ALL_TIMES)           # Maximum Reflectivity [dBZ] 
# # new to add
# cape_2d      = getvar(wrf_file, "cape_2d",timeidx=ALL_TIMES)           #!!!NEW   	Returns 2D fields mcape/mcin/lcl/lfc
# cape_3d      = getvar(wrf_file, "cape_3d",timeidx=ALL_TIMES)           #!!!NEW  Returns 3D fields cape and cin (J kg-1)
# temp_wb      = getvar(wrf_file, "twb",timeidx=ALL_TIMES)           #!!!NEW  Wet Bulb Temperature (K)
# temp_wb      = temp_wb-273.15                                      #!!! convert temp (K --> C ) (ONLY FOR TEMPERATURE)
# T            = getvar(wrf_file, "T", timeidx=ALL_TIMES)           #!!!NEW Temperature IN K
# temp_td      = getvar(wrf_file, "td",timeidx=ALL_TIMES)            #!!!NEW Dew Point Temperature in Celsius
# rh           = getvar(wrf_file, "rh", timeidx=ALL_TIMES)           #!!!NEW  Relative humidity (%)
''' Surface Parameters '''
wrf_rh2 = getvar(wrf_file, "rh2", timeidx=ALL_TIMES)                  # Relative humidity at 2m (%)  
wind_wrf_U10 = getvar(wrf_file, "U10", timeidx=ALL_TIMES)  
wind_wrf_V10 = getvar(wrf_file, "V10", timeidx=ALL_TIMES)
wrf_temp_2m = getvar(wrf_file, "T2", timeidx=ALL_TIMES)               # temperature at 2m(K)
wrf_temp_2m = wrf_temp_2m-273.15                                      # convert temp (K --> C ) (ONLY FOR TEMPERATURE)
wrf_snow_wat_equiv = getvar(wrf_file, "SNOW", timeidx=ALL_TIMES)       #!!!NEW   Snow water equivalent ( kg m-2  )
# wrf_snowh = getvar(wrf_file, "SNOWH", timeidx=ALL_TIMES)              # Physical snow depth (m)
# wrf_snowh = wrf_snowh*100                                             # convert snow (m --> cm ) (ONLY FOR SNOW)
# # new to add
# wrf_snow_daily = getvar(wrf_file, "SNOWNC", timeidx=ALL_TIMES)         #!!!NEW   Daily total snow and ice (mm )         
# wrf_sst = getvar(wrf_file, "SST", timeidx=ALL_TIMES)                   #!!!NEW   Sea surface temperature (K)               
# wrf_sst = wrf_sst-273.15                                               #!!!NEW   convert temp (K --> C ) (ONLY FOR TEMPERATURE)
# wrf_mix_rat_2m = getvar(wrf_file, "Q2", timeidx=ALL_TIMES)             #!!!NEW   Water vapor mixing ratio at 2m  (kg kg-1)
# wrf_hail_daily = getvar(wrf_file, "HAILNC", timeidx=ALL_TIMES)         #!!!NEW   Daily total hail (mm)  
# wrf_td_2m = getvar(wrf_file, "td2", timeidx=ALL_TIMES)                 #!!!NEW   2m Dew Point Temperature in Celsius
# Rain
wrf_RAINC  = getvar(wrf_file, "RAINC", timeidx=ALL_TIMES)           # Daily total convective precipitation (mm) (/25.4-->from ncl script)
wrf_RAINNC = getvar(wrf_file, "RAINNC", timeidx=ALL_TIMES)        # Daily total non-convective precipitation(mm) (/25.4-->from ncl script)
wrf_total_precip  = wrf_RAINC + wrf_RAINNC
# In order to check for each wrf_variable, which dimension it's based on we give >> len(wrf_variable)  

# Define the pressure levels for interpolation
interp_levels = [500,700,850]

#%%  Set paths for every variable
print("\n  Set paths for every variable \n")



d01_basic_path="/path_to_WRF_visualization/d01_temp_maps"   #!!!

wind_500_dir = "/wind_500"
temp_500_dir = "/temp_500"
wind_700_dir = "/wind_700"
temp_700_dir = "/temp_700"
wind_850_dir = "/wind_850"
temp_850_dir = "/temp_850"
# cloud_fraction_low_clouds_dir     ="/cloud_fraction/low_clouds"
# cloud_fraction_medium_clouds_dir  ="/cloud_fraction/medium_clouds"
# cloud_fraction_high_clouds_dir    ="/cloud_fraction/high_clouds"
cloud_fraction_all_clouds_dir     ="/cloud_fraction"
wind_10_slp_dir = "/wind_10_slp"
t2m_dir         = "/t2m"
snow_dir        = "/snow"
rh2_dir         = "/rh2"
rain_dir        = "/rain"
Max_refl_dir    = "/Max_refl"

variable_directory = [wind_500_dir, temp_500_dir, wind_700_dir, temp_700_dir,
                      wind_850_dir, temp_850_dir, cloud_fraction_all_clouds_dir, 
                      wind_10_slp_dir, t2m_dir, snow_dir, rh2_dir, rain_dir, Max_refl_dir]

#%% We set the cities-villages we want to check
print("\n We set the cities-villages we want to check \n")

cities    = ['put the desired cities on your domain.......']

city_lat  = ["add the corresponding lats"]

city_lon  = ["add the corresponding lons"]
#%% COLORBARS FOR clouds-rain-snow 
nws_precip_colors = [ "#fdfdfd", "#b3f8f7", "#b3f8f7", "#68f1f0", "#36edeb",
                      "#04e9e7", "#03d1cf", "#36b9ed", "#04a8e9", "#0397d1",
                      "#0386ba", "#0275a3", "#fffa99", "#e5e189", "#ccc87a",
                      "#ffda99", "#e5c489", "#f48536", "#db7730", "#f43636", 
                      "#db3030", "#f43665", "#db305a", "#FF7C00", "#FF5300",
                      "#FF00EA", "#e500d2", "#cc00bb", "#c400cc", "#b200cc"] 
precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

snow_colors = ["#fdfdfd", "#c1da46", "#bbd632", "#b2d119", "#aacc00", "#92d119",
               "#86cc00", "#41d119", "#2ccc00", "#00cc6d", "#00b762", "#00a357",
               "#008e4c", "#007a41", "#006636", "#00512b", "#003d20"]
snow_colormap = matplotlib.colors.ListedColormap(snow_colors)

high_cloud_colors = ["#ffffff", "#d7e3d1", "#c3d5bb","#afc8a4", "#9bba8e",
                     "#87ac77", "#739f60", "#5f914a", "#4b8333", "#38761d",
                     "#326a1a", "#2c5e17", "#214611", "#1c3b0e", "#162f0b",
                     "#102308"]
high_cloud_colormap = matplotlib.colors.ListedColormap(high_cloud_colors)

med_cloud_colors = ["#ffffff", "#fcd9d6", "#fbc6c2", "#fab3ae", "#f9a19a" , 
                    "#f88e86", "#f77b72", "#f6685e", "#f5554a", "#f44336",
                    "#db3c30", "#c3352b",  "#aa2e25", "#922820","#7a211b",
                    "#611a15"]
med_cloud_colormap = matplotlib.colors.ListedColormap(med_cloud_colors)

low_cloud_colors = ["#ffffff", "#bedaef", "#a9ceea", "#94c2e5",  "#7eb6e0", 
                    "#69aadb", "#3e92d1", "#2986cc", "#2478b7",  "#206ba3",
                    "#1c5d8e", "#18507a", "#144366", "#103551",  "#0c283d", 
                    "#081a28"]
low_cloud_colormap = matplotlib.colors.ListedColormap(low_cloud_colors)

#%% We keep the min-max values
print("\n We keep the min-max values of every parameter FOR THIS WRF RUN \n")
# We keep the min-max values of every parameter in lists in order to set the 
# contour levels for each one of them FOR THIS WRF RUN

wind_10_slp_max_list = []
wind_10_slp_min_list = []
#
wind_500_max_list = []
wind_500_min_list = []
temp_500_max_list = []
temp_500_min_list = []
z_500_max_list    = []
z_500_min_list    = []
wind_700_max_list = []
wind_700_min_list = []
temp_700_max_list = []
temp_700_min_list = []
z_700_max_list    = []
z_700_min_list    = []
wind_850_max_list = []
wind_850_min_list = []
temp_850_max_list = []
temp_850_min_list = []
z_850_max_list    = []
z_850_min_list    = []
#
slp_max_list      = []
slp_min_list      = []
t2m_max_list      = []
t2m_min_list      = []
# snow_max_list     = []
# snow_min_list     = []
rain_max_list     = []
# Max_refl_max_list = []
# Max_refl_min_list = []

for i in tqdm(range(12,len(time))): 
     #  Vertical Interpolation        
    """==== [Vertical Interpolation for some atmospheric parameters] =====
      ------  Interpolation for 500,700,850 hPa  pressure levels   -------"""  
      
    # Interpolate z to pressure levels (hPa)
    interp_z = vinterp(wrf_file,
                   field=z[i],
                   vert_coord="pressure",
                   interp_levels=interp_levels,
                   extrapolate=True,
                   field_type="pressure_hpa",
                   log_p=True)
    
    # Interpolate wspd to pressure levels (hPa)
    interp_wspd = vinterp(wrf_file,
                   field=wspd[i],
                   vert_coord="pressure",
                   interp_levels=interp_levels,
                   extrapolate=True,
                   field_type="pressure_hpa",
                   log_p=True)
    
    # Interpolate ua to pressure levels (hPa)
    interp_ua = vinterp(wrf_file,
                   field=ua[i],
                   vert_coord="pressure",
                   interp_levels=interp_levels,
                   extrapolate=True,
                   field_type="pressure_hpa",
                   log_p=True)
    
    # Interpolate va to pressure levels (hPa)
    interp_va = vinterp(wrf_file,
                   field=va[i],
                   vert_coord="pressure",
                   interp_levels=interp_levels,
                   extrapolate=True,
                   field_type="pressure_hpa",
                   log_p=True)    
    
    # Interpolate temp to pressure levels (hPa)
    interp_temp = vinterp(wrf_file,
                   field=temp[i],
                   vert_coord="pressure",
                   interp_levels=interp_levels,
                   extrapolate=True,
                   field_type="pressure_hpa",
                   log_p=True)    
    
    # For 500,700,850 hPa interpolation:    
    # ============================
    z_500 = interp_z[0]
    z_700 = interp_z[1]
    z_850 = interp_z[2]
    
    u_500 =  interp_ua[0]
    u_700 =  interp_ua[1]
    u_850 =  interp_ua[2]

    v_500 =  interp_va[0]    
    v_700 =  interp_va[1]    
    v_850 =  interp_va[2]       
    
    wspd_500 = interp_wspd[0]
    wspd_700 = interp_wspd[1]
    wspd_850 = interp_wspd[2]
    
    temp_500 = interp_temp[0]
    temp_700 = interp_temp[1]
    temp_850 = interp_temp[2]
    # ==========================================   
    
    # ===========================
    # For wind at 10m 
    # we give units to our data
    U10_wrf = wind_wrf_U10[i].values  * units('m/s')
    #   ---------   Wind V10  ------------------------
    # we give units to our data
    V10_wrf = wind_wrf_V10[i].values  * units('m/s')    
    # -----------  Wind speed construction of U10,V10 WRF --------
    wrf_WSpeed_10m = mpcalc.wind_speed(U10_wrf, V10_wrf)
    wrf_WSpeed_10m =(wrf_WSpeed_10m*1.9438444924406)* units('s/m*knots/knots') 
    wrf_WSpeed_10m     
    #====================
    # keep min-max values
    wind_10_slp_max_list.append(wrf_WSpeed_10m.max())
    wind_10_slp_min_list.append(wrf_WSpeed_10m.min())    
    #
    wind_500_max_list.append(wspd_500.max().values)
    wind_500_min_list.append(wspd_500.min().values)
    
    wind_700_max_list.append(wspd_700.max().values)
    wind_700_min_list.append(wspd_700.min().values)
    
    wind_850_max_list.append(wspd_850.max().values)
    wind_850_min_list.append(wspd_850.min().values)
    
    temp_500_max_list.append(temp_500.max().values)
    temp_500_min_list.append(temp_500.min().values)
    
    temp_700_max_list.append(temp_700.max().values)
    temp_700_min_list.append(temp_700.min().values)
    
    temp_850_max_list.append(temp_850.max().values)
    temp_850_min_list.append(temp_850.min().values)
    
    z_500_max_list.append(z_500.max().values)
    z_500_min_list.append(z_500.min().values)
    
    z_700_max_list.append(z_700.max().values)
    z_700_min_list.append(z_700.min().values)
    
    z_850_max_list.append(z_850.max().values)
    z_850_min_list.append(z_850.min().values)
    #    
    slp_max_list.append((slp[i].max()).values)
    slp_min_list.append((slp[i].min()).values)
    
    t2m_max_list.append((wrf_temp_2m[i].max()).values)
    t2m_min_list.append((wrf_temp_2m[i].min()).values)
    
    # snow_max_list.append((wrf_snow_wat_equiv[i].max()).values)
    # snow_min_list.append((wrf_snow_wat_equiv[i].min()).values)
    
    # Max_refl_max_list.append((max_reflect[i].max()).values)
    # Max_refl_min_list.append((max_reflect[i].min()).values)
    
    rain_max_list.append((wrf_total_precip[i].max()).values)

# We name the min-max values as variables in order to use them in contour
# levels on our plots   
slp_max = int(max(slp_max_list))
slp_min = int(min(slp_min_list))    

wind_10_slp_max = int(max(wind_10_slp_max_list))
wind_10_slp_min = int(min(wind_10_slp_min_list))
#
wind_500_max = int(max(wind_500_max_list))
wind_500_min = int(min(wind_500_min_list))

wind_700_max = int(max(wind_700_max_list))
wind_700_min = int(min(wind_700_min_list))

wind_850_max = int(max(wind_850_max_list))
wind_850_min = int(min(wind_850_min_list))

temp_500_max = int(max(temp_500_max_list))
temp_500_min = int(min(temp_500_min_list))

temp_700_max = int(max(temp_700_max_list))
temp_700_min = int(min(temp_700_min_list))

temp_850_max = int(max(temp_850_max_list))
temp_850_min = int(min(temp_850_min_list))

z_500_max = int(max(z_500_max_list))
z_500_min = int(min(z_500_min_list))

z_700_max = int(max(z_700_max_list))
z_700_min = int(min(z_700_min_list))

z_850_max = int(max(z_850_max_list))
z_850_min = int(min(z_850_min_list))
#
wrf_temp_2m_max = int(max(t2m_max_list))
wrf_temp_2m_min = int(min(t2m_min_list))

# snow_max = int(max(snow_max_list))
snow_min = 0

rh2_max = 100
rh2_min = 0

rain_max = int(max(rain_max_list))
rain_min =0

# Max_refl_max = int(max(Max_refl_max_list))
# Max_refl_min = int(min(Max_refl_min_list))

### the next one overcomes the memory problems with plots 
plt.rcParams.update({'figure.max_open_warning': 0}) 

#%% Delete previous files (plots)
print(" \n ------  Delete maps from previous wrf runs   -------------\n")

""" As I get into a directory, in which I have to save the new files (plots)
        of the latest wrf run, I have to delete the previous files"""
        
for i in range (0,len(variable_directory)):
    files = glob.glob(d01_basic_path + variable_directory[i]+'/*')   #!!!
    for f in files:
        os.remove(f) 
    

print("\n ------- Construction of wrfout maps for d01 domain  ----------\n")

#%% Make map plots in a loop of runs based on time
for i in tqdm(range(12,len(time))):  
     
    #%%   Vertical Interpolation        
    """==== [Vertical Interpolation for some atmospheric parameters] =====
      ------  Interpolation for 500,700,850 hPa  pressure levels   -------"""  
      
    # Interpolate z to pressure levels (hPa)
    interp_z = vinterp(wrf_file,
                    field=z[i],
                    vert_coord="pressure",
                    interp_levels=interp_levels,
                    extrapolate=True,
                    field_type="pressure_hpa",
                    log_p=True)
    
    # Interpolate wspd to pressure levels (hPa)
    interp_wspd = vinterp(wrf_file,
                    field=wspd[i],
                    vert_coord="pressure",
                    interp_levels=interp_levels,
                    extrapolate=True,
                    field_type="pressure_hpa",
                    log_p=True)
    
    # Interpolate ua to pressure levels (hPa)
    interp_ua = vinterp(wrf_file,
                    field=ua[i],
                    vert_coord="pressure",
                    interp_levels=interp_levels,
                    extrapolate=True,
                    field_type="pressure_hpa",
                    log_p=True)
    
    # Interpolate va to pressure levels (hPa)
    interp_va = vinterp(wrf_file,
                    field=va[i],
                    vert_coord="pressure",
                    interp_levels=interp_levels,
                    extrapolate=True,
                    field_type="pressure_hpa",
                    log_p=True)    
    
    # Interpolate temp to pressure levels (hPa)
    interp_temp = vinterp(wrf_file,
                    field=temp[i],
                    vert_coord="pressure",
                    interp_levels=interp_levels,
                    extrapolate=True,
                    field_type="pressure_hpa",
                    log_p=True)
    
    
    # For 500,700,850 hPa interpolation:    
    # ============================
    z_500 = interp_z[0]
    z_700 = interp_z[1]
    z_850 = interp_z[2]
    
    u_500 =  interp_ua[0]
    u_700 =  interp_ua[1]
    u_850 =  interp_ua[2]

    v_500 =  interp_va[0]    
    v_700 =  interp_va[1]    
    v_850 =  interp_va[2]       
    
    wspd_500 = interp_wspd[0]
    wspd_700 = interp_wspd[1]
    wspd_850 = interp_wspd[2]
    
    temp_500 = interp_temp[0]
    temp_700 = interp_temp[1]
    temp_850 = interp_temp[2]
    # ===========================
     
    #%% (1) 500 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])
    """ 500 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])"""
    
    os.chdir(d01_basic_path + wind_500_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
   
    # ---------------  Construct the interpolated map  ----------------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_500)      
    # Get the map projection information
    cart_proj = get_cartopy(z_500)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2) 
    #  ----------------------------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    
    # Add the 500 hPa geopotential height contours
    levels = np.arange(z_500_min-1, z_500_max+1, 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_500),
                           levels=levels, colors="white",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
    
    # Add the wind speed contours
    levels = np.arange(wind_500_min, wind_500_max+1, 2.)
    wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_500),
                                 levels=levels,cmap=get_cmap("rainbow"),
                                 transform=ccrs.PlateCarree())
    cb = plt.colorbar(wspd_contours, ax=ax, orientation="vertical", shrink=0.7)
    cb.set_label('Wind Speed  500hPa (kt)',size=15)
    # Add the 500 hPa wind barbs
    plt.barbs(to_np(lons[::5,::5]), to_np(lats[::5,::5]),
              to_np(u_500[::5, ::5]), to_np(v_500[::5, ::5]),
              transform=ccrs.PlateCarree(), length=6)
     
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_500))
    ax.set_ylim(cartopy_ylim(z_500))    
    # ax.gridlines()  
    # Titles
    plt.title('500 hPa Height (dm), Wind Speed (kt), Barbs (kt) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('500_hPa_Wind_Speed__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()  
    
     #%%   (2) 500 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])
    """ 500 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])"""

    os.chdir(d01_basic_path + temp_500_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT  
      
    # ---------------  Construct the interpolated map  ----------------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_500)      
    # Get the map projection information
    cart_proj = get_cartopy(z_500)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2) 
    #  ----------------------------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
   
    # Now we build the contours for geopotential height at 500hPa
    levels = np.arange(z_500_min-1, z_500_max+1, 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_500),
                           levels=levels, colors="m",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Now we build the contour & contourf for temperature at 500hpa 
    temp_500_levels = np.arange(temp_500_min,temp_500_max+1,1) # specify the step of contours
    temp_500_contour =ax.contour(to_np(lons), to_np(lats),temp_500,colors='k', levels=temp_500_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())       
    temp_500_contourf = ax.contourf(to_np(lons), to_np(lats), temp_500,levels=temp_500_levels,
                                cmap="jet", zorder=1, transform = ccrs.PlateCarree())    
       
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)    
    plt.clabel(temp_500_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(temp_500_contourf, shrink=0.5) # shrink= size of the colorbar
    cb.set_label('Temperature (\u2103)',size=15)    
    # add gridlines
    ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')
    # Add the 500 hPa wind barbs
    plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
              to_np(u_500[::10, ::10]), to_np(v_500[::10, ::10]),
              transform=ccrs.PlateCarree(), length=6)
     
     # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_500))
    ax.set_ylim(cartopy_ylim(z_500))    
    # ax.gridlines() 
    # Titles
    plt.title('500 hPa Height (dm), Temperature (\u2103), Barbs (kt) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('500_hPa_Temperature__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
    #%% (3) 700 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])"
    """ 700 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])"""

    os.chdir(d01_basic_path + wind_700_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        

    # ---------------  Construct the interpolated map  ----------------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_700)      
    # Get the map projection information
    cart_proj = get_cartopy(z_700)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2)
    #  ----------------------------------------------------------------------
     # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    
    # Add the 700 hPa geopotential height contours
    levels = np.arange(z_700_min-1, z_700_max+1, 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_700),
                           levels=levels, colors="white",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
    
    # Add the wind speed contours
    levels = np.arange(wind_700_min, wind_700_max+1, 2.)
    wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_700),
                                 levels=levels,cmap=get_cmap("rainbow"),
                                 transform=ccrs.PlateCarree())
    cb = plt.colorbar(wspd_contours, ax=ax, orientation="vertical", shrink=0.7)
    cb.set_label('Wind Speed  700hPa (kt)',size=15)
    # Add the 700 hPa wind barbs
    plt.barbs(to_np(lons[::5,::5]), to_np(lats[::5,::5]),
              to_np(u_700[::5, ::5]), to_np(v_700[::5, ::5]),
              transform=ccrs.PlateCarree(), length=6)
     
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_700))
    ax.set_ylim(cartopy_ylim(z_700))    
    # ax.gridlines() 
    # Titles
    plt.title('700 hPa Height (dm), Wind Speed (kt), Barbs (kt)  \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('700_hPa_Wind_Speed__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
      #%%  (4) 700 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])
    """ 700 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])"""

    os.chdir(d01_basic_path + temp_700_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    # ---------------  Construct the interpolated map  ----------------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_700)      
    # Get the map projection information
    cart_proj = get_cartopy(z_700)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2)
    #  ----------------------------------------------------------------------
     # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())    
   
    # Now we build the contours for geopotential height at 700hPa
    levels = np.arange(z_700_min-1, z_700_max+1, 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_700),
                           levels=levels, colors="m",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Now we build the contour & contourf for temperature at 700hpa 
    temp_700_levels = np.arange(temp_700_min-1,temp_700_max+1,1) # specify the step of contours
    temp_700_contour  = ax.contour(to_np(lons), to_np(lats),temp_700,colors='k', levels=temp_700_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())       
    temp_700_contourf = ax.contourf(to_np(lons), to_np(lats), temp_700,levels=temp_700_levels,
                                cmap="jet", zorder=1, transform = ccrs.PlateCarree())    
       
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)    
    plt.clabel(temp_700_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(temp_700_contourf, shrink=0.5) # shrink= size of the colorbar
    cb.set_label('Temperature (\u2103)',size=15)    
    # add gridlines
    ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')
    # Add the 700 hPa wind barbs
    plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
              to_np(u_700[::10, ::10]), to_np(v_700[::10, ::10]),
              transform=ccrs.PlateCarree(), length=6)
     
     # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_700))
    ax.set_ylim(cartopy_ylim(z_700))    
    # ax.gridlines() 
    # Titles
    plt.title('700 hPa Height (dm), Temperature (\u2103), Barbs (kt)\n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('700_hPa_Temperature__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
       
    
    #%%  (5) 850 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])
    """ 850 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed])"""
 
    os.chdir(d01_basic_path + wind_850_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    # ---------------  Construct the interpolated map  ----------------------
    # Smooth the sea level pressure/geopotential height etc since they tend to be noisy near the
    # mountains
    z_850 = smooth2d(z_850, 3, cenweight=4)    
     # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_850)      
    # Get the map projection information
    cart_proj = get_cartopy(z_850)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2)  
    #  ----------------------------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='m' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())    
    
    # Add the 850 hPa geopotential height contours
    levels = np.arange(z_850_min-1, z_850_max+1, 1.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_850),
                           levels=levels, colors="white",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
    
    # Add the wind speed contours
    levels = np.arange(wind_850_min, wind_850_max+1, 2.)
    wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_850),
                                 levels=levels,cmap=get_cmap("rainbow"),
                                 transform=ccrs.PlateCarree())
    cb = plt.colorbar(wspd_contours, ax=ax, orientation="vertical", shrink=0.7)
    cb.set_label('Wind Speed  850hPa (kt)',size=15)
    # Add the 850 hPa wind barbs
    plt.barbs(to_np(lons[::5,::5]), to_np(lats[::5,::5]),
              to_np(u_850[::5, ::5]), to_np(v_850[::5, ::5]),
              transform=ccrs.PlateCarree(), length=6)
     
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_850))
    ax.set_ylim(cartopy_ylim(z_850))    
    # ax.gridlines() 
    # Titles
    plt.title('850 hPa Height (dm), Wind Speed (kt), Barbs (kt) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)     
    # Save the plot
    plt.savefig('850_hPa_Wind_Speed__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
    #%%  (6) 850 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])
    """ 850 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature])"""

    os.chdir(d01_basic_path + temp_850_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT 

    # ---------------  Construct the interpolated map  ----------------------
    # Smooth the sea level pressure/geopotential height etc since they tend to be noisy near the
    # mountains
    z_850 = smooth2d(z_850, 3, cenweight=4)    
     # Get the lat/lon coordinates
    lats, lons = latlon_coords(z_850)      
    # Get the map projection information
    cart_proj = get_cartopy(z_850)    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2)  
    #  ----------------------------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
   
    # Now we build the contours for geopotential height at 850hPa
    levels = np.arange(z_850_min-1, z_850_max+1, 1.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_850),
                           levels=levels, colors="m",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Now we build the contour & contourf for temperature at 850hpa 
    temp_850_levels = np.arange(temp_850_min-4,temp_850_max+4,1) # specify the step of contours
    temp_850_contour =ax.contour(to_np(lons), to_np(lats),temp_850,colors='k', levels=temp_850_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())       
    temp_850_contourf = ax.contourf(to_np(lons), to_np(lats), temp_850,levels=temp_850_levels,
                                cmap="jet", zorder=1, transform = ccrs.PlateCarree())    
       
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)    
    plt.clabel(temp_850_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(temp_850_contourf, shrink=0.5) # shrink= size of the colorbar
    cb.set_label('Temperature (\u2103)',size=15)    
    # add gridlines
    ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')
    # Add the 850 hPa wind barbs
    plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
              to_np(u_850[::10, ::10]), to_np(v_850[::10, ::10]),
              transform=ccrs.PlateCarree(), length=6)
     
     # Set the map bounds
    ax.set_xlim(cartopy_xlim(z_850))
    ax.set_ylim(cartopy_ylim(z_850))    
    # ax.gridlines() 
    # Titles
    plt.title('850 hPa Height (dm), Temperature (\u2103), Barbs (kt) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('850_hPa_Temperature__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()   
      
    #%% (7,8,9) CLOUD FRACTION FOR ALL  -->(low,medium,high) clouds
    """ CLOUD FRACTION FOR ALL  -->(low,medium,high) clouds """

    os.chdir(d01_basic_path + cloud_fraction_all_clouds_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    # Set Cloud fraction for (low,medium,high) clouds
    """ Cloud fraction plots for (low,medium,high) clouds   """    
    clfr_low_clouds    = cloud_fraction[0,i,:,:]*100  #  low_thresh:  300.0 m
    clfr_medium_clouds = cloud_fraction[1,i,:,:]*100  #  mid_thresh: 2000.0 m
    clfr_high_clouds   = cloud_fraction[2,i,:,:]*100 #  high_thresh: 6000.0 m
    max_clfr = 100
    min_clfr = 0
    
    #------------ What We Need For A Plot---------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z[i])      
    # Get the map projection information
    cart_proj = get_cartopy(z[i])    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    # ---------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
        
        
    # Add the Sea Level Pressure  contours
    levels = np.arange(slp_min-1, slp_max+1, 2.)
    slp_contours = plt.contour(to_np(lons), to_np(lats), to_np(slp[i]),
                           levels=levels, colors="orange",
                           transform=ccrs.PlateCarree())
    plt.clabel(slp_contours, inline=1, fontsize=10, fmt="%i")    
    
    """ CLOUD FRACTION FOR --> HIGH CLOUDS """
   
    clfr_high_clouds_levels = np.arange(min_clfr,max_clfr+1,5) # specify the step of contours
    clfr_high_clouds_contour =ax.contour(to_np(lons), to_np(lats),clfr_high_clouds,colors='k',
                                           levels=clfr_high_clouds_levels, linewidths=0.1,
                                           zorder=3,transform = ccrs.PlateCarree())    
    clfr_high_clouds_contourf = ax.contourf(to_np(lons), to_np(lats), clfr_high_clouds,alpha=0.5,
                                              levels=clfr_high_clouds_levels,cmap=high_cloud_colormap,
                                             zorder=1, transform = ccrs.PlateCarree())    
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(clfr_high_clouds_contour, inline=True, fmt='%1i', fontsize=12)    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(clfr_high_clouds_contourf, shrink=0.5, fraction=.02) # shrink= size of the colorbar
    cb.set_label('Cloud fraction HIGH CLOUDS (%)',size=15)      
    ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')  # add gridlines  
    
    """ CLOUD FRACTION FOR --> MEDIUM CLOUDS """

    clfr_medium_clouds_levels = np.arange(min_clfr,max_clfr+1,5) # specify the step of contours
    clfr_medium_clouds_contour =ax.contour(to_np(lons), to_np(lats),clfr_medium_clouds,colors='k',
                                           levels=clfr_medium_clouds_levels, linewidths=0.1,
                                           zorder=3,transform = ccrs.PlateCarree())    
    clfr_medium_clouds_contourf = ax.contourf(to_np(lons), to_np(lats), clfr_medium_clouds,alpha=0.5,
                                              levels=clfr_medium_clouds_levels,cmap=med_cloud_colormap,
                                             zorder=1, transform = ccrs.PlateCarree())    
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(clfr_medium_clouds_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb2  = plt.colorbar(clfr_medium_clouds_contourf, shrink=0.5, fraction=.02) # shrink= size of the colorbar
    cb2.set_label('Cloud fraction MEDIUM CLOUDS  (%)',size=15)    
    # ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')      # add gridlines
    
    """ CLOUD FRACTION FOR --> LOW CLOUDS """
  
    clfr_low_clouds_levels = np.arange(min_clfr,max_clfr+1,5) # specify the step of contours
    clfr_low_clouds_contour =ax.contour(to_np(lons), to_np(lats),clfr_low_clouds,colors='k', levels=clfr_low_clouds_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())    
    clfr_low_clouds_contourf = ax.contourf(to_np(lons), to_np(lats), clfr_low_clouds,alpha=0.5,
                                           levels=clfr_low_clouds_levels,cmap=low_cloud_colormap, zorder=1,
                                           transform = ccrs.PlateCarree())    
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(clfr_low_clouds_contour, inline=True, fmt='%1i', fontsize=12)    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb3  = plt.colorbar(clfr_low_clouds_contourf, shrink=0.5,fraction=.02) # shrink= size of the colorbar
    cb3.set_label('Cloud fraction LOW CLOUDS (%)',size=15)    
    # ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')   # add gridlines    
    # Titles
    plt.title('Cloud fraction for high (>6000 m), medium (2.000-6.000 m) and low (300-2.000 m) clouds \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot  
    plt.savefig('Cloud_fraction_ALL_clouds__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')    
    plt.close()     
              
    #%% (10) -------------  Construct a 10m wind map  ------------------ 
    """10m Wind Speed and Sea Level Pressure """

    os.chdir(d01_basic_path + wind_10_slp_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    #   ---------   Wind U10m ------------------------    
    # we give units to our data
    U10_wrf = wind_wrf_U10[i].values  * units('m/s')
    #   ---------   Wind V10  ------------------------
    # we give units to our data
    V10_wrf = wind_wrf_V10[i].values  * units('m/s')    
    # -----------  Wind speed construction of U10,V10 WRF --------
    wrf_WSpeed_10m = mpcalc.wind_speed(U10_wrf, V10_wrf)
    wrf_WSpeed_10m =(wrf_WSpeed_10m*1.9438444924406)* units('s/m*knots/knots') 
    wrf_WSpeed_10m 
    
    # ---------------- Construct wind plot ------------------------
     # Get the lat/lon coordinates
    lats, lons = latlon_coords(wind_wrf_U10[i])     
    # Get the map projection information
    cart_proj = get_cartopy(wind_wrf_U10[i])    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    # Download and add the coastlines
    ax.coastlines('10m', linewidth=2)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'),linewidth=2) 
    # -------------------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    
    # Add the Sea Level Pressure  contours
    levels = np.arange(slp_min-1, slp_max+1, 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(slp[i]),
                           levels=levels, colors="white",
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
    
    # Add the wind speed contours
    levels = np.arange(wind_10_slp_min, wind_10_slp_max+1, 2.)
    wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wrf_WSpeed_10m),
                                 levels=levels,cmap=get_cmap("rainbow"),
                                 transform=ccrs.PlateCarree())      
    cb = plt.colorbar(wspd_contours, ax=ax, orientation="vertical", shrink=0.7)
    cb.set_label('Wind Speed (kt)',size=15)
    
    U10_wrf = wind_wrf_U10.values
    V10_wrf = wind_wrf_V10.values
    # Add the 850 hPa wind barbs, only plotting every 5th data point.
    plt.barbs(to_np(lons[::5,::5]), to_np(lats[::5,::5]),
              to_np(U10_wrf[i][::5, ::5]), to_np(V10_wrf[i][::5, ::5]),
              transform=ccrs.PlateCarree(), length=6)   
    
   
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(wind_wrf_U10[i]))
    ax.set_ylim(cartopy_ylim(wind_wrf_U10[i]))    
    # ax.gridlines() 
    # Titles
    plt.title('Wind Speed at 10m Height (kt),Wind Barbs (kt) & Sea Level Pressure (hPa) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot   
    plt.savefig('Wind Speed_10m__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()

    #%% (11) AIR TEMPERATURE AT 2M
    """=========(AIR TEMPERATURE AT 2M PLOT ) ==========="""    
 
    os.chdir(d01_basic_path + t2m_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
        
    #------------ What We Need For A Plot---------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z[i])      
    # Get the map projection information
    cart_proj = get_cartopy(z[i])    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    # ---------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    lat =lat_made
    lon = lon_made
    wrf_temp_2m_levels = np.arange(wrf_temp_2m_min-1,wrf_temp_2m_max+3,1) # specify the step of contours
    wrf_temp_2m_contour =ax.contour(to_np(lons), to_np(lats),wrf_temp_2m[i],colors='k', levels=wrf_temp_2m_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())
    
    wrf_temp_2m_contour_extra = ax.contourf(to_np(lons), to_np(lats), wrf_temp_2m[i],levels=wrf_temp_2m_levels,
                                cmap="jet", zorder=1, transform = ccrs.PlateCarree())
    
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    plt.clabel(wrf_temp_2m_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(wrf_temp_2m_contour_extra, shrink=0.5) # shrink= size of the colorbar
    cb.set_label('Temperature (\u2103)',size=15)        
    ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')  # add gridlines  
    # Titles
    plt.title('Temperature (\u2103) at 2m \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot
    plt.savefig('Temperature_2m__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
    #%% (12)   """ PLOT THE SNOW HEIGHT """  perhaps it is right too as a snow parameter...
    # """ PLOT THE SNOW HEIGHT """
 
    # os.chdir(d01_basic_path + snow_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    # #------------ What We Need For A Plot---------------
    # # Get the lat/lon coordinates
    # lats, lons = latlon_coords(z[i])      
    # # Get the map projection information
    # cart_proj = get_cartopy(z[i])    
    # # Create the figure
    # fig = plt.figure(figsize=(20,10))
    # ax = plt.axes(projection=cart_proj)    
    # ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    # ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    # # ---------------------------------------------------
    # # Add cities annotations
    # for w in range(0, len(cities)):
    #     ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
    #     ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
        
    # lat =lat_made
    # lon = lon_made
    # wrf_snowh_levels = np.arange(snow_min,snow_max+3,1) # specify the step of contours
    # wrf_snowh_contour =ax.contour(to_np(lons), to_np(lats),wrf_snowh[i],colors='blue', levels=wrf_snowh_levels, 
    #                           linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())
    # wrf_snowh_contourf = ax.contourf(to_np(lons), to_np(lats), wrf_snowh[i],levels=wrf_snowh_levels,
    #                             cmap="Blues", zorder=1, transform = ccrs.PlateCarree())    
    
    # # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(wrf_snowh_contour, inline=True, fmt='%1i', fontsize=12)
    
    # #Create a wrf_variable colorbar and shrink it down a bit.
    # cb  = plt.colorbar(wrf_snowh_contourf, shrink=0.5) # shrink= size of the colorbar
    # cb.set_label('Snow Height (cm)',size=15)        
    # ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--') # add gridlines    
    
    # plt.title('Snow Height (cm) \n\n\n ',size=15)
    # plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    # plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    
    
    # # plt.show()   
    # plt.savefig('Snow Height__'+filenames[0]+' Valid:_'+str(wrf_total_precip[i].Time.values)[:19]+'__.png',
    #             dpi=150,bbox_inches='tight')
    # plt.close()  

    ##########################################################################
     #%% (12bbbb)   """ PLOT THE RIGHT SNOW HEIGHT PARAMETER"""
    """ PLOT THE RIGHT SNOW HEIGHT PARAMETER 
    
    *** WE REPLACE THIS PARAMETER WITH THE Current Amount of Rain & Snow water equivalent (mm)"""
 
    # os.chdir(d01_basic_path + snow_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT   
         
    # #------------ What We Need For A Plot---------------
    # # Get the lat/lon coordinates
    # lats, lons = latlon_coords(z[i])      
    # # Get the map projection information
    # cart_proj = get_cartopy(z[i])    
    # # Create the figure
    # fig = plt.figure(figsize=(20,10))
    # ax = plt.axes(projection=cart_proj)    
    # ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    # ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    # # ---------------------------------------------------
    # # Add cities annotations
    # for w in range(0, len(cities)):
    #     ax.plot(city_lon[w], city_lat[w], 'bo',color='grey' ,markersize=5, transform=ccrs.Geodetic())
    #     ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
        
    # lat =lat_made
    # lon = lon_made   
    # wrf_snowh_levels = np.arange(snow_min,snow_max,3) # specify the step of contours
    # # wrf_snowh_levels = np.arange(snow_min,30,1) # specify the step of contours
    # # nn = len(snow_colors)              #number of steps
    # # wrf_snowh_levels = 1.2**np.arange(0, nn, 1)-1 # specify the step of contours

    # wrf_snowh_contour =ax.contour(to_np(lons), to_np(lats),wrf_snow_wat_equiv[i],colors='blue',
    #                               levels=wrf_snowh_levels,linewidths=0.1, zorder=3,
    #                               transform = ccrs.PlateCarree())
    # wrf_snowh_contourf = ax.contourf(to_np(lons), to_np(lats), wrf_snow_wat_equiv[i],levels=wrf_snowh_levels,
    #                             cmap=snow_colormap, zorder=1, transform = ccrs.PlateCarree())    
    
    # # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(wrf_snowh_contour, inline=True, fmt='%1i', fontsize=12)
    
    # #Create a wrf_variable colorbar and shrink it down a bit.
    # cb  = plt.colorbar(wrf_snowh_contourf, shrink=0.5) # shrink= size of the colorbar
    # cb.set_label('Snow water equivalent (mm)',size=15)     
      
    # ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--') # add gridlines    
    
    # plt.title('Snow Accumulation(mm) \n\n\n ',size=15)
    # plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    # plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # # Save the plot   
    # plt.savefig('Snow_Accumulation__'+filenames[0]+' Valid:_'+str(wrf_total_precip[i].Time.values)[:19]+'__.png',
    #             dpi=150,bbox_inches='tight')
    # plt.close()  

    # ##########################################################################
    #%%   (13)   """ PLOT THE RELATIVE HUMIDITY AT 2M """
    """ RELATIVE HUMIDITY AT 2M """
 
    os.chdir(d01_basic_path + rh2_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT        
    
    #------------ What We Need For A Plot---------------
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(z[i])      
    # Get the map projection information
    cart_proj = get_cartopy(z[i])    
    # Create the figure
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection=cart_proj)    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    # ---------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='red' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    lat =lat_made
    lon = lon_made
    wrf_rh2_levels = np.arange(rh2_min,rh2_max+1,10) # specify the step of contours
    wrf_rh2_contour =ax.contour(to_np(lons), to_np(lats),wrf_rh2[i],colors='k', levels=wrf_rh2_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())
    wrf_rh2_contourf = ax.contourf(to_np(lons), to_np(lats), wrf_rh2[i],levels=wrf_rh2_levels,
                                cmap="cividis", zorder=1, transform = ccrs.PlateCarree())    
    
    # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    plt.clabel(wrf_rh2_contour, inline=True, fmt='%1i', fontsize=12)
    
    #Create a wrf_variable colorbar and shrink it down a bit.
    cb  = plt.colorbar(wrf_rh2_contourf, shrink=0.5) # shrink= size of the colorbar
    cb.set_label('Relative Humidity (%)',size=15)    
    ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--')    # add gridlines 
    # Titles
    plt.title('Relative Humidity (%) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot   
    plt.savefig('Relative Humidity__'+filenames[0]+' Valid:_'+str(wrf_total_precip[i].Time.values)[:19]+'__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
    
   #%% (14)  Current Amount of Rain & Snow water equivalent plot + Total Precipitation/Snow plots
   
    """ Current Amount of Rain & Snow water equivalent + Total Precipitation/Snow plots"""
     
    os.chdir(d01_basic_path + rain_dir) # SOS:WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT   
     
    #------------ What We Need For A Plot---------------
    # Get the lat/lon coordinates    
    lats, lons = latlon_coords(z[i])       
    # Get the map projection information
    cart_proj = get_cartopy(z[i])    
    # Create the figure
    fig = plt.figure(figsize=(22,12))
    ax = plt.axes(projection=cart_proj)    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))  
    # ax.add_feature(cfeature.LAND.with_scale('10m'), edgecolor='gray')
    # ax.add_feature(cfeature.OCEAN.with_scale('10m')) ## here it stucks
    # ax.add_feature(cfeature.LAKES.with_scale('10m'), alpha=0.5)
    # ax.add_feature(cfeature.RIVERS)          
    # ---------------------------------------------------
    # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='red' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    lat =lat_made
    lon = lon_made      
    
    '''Snow variable construction '''
    snow = wrf_snow_wat_equiv[i] # total ppn / mm
    currsnow = np.array(snow.shape)
    if i == 0: # initial amount
        currsnow = snow
    else: # current amount
        prev =wrf_snow_wat_equiv[i-1]
        currsnow = snow-prev    
    # Snow Colormap   
    cmap_snow = snow_colormap
    # cmap_snow.set_over('white')
    # cmap_snow.set_under('white')    
    wrf_snowh_levels = np.arange(0,12,1) # specify the step of contours    
    labels = wrf_snowh_levels[1:-1]
    norm = matplotlib.colors.BoundaryNorm(wrf_snowh_levels, cmap_snow.N)      
    wrf_snowh_contour =ax.contour(to_np(lons), to_np(lats),currsnow,colors='k',
                                  levels=wrf_snowh_levels,linewidths=0.1, zorder=3,
                                  transform = ccrs.PlateCarree())
    wrf_snowh_contourf = ax.contourf(to_np(lons), to_np(lats), currsnow,levels=wrf_snowh_levels,
                                cmap=snow_colormap, zorder=1, alpha=0.5,
                                transform = ccrs.PlateCarree(),norm=norm, extend='both')         
    # SNOW ''' Plot contour labels for the SNOW vatiable, leaving a break in the contours for the text (inline=True) '''
    plt.clabel(wrf_snowh_contour, inline=True, fmt='%1i', fontsize=12)    
    #Create the snow colorbar and shrink it down a bit.
    cb  = plt.colorbar(wrf_snowh_contourf, shrink=0.5,fraction=.02) # shrink= size of the colorbar
    cb.set_label('Snow water equivalent (mm)',size=15)  
        
    '''Current Rain variable construction ''' 
     #  Define the current amount of rain      
    ppn = wrf_total_precip[i] # total ppn / mm
    currppn = np.array(ppn.shape)
    if i == 0: # initial amount
        currppn = ppn
    else: # current amount
        prev =wrf_total_precip[i-1]
        currppn = ppn-prev        
     # Rain Colormap   
    cmap_rain = precip_colormap
    # cmap_rain.set_over('white')
    # cmap_rain.set_under('white')     
    nn = 18                     #number of steps
    wrf_total_precip_levels = 1.2**np.arange(0, nn, 1)-1 # specify the step of contours
    
    labels = wrf_total_precip_levels[1:-1]
    norm_rain = matplotlib.colors.BoundaryNorm(wrf_total_precip_levels, cmap_rain.N)
    # rain_step = rain_max/20
    # wrf_total_precip_levels = np.arange(rain_min,rain_max,rain_step) # specify the step of contours
    # wrf_total_precip_contour =ax.contour(to_np(lons), to_np(lats), currppn, colors='k', levels=wrf_total_precip_levels, 
    #                           linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())
    wrf_curr_precip_contourf = ax.contourf(to_np(lons), to_np(lats), currppn,levels=wrf_total_precip_levels,                                            
                                cmap=precip_colormap, zorder=1,alpha=0.6, transform = ccrs.PlateCarree(),
                                norm=norm_rain, extend='both')      
    #RAIN ''' Plot contour labels for the RAIN, leaving a break in the contours for the text (inline=True) '''
    # plt.clabel(wrf_total_precip_contour, inline=True, fmt='%1i', fontsize=12)    
    #Create THE RAIN colorbar and shrink it down a bit.
    cb2  = plt.colorbar(wrf_curr_precip_contourf, shrink=0.5,fraction=.02) # shrink= size of the colorbar
    cb2.set_label('Current Amount of Rain (mm)',size=15)        
    # ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--') # add gridlines     
    
    # Titles    
    plt.title('Current Amount of Rain & Snow water equivalent (mm) \n\n\n ',size=15)
    plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    # Save the plot  
    plt.savefig('Current_Rain_&_Snow__'+filenames[0][:24]+'UTC'+'__Valid:_'+str(wrf_total_precip[i].Time.values)[:13]+'UTC__.png',
                dpi=150,bbox_inches='tight')
    plt.close()
    
    
    '''Total Precipitation variable construction FOR THE WHOLE RUN 
        ( plot total precip accumulation)''' 
                
    #------------ What We Need For A Plot---------------
    # Get the lat/lon coordinates    
    lats, lons = latlon_coords(z[i])       
    # Get the map projection information
    cart_proj = get_cartopy(z[i])    
    # Create the figure
    fig = plt.figure(figsize=(22,12))
    ax = plt.axes(projection=cart_proj)    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))    
    # ---------------------------------------------------
     # Add cities annotations
    for w in range(0, len(cities)):
        ax.plot(city_lon[w], city_lat[w], 'bo',color='red' ,markersize=5, transform=ccrs.Geodetic())
        ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
    lat =lat_made
    lon = lon_made 
    wrf_total_precip_levels = 1.32**np.arange(0, nn, 1)-1 # specify the step of contours

    # wrf_total_precip[i] =wrf_total_precip[i].copy()
    wrf_total_precip_contour =ax.contour(to_np(lons), to_np(lats),wrf_total_precip[i],colors='k', levels=wrf_total_precip_levels, 
                              linewidths=0.1, zorder=3,transform = ccrs.PlateCarree())
    wrf_total_precip_contourf = ax.contourf(to_np(lons), to_np(lats), wrf_total_precip[i],levels=wrf_total_precip_levels,
                                cmap=precip_colormap, zorder=1, transform = ccrs.PlateCarree()) 
# Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
plt.clabel(wrf_total_precip_contour, inline=True, fmt='%1i', fontsize=12)
    
#Create a wrf_variable colorbar and shrink it down a bit.
cb3  = plt.colorbar(wrf_total_precip_contourf, shrink=0.5) # shrink= size of the colorbar
cb3.set_label('Total precipitation (mm)',size=15)        
ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--') # add gridlines  

# Titles   
plt.title('Total precipitation (mm) for the run \n\n\n ',size=15)
plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
# Save the plot
plt.savefig('Total_precipitation For the run__'+filenames[0][:24]+'UTC',dpi=150,bbox_inches='tight')
plt.close()

'''Total Snow accumulation variable construction FOR THE WHOLE RUN 
        ( plot total snow accumulation)''' 
#------------ What We Need For A Plot---------------
# Get the lat/lon coordinates    
lats, lons = latlon_coords(z[i])       
# Get the map projection information
cart_proj = get_cartopy(z[i])    
# Create the figure
fig = plt.figure(figsize=(22,12))
ax = plt.axes(projection=cart_proj)    
ax.add_feature(cfeature.BORDERS.with_scale('10m'))
ax.add_feature(cfeature.COASTLINE.with_scale('10m'))    
# ---------------------------------------------------
 # Add cities annotations
for w in range(0, len(cities)):
    ax.plot(city_lon[w], city_lat[w], 'bo',color='red' ,markersize=5, transform=ccrs.Geodetic())
    ax.text(city_lon[w], city_lat[w], cities[w], transform=ccrs.Geodetic())
lat =lat_made
lon = lon_made 
wrf_snowh_levels = np.arange(1,32,2) # specify the step of contours    
for i in tqdm(range(12,len(time))):    
    wrf_snowh_contour =ax.contour(to_np(lons), to_np(lats),wrf_snow_wat_equiv[i],cmap=snow_colormap,
                                  levels=wrf_snowh_levels,linewidths=0.1, zorder=3,
                                  transform = ccrs.PlateCarree())
    wrf_snowh_contourf = ax.contourf(to_np(lons), to_np(lats), wrf_snow_wat_equiv[i],levels=wrf_snowh_levels,
                                cmap=snow_colormap, zorder=1, alpha=0.5,
                                transform = ccrs.PlateCarree(),norm=norm, extend='both')         
# SNOW ''' Plot contour labels for the SNOW vatiable, leaving a break in the contours for the text (inline=True) '''
plt.clabel(wrf_snowh_contour, inline=True, fmt='%1i', fontsize=12)    
#Create the snow colorbar and shrink it down a bit.
cb4  = plt.colorbar(wrf_snowh_contourf, shrink=0.5,fraction=.02) # shrink= size of the colorbar
cb4.set_label(' Total Snow water equivalent (mm)',size=15)     
     
# ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='ivory', alpha=0.1, linestyle='--') # add gridlines  

# Titles   
plt.title('Total snow water equivalent (mm) for the run \n\n\n ',size=15)
plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
# Save the plot
plt.savefig('Total_Snow For the run__'+filenames[0][:24]+'UTC',dpi=150,bbox_inches='tight')
plt.close()
    

    #%%  #%% (15) """ MAXIMUM REFLECTIVITY """
    # """ MAXIMUM REFLECTIVITY """
    # os.chdir(d01_basic_path + Max_refl) # SOS: WE CHANGE DIRECTORY TO SAVE THE PLOT WHERE WE WANT
    #  # Get the lat/lon coordinates
    # lats, lons = latlon_coords(z[i])     
    # # Get the map projection information
    # cart_proj = get_cartopy(z[i])
    #
    # # Create the figure
    # fig = plt.figure(figsize=(20,10))
    # ax = plt.axes(projection=cart_proj)    
    # ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    # ax.add_feature(cfeature.COASTLINE.with_scale('10m'))    
    # max_reflect_highest = 70
    # min_reflect_lowest = 0
    # max_reflect_levels  = np.arange(Max_refl_max-1,Max_refl_max+1,5) # specify the step of contours
    # max_reflect_contour = ax.contour(to_np(lons), to_np(lats),max_reflect[i],colors='k',
    #                                        levels=max_reflect_levels, linewidths=0.1,
    #                                        zorder=3,transform = ccrs.PlateCarree())
    #
    # max_reflect_contourf = ax.contourf(to_np(lons), to_np(lats), max_reflect[i],
    #                                           levels=max_reflect_levels,cmap="cool",
    #                                          zorder=1, transform = ccrs.PlateCarree())
    #
    # # Plot contour labels for the wrf_variable, leaving a break in the contours for the text (inline=True)
    # plt.clabel(max_reflect_contour, inline=True, fmt='%1i', fontsize=12)
    #
    # #Create a wrf_variable colorbar and shrink it down a bit.
    # cb  = plt.colorbar(max_reflect_contourf, shrink=0.5) # shrink= size of the colorbar
    # cb.set_label('Maximum Reflectivity (dbz)',size=15)       
    # ax.gridlines(linewidth=1, color='ivory', alpha=0.1, linestyle='--')   # add gridlines    
    #
    # plt.title('Maximum Reflectivity (dbz) \n\n\n ',size=15)
    # plt.title('Valid: '+str(wrf_total_precip[i].Time.values)[:13]+ 'UTC \n', loc='left',size=15)
    # plt.title('WRF-ARW (V4.1.5) '+wrf_grid_resolution+ '\n Init: '+filenames[0][:24]+'UTC \n', loc='right',size=12)
    #
    # # plt.show()   
    # plt.savefig('MAXIMUM REFLECTIVITY__'+filenames[0]+'  Valid: '+str(wrf_total_precip[i].Time.values)[:19]+'__.png',
    #             dpi=150,bbox_inches='tight')
    # plt.close()    
#%%
# !!! SOS:NOW WE RETURN TO THE DIRECTORY WITH THE WRF OUTPUT FILES 
os.chdir(pwd_directory)     
    
end_time = datetime.now()
print ("End_time = ", end_time)

dif = end_time - start_time  
print (" (d01) Visualitation:> Total execution time = ",dif)
            