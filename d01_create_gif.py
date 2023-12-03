#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:16:18 2021

@author: kostas
"""

"""
# IF ANYONE NEEDS THE PRESENT CODE TO ADJUST IT IN HIS MACHINE HE ONLY HAS TO CHANGE
# THE #!!! d01_basic_path LINE (line 21)   WITH THE RIGHT PATH.
# PERHAPS HE WILL NEED TO CHANGE THE REST OF THE PATHS TOO.
# THEN THE CODE WILL USE THE .png FILES THAT THEY EXIST IN THESE PATHS AND IT WILL 
# CONSTRUST A GIF FOR EACH ONE OF THEM.
"""
#%%  Set paths for every variable
print("\n  Set paths for every variable \n")



d01_basic_path ="/path_to_WRF_visualization/d01_temp_maps"  #!!!

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
# snow_dir        = "/snow"
rh2_dir         = "/rh2"
rain_dir        = "/rain"
Max_refl_dir    = "/Max_refl"

variable_directory = [wind_500_dir, temp_500_dir, wind_700_dir, temp_700_dir,
                      wind_850_dir, temp_850_dir, cloud_fraction_all_clouds_dir,
                      wind_10_slp_dir, t2m_dir,rh2_dir, rain_dir,
                      Max_refl_dir]




"""  This code creates GIFS from the wrf variables we plotted to maps"""
# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif

from datetime import datetime
import glob
from PIL import Image

start_time = datetime.now()
print ("\n GIFS for d01 starts at: ",start_time)

print(" \n --------    GIF CREATION FOR d01 DOMAIN   --------   \n")


"""
# BE CAREFUL: IF THE CODE DOESN'T FIND AT LEAST 1 .png FILE IN ANY DIRECTORY OF THE ABOVE
# IT WILL GIVE AN ERROR!
"""

#%%   (1)
# """ Maximum Reflectivity GIF """
# # ----------------------------------------------------------------------------
# # filepaths
# fp_in  = d01_basic_path + Max_refl_dir+"/MAXIMUM REFLECTIVITY__wrfout*.png"
# fp_out = d01_basic_path + Max_refl_dir+"/MAXIMUM_REFLECTIVITY_image.gif"

# img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
# img.save(fp=fp_out, format='GIF', append_images=imgs,
#           save_all=True, duration=1600, loop=0)

#%%   (2)
""" 500 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + wind_500_dir+"/500_hPa_Wind_Speed__wrfout*.png"
fp_out = d01_basic_path + wind_500_dir+"/500_hPa_Wind_Speed__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (3)
""" 500 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + temp_500_dir+"/500_hPa_Temperature__wrfout*.png"
fp_out = d01_basic_path + temp_500_dir+"/500_hPa_Temperature__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (4)
""" 700 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + wind_700_dir+"/700_hPa_Wind_Speed__wrfout*.png"
fp_out = d01_basic_path + wind_700_dir+"/700_hPa_Wind_Speed__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (5)
""" 700 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + temp_700_dir+"/700_hPa_Temperature__wrfout*.png"
fp_out = d01_basic_path + temp_700_dir+"/700_hPa_Temperature__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)

#%%   (6)
""" 850 hPa LEVEL  (wind barbs, [contour>dm] & [contourf>wind speed]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + wind_850_dir+"/850_hPa_Wind_Speed__wrfout*.png"
fp_out = d01_basic_path + wind_850_dir+"/850_hPa_Wind_Speed__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (7)
""" 850 hPa LEVEL  (wind speed, [contour>dm] & contourf> temperature]) GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + temp_850_dir+"/850_hPa_Temperature__wrfout*.png"
fp_out = d01_basic_path + temp_850_dir+"/850_hPa_Temperature__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (8)
# """ CLOUD FRACTION FOR --> LOW CLOUDS GIF"""
# # ----------------------------------------------------------------------------
# # # filepaths
# # fp_in  = d01_basic_path + cloud_fraction_low_clouds_dir+"/Cloud_fraction_low_clouds__wrfout*.png"
# # fp_out = d01_basic_path + cloud_fraction_low_clouds_dir+"/Cloud_fraction_low_clouds__image.gif"

# # img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
# # img.save(fp=fp_out, format='GIF', append_images=imgs,
# #          save_all=True, duration=1600, loop=0)

# #%%   (9)
# """ CLOUD FRACTION FOR --> MEDIUM CLOUDS GIF"""
# # ----------------------------------------------------------------------------
# # # filepaths
# # fp_in  = d01_basic_path + cloud_fraction_medium_clouds_dir+"/Cloud_fraction_medium_clouds__wrfout*.png"
# # fp_out = d01_basic_path + cloud_fraction_medium_clouds_dir+"/Cloud_fraction_medium_clouds__image.gif"

# # img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
# # img.save(fp=fp_out, format='GIF', append_images=imgs,
# #          save_all=True, duration=1600, loop=0)

#%%   (8,9,10)
""" CLOUD FRACTION FOR ALL  -->(low,medium,high) clouds """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + cloud_fraction_all_clouds_dir+"/Cloud_fraction_ALL_clouds__wrfout*.png"
fp_out = d01_basic_path + cloud_fraction_all_clouds_dir+"/Cloud_fraction_ALL_clouds__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (11)
"""10m Wind Speed and Sea Level Pressure GIF"""
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + wind_10_slp_dir+"/Wind Speed_10m__wrfout*.png"
fp_out = d01_basic_path + wind_10_slp_dir+"/Wind Speed_10m__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%   (12)
""" Temperature at 2m GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + t2m_dir+"/Temperature_2m__wrfout*.png"
fp_out = d01_basic_path + t2m_dir+"/Temperature_2m_image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


# %%   (13)

""" Snow Height GIF """
# ----------------------------------------------------------------------------
# # filepaths
# fp_in  = d01_basic_path + snow_dir+"/Snow Height__wrfout*.png"
# fp_out = d01_basic_path + snow_dir+"/Snow_Height__image.gif"

# img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
# img.save(fp=fp_out, format='GIF', append_images=imgs,
#           save_all=True, duration=1600, loop=0)


#%%  (14)

""" RELATIVE HUMIDITY AT 2M  GIF"""
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + rh2_dir+"/Relative Humidity__wrfout*.png"
fp_out = d01_basic_path + rh2_dir+"/Relative Humidity__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
          save_all=True, duration=1600, loop=0)


#%%   (15)

""" Current Amount of Rain & Snow water equivalent GIF """
# ----------------------------------------------------------------------------
# filepaths
fp_in  = d01_basic_path + rain_dir+"/Current_Rain_&_Snow__wrfout*.png"
fp_out = d01_basic_path + rain_dir+"/Current_Rain_&_Snow__image.gif"

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=1600, loop=0)


#%%

end_time = datetime.now()
print ("\n (d01 GIFS) >> End_time = ", end_time)

dif = end_time - start_time  
print ("\n (d01 GIFS) >> Total execution time = ",dif)

















