#!/bin/bash


# Environments:
WRFOUTDIR=/path_to_run_directory_of_WRF  
SCRIPTSDIR=/path_to_scripts_directory  


# PYTHON Visualization
echo "                               "
echo "        --------------------------------------------------         "
echo "  ------- running WRF Python Visualization for ALL DOMAINS!!! ------------     "
echo "        --------------------------------------------------         "

# DOMAIN 01
cd $SCRIPTSDIR
echo "                               "
echo "  -----   running WRF d01 Python Visualization    -------  "
echo "                               "
cp $SCRIPTSDIR/d01_maps_15_variables.py $WRFOUTDIR/d01_maps_15_variables.py
cd $WRFOUTDIR
./d01_maps_15_variables.py
echo "                               "
cd $SCRIPTSDIR
./d01_create_gif.py
rm -rf $WRFOUTDIR/d01_maps_15_variables.py
echo "                               "
echo " -------    WRF d01 Python Visualization DONE -----------"
echo "     --------------------------------------------------     "

# DOMAIN 02
cd $SCRIPTSDIR
echo "                               "
echo "  -----   running WRF d02 Python Visualization    -------  "
echo "                               "
cp $SCRIPTSDIR/d02_maps_15_variables.py $WRFOUTDIR/d02_maps_15_variables.py
cd $WRFOUTDIR
./d02_maps_15_variables.py
echo "                               "
cd $SCRIPTSDIR
./d02_create_gif.py
rm -rf $WRFOUTDIR/d02_maps_15_variables.py
echo "                               "
echo " -------    WRF d02 Python Visualization DONE -----------"
echo "     --------------------------------------------------     "

# DOMAIN 03
cd $SCRIPTSDIR
echo "                               "
echo "  -----   running WRF d03 Python Visualization    -------  "
echo "                               "
cp $SCRIPTSDIR/d03_maps_15_variables.py $WRFOUTDIR/d03_maps_15_variables.py
cd $WRFOUTDIR
./d03_maps_15_variables.py
echo "                               "
cd $SCRIPTSDIR
./d03_create_gif.py
rm -rf $WRFOUTDIR/d03_maps_15_variables.py
echo "                               "
echo " -------    WRF d03 Python Visualization DONE -----------"
echo "     --------------------------------------------------     "

# END OF WRF VISUALIZATION 
echo "        --------------------------------------------------         "
echo "  ------  WRF Python Visualization DONE for ALL DOMAINS  --------   "   
echo "        --------------------------------------------------         "




