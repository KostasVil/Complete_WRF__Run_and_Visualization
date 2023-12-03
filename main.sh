#!/bin/bash

STARTYEAR=` date +'%Y'`
STARTMONTH=`date +'%m'`
STARTDAY=`date +'%d'`
STARTHOUR=00
RUNLENGTH=72


# Environment:
WRFROOT=/path_to_installed_WRF
WPSDIR=$WRFROOT/WPS-4.1
WRFDIR=$WRFROOT/WRF-4.1.5
SCRIPTSDIR=$WRFROOT/scripts_path
CONFIGDIR=$WRFROOT/config


# Cleaning:
echo "                               "
echo "   -- Cleaning:--> Delete old files startin from:___   FILE,met_em,ungrib.log,  ---------  "
echo "   ----- metgrid.log, GRIBFILE, namelist.wps ___ IN THE  WPS DIRECTORY    --------     "
echo "                               "
rm -rf $WPSDIR/FILE:*
rm -rf $WPSDIR/met_em*
rm -rf $WPSDIR/ungrib.log
rm -rf $WPSDIR/metgrid.log
rm -rf $WPSDIR/GRIBFILE.*
rm -rf $WPSDIR/namelist.wps


echo "                               "
echo "   -- Cleaning:--> Delete old files startin from:___   wrfout, wrfbdy, wrfinput,  ---------  "
echo "   ----- rsl, met_em, namelist.input  ____ IN THE  WRF DIRECTORY    --------     "
echo "                               "
rm -rf $WRFDIR/run/wrfout*
rm -rf $WRFDIR/run/wrfbdy*
rm -rf $WRFDIR/run/wrfinput*
rm -rf $WRFDIR/run/rsl.*
rm -rf $WRFDIR/run/met_em*
rm -rf $WRFDIR/run/namelist.input

# As we set the namelist.wps & namelist.input  files in WPS & WRF respectively
# we run the geogrid.exe in WPS file
cp $CONFIGDIR/namelist.wps $WPSDIR/namelist.wps
cd $WPSDIR
echo "                               "
echo "   ------------- geogrid.exe IS RUNNING  ---------------------  "
echo "                               "
./geogrid.exe
echo "                               "

# Downloading GFS
echo "                               "
echo "  -----------     Downloading GFS   ------------------------"
bash $SCRIPTSDIR/download_gfs.sh $STARTYEAR $STARTMONTH $STARTDAY $STARTHOUR

# setting namelist.wps
echo "                               "
echo "   ----------   setting namelist.wps    ---------------"
echo "                               "
cp $CONFIGDIR/namelist.wps $WPSDIR/namelist.wps
sed -i "s/STARTYEAR/$STARTYEAR/g" $WPSDIR/namelist.wps
sed -i "s/STARTMONTH/$STARTMONTH/g" $WPSDIR/namelist.wps
sed -i "s/STARTDAY/$STARTDAY/g" $WPSDIR/namelist.wps
sed -i "s/STARTHOUR/$STARTHOUR/g" $WPSDIR/namelist.wps
# Here we copy namelist.wps from CONFIGDIR to WPSDIR and 
# we get into WPSDIR and we change START**** with the value of $START***  to namelist.wps
#  "sed" is the command for replacing every START**** with the  value of $START**

# we use -u flag in order to take the hour in UTC
START="$STARTYEAR-$STARTMONTH-$STARTDAY $STARTHOUR:00:00Z"
END=`date -u --date="$START + $RUNLENGTH hours" +'%Y-%m-%d %H:00:00Z'`
ENDYEAR=`date -u --date="$END" +'%Y'`
ENDMONTH=`date -u --date="$END" +'%m'`
ENDDAY=`date -u --date="$END" +'%d'`
ENDHOUR=`date -u --date="$END" +'%H'`

sed -i "s/ENDYEAR/$ENDYEAR/g" $WPSDIR/namelist.wps
sed -i "s/ENDMONTH/$ENDMONTH/g" $WPSDIR/namelist.wps
sed -i "s/ENDDAY/$ENDDAY/g" $WPSDIR/namelist.wps
sed -i "s/ENDHOUR/$ENDHOUR/g" $WPSDIR/namelist.wps


# WPS programs
echo "                               "
echo "   ------------- running WPS programs  ---------------------  "
echo "                               "
cd $WPSDIR
./link_grib.csh /home/kostas/WRF/GFS/


echo "   ------------- ungrib.exe IS RUNNING  ---------------------  "
echo "                               "
./ungrib.exe
echo "                               "
echo "   ------------- metgrid.exe IS RUNNING  ---------------------  "
./metgrid.exe
echo "                               "


# setting namelist.input
echo "                               "
echo "  -------------  setting namelist.input   -------------------- "
echo "                               "
cp $CONFIGDIR/namelist.input $WRFDIR/run/namelist.input

sed -i "s/STARTYEAR/$STARTYEAR/g" $WRFDIR/run/namelist.input
sed -i "s/STARTMONTH/$STARTMONTH/g" $WRFDIR/run/namelist.input
sed -i "s/STARTDAY/$STARTDAY/g" $WRFDIR/run/namelist.input
sed -i "s/STARTHOUR/$STARTHOUR/g" $WRFDIR/run/namelist.input

sed -i "s/ENDYEAR/$ENDYEAR/g" $WRFDIR/run/namelist.input
sed -i "s/ENDMONTH/$ENDMONTH/g" $WRFDIR/run/namelist.input
sed -i "s/ENDDAY/$ENDDAY/g" $WRFDIR/run/namelist.input
sed -i "s/ENDHOUR/$ENDHOUR/g" $WRFDIR/run/namelist.input


# WRF programs
echo "                               "
echo "  ---------------  running WRF programs (please wait!!!)   ------------------"
echo "                               "
cd $WRFDIR/run
ln -s $WPSDIR/met_em* .
echo " -------- Run the real.exe file  ------------              "
mpirun -n 3 ./real.exe
echo "                               "
echo " -------- Run the wrf.exe file  ------------              "
mpirun -n 10 ./wrf.exe
echo "                               "
echo "running WRF programs DONE"
echo "-------------------------------------------------- "



# PYTHON Visualization
cd $SCRIPTSDIR
./wrf_visualization_ALL_DOMAINS.sh


