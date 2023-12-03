#!/bin/bash

inputdir=/path_to_GFS_data
rm -rf $inputdir
mkdir $inputdir

# Get the current year,month, day 
year=` date +'%Y'`
month=`date +'%m'`
day=`date +'%d'`
cycle=00


# for 72hour run set-->   i<=72
for ((i=000; i<=72; i+=3))
do
    ftime=`printf "%03d\n" "${i}"`

    server=https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod
    directory=gfs.${year}${month}${day}/${cycle}/atmos
    file=gfs.t${cycle}z.pgrb2.0p50.f${ftime}

    url=${server}/${directory}/${file}

    echo $url

    wget -O ${inputdir}/${file} ${url}

done 
