#! /bin/sh

# Using proj4 from given parameters


XP=418.25            # X grid coordinate of north pole
YP=257.25            # Y grid coordinate of north pole
DX=10000.0           # Grid size [m] 
YLON=58.0            # Longitude parallel to Y-axis

LAT_TS=60.0          # Latitude of true scale
RADIUS=6371000       # Earth radius [m] (spherical case)

PROJ=/usr/bin/proj   # proj command

INFILE=lonlat.dat    # file with some longitude latitude pairs

# Multiply XP and YP with DX
XPDX=$(echo "($XP*$DX)" | bc)
YPDX=$(echo "($YP*$DX)" | bc)

# ==============
# Sphere case
# ==============

echo " --- Sphere ---"

$PROJ -m 1:$DX +proj=stere +R=$RADIUS +lat_0=90.0 +lat_ts=$LAT_TS +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON $INFILE

# =============
# WGS84 case
# =============

echo " --- WGS84 ---"

$PROJ -m 1:$DX +proj=stere  +ellps=WGS84 +lat_0=90.0 +lat_ts=$LAT_TS +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON $INFILE


