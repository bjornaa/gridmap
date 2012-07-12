#! /bin/sh



# Using proj4 from given parameters


XP=418.25            # X grid coordinate of north pole
YP=257.25            # Y grid coordinate of north pole
DX=10000.0           # Grid spacing [m] 
YLON=58.0            # Longitude parallel to Y-axis

LAT_TS=60.0          # Latitude of true scale
RADIUS=6371000       # Earth radius [m] (spherical case)

PROJ=/usr/bin/proj   # proj command

# Multiply XP and YP with DX
XPDX=$(echo "($XP*$DX)" | bc)
YPDX=$(echo "($YP*$DX)" | bc)

LON=2
LAT=66
echo "Forward: lon = $LON, lat=$LAT"
$PROJ -m 1:$DX -f %10.6f +proj=stere +R=$RADIUS +lat_0=90.0 +lat_ts=$LAT_TS +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON <<EOF
  2.0  66.0
EOF
#  $LON $LAT
#EOF

X=0
Y=0
echo "Inverse: x = $X, y = $Y"
$PROJ -I -m 1:$DX -f %9.6f +proj=stere +R=$RADIUS +lat_0=90.0 +lat_ts=$LAT_TS +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON <<EOF
0.0 0
EOF


