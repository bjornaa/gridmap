#!/bin/bash

# Read polar stereographic parameters and grid size
# from a ROMS grid file (with polar stereographic info)
#
# --------------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-06-30
# -------------------------------------

# TODO: finn godt navn, og ta GRIDFILE fra kommando-linjen
# TODO: Ta alltid Lm, Mm, stopp med feil melding
# OBS: Får ikke .25 på xp og yp

GRIDFILE=demo10km_grid.nc

# Temporary CDL-file
CDL_FILE=a${RANDOM}.cdl

ncdump -h $GRIDFILE > $CDL_FILE

line=( `grep "xi_rho = " $CDL_FILE` )
Lp=${line[2]}
let Lm=Lp-2

line=( `grep "eta_rho = " $CDL_FILE` )
Mp=${line[2]}
let Mm=Mp-2

line=( `grep grid_mapping:dx $CDL_FILE` )
dx=${line[2]}

line=( `grep grid_mapping:false_easting $CDL_FILE` )
xp=${line[2]}
xp=$(echo  "scale=3; $xp/$dx" | bc)

line=( `grep grid_mapping:false_northing $CDL_FILE` )
yp=${line[2]}
yp=$(echo "scale=3; $yp/$dx" | bc)

line=( `grep grid_mapping:straight_vertical_longitude_from_pole $CDL_FILE` )
ylon=${line[2]}


# Clean up
rm $CDL_FILE

echo $GRIDFILE
echo Lm, Mm = $Lm, $Mm
echo xp, yp, dx, ylon = $xp, $yp, $dx, $ylon
