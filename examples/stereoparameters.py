"""polarparameters.py

Script to find parameters for a polar stereographic projection from a ROMS file.

Known parameter may be explicitly prescribed initially. 
If the pm and/or angle variables is missing, the result may be inaccurate
for instance dx = 800.007 and ylon = 70.001 for NorKyst800 (version 2)
Put in reasonable nice rounded values for dx and ylon and rerun the script.

If the ellipsoid is unknown, try both.

"""

# --------------------------------
# Bjørn Ådlandsvik <bjorn@hi.no>
# Institute of Marine Research
# 2022-06-18
# --------------------------------


from math import pi
from pathlib import Path

import pyproj
from netCDF4 import Dataset

# --- Settings ---

datadir = Path("/home/bjorn/data/osea/scratch/ROMS/NorKyst-800m_Lus2/his/2020")
# datadir = Path("gpfs/gpfs0/osea/scratch/ROMS/NorKyst-800m_Lus2/his/2020")
romsfile = "norkyst_800m_his.nc4_2020122801-2020122900"


# Parameter settings
# ellipsis = "WGS84"
ellipsis = "sphere"

dx = None
ylon = None
xp = None
yp = None

# -----------

DEG = 180.0 / pi

proj_params = dict(proj="stere", lat_0=90, lat_ts=60)
if ellipsis == "WGS84":
    proj_params["ellps"] = "WGS84"
else:
    proj_params["R"] = 6371000  # Spherical with standard earth radius
proj = pyproj.Proj(proj_params)

ncid = Dataset(datadir / romsfile)
lon0, lat0 = ncid.variables["lon_rho"][0, 0], ncid.variables["lat_rho"][0, 0]

# Find dx
if dx is None:
    try:  # Easier and more accurate solution if pm variable is present
        pm0 = ncid.variables["pm"][0, 0]
        dx0 = 1 / pm0
    except KeyError:  # No pm variable
        lon1 = ncid.variables["lon_rho"][1, 0]
        lat1 = ncid.variables["lat_rho"][1, 0]
        geod = pyproj.Geod(ellps=proj_params["ellps"])
        _, _, dx0 = geod.inv(lon0, lat0, lon1, lat1)

    scale = proj.get_factors(lon0, lat0).parallel_scale
    dx = round(scale * dx0, 3)
print(f"dx   = {dx:9.3f}")

# Find ylon
if ylon is None:
    try:  # Easier and more accurate solution if angle variable is present
        angle0 = ncid.variables["angle"][0, 0]
        ylon = round(float(lon0 + angle0 * DEG), 3)
    except KeyError:  # No angle variable
        lon1 = ncid.variables["lon_rho"][1, 0]
        lat1 = ncid.variables["lat_rho"][1, 0]
        geod = pyproj.Geod(ellps=proj_params["ellps"])
        a0, _, _ = geod.inv(lon0, lat0, lon1, lat1)
        ylon = round(lon0 - a0, 3)

print(f"ylon = {ylon:9.3f}")
# Update the projection
proj_params["lon_0"] = ylon
proj = pyproj.Proj(proj_params)

# Find xp, yp
if xp is None or yp is None:
    x0, y0 = proj(lon0, lat0)
    xp = -round(x0 / dx, 3)
    yp = -round(y0 / dx, 3)

print(f"xp   = {xp:9.3f}")
print(f"yp   = {yp:9.3f}")
# Update the projection
proj_params["x_0"] = xp * dx
proj_params["y_0"] = yp * dx
proj = pyproj.Proj(proj_params)
print("proj4string = ", proj.to_proj4())

# Check
# Opposite corner
jmax, imax = ncid.variables["lon_rho"].shape
lon1 = ncid.variables["lon_rho"][jmax - 1, imax - 1]
lat1 = ncid.variables["lat_rho"][jmax - 1, imax - 1]
x1, y1 = [v / dx for v in proj(lon1, lat1)]
print()
print(" --- test ---")
# The values should be equal
print(f"X test: {imax - 1}, {x1:.3f}")
print(f"Y test: {jmax - 1}, {y1:.3f}")
