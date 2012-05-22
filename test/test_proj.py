
# -*- coding: utf-8 -*-

import sys
import subprocess
from gridmap import *

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_c)   [m]
ylon  = 58.0          # angle of y-axis        [deg]


lon, lat = 2, 66   # Station M

gmap = PolarGridMap(xp, yp, dx, ylon)

#print projstring1
#os.system('proj ' + projstring1)
#subprocess.call('proj', projstring, shell=True)

command1 = 'proj -f %17.13f ' + gmap.projstring1()
command2 = 'proj -f %17.13f ' + gmap.projstring2()

print "sphere"

print command1
p1 = subprocess.Popen(command1, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

p1.stdin.write('%g %g\n' % (lon,lat))
out1, err1 = p1.communicate()

print command2
p2 = subprocess.Popen(command2, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

p2.stdin.write('%g %g\n' % (lon,lat))
out2, err2 = p2.communicate()

print "projstring1  : ", out1.strip()
print "projstring2  : ", [float(w)/gmap.dx for w in out2.split()]
print "gmap.ll2grid : ", gmap.ll2grid(lon, lat)

import sys
sys.exit(0)

print "WGS84"

gmap = PolarGridMap(xp, yp, dx, ylon, ellipsoid=WGS84())

#print projstring1
#os.system('proj ' + projstring1)
#subprocess.call('proj', projstring, shell=True)

command1 = 'proj -f %17.13f ' + gmap.projstring1()
command2 = 'proj -f %13.9f ' + gmap.projstring2()


print command1
p1 = subprocess.Popen(command1, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

p1.stdin.write('%g %g\n' % (lon,lat))
out1, err1 = p1.communicate()

print command2
p2 = subprocess.Popen(command2, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

p2.stdin.write('%g %g\n' % (lon,lat))
out2, err2 = p2.communicate()

print "projstring1  : ", out1.strip()
print "projstring2  : ", [float(w)/gmap.dx for w in out2.split()]
print "gmap.ll2grid : ", gmap.ll2grid(lon, lat)





