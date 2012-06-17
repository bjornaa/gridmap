# -*- coding: utf-8 -*-

"""Plot a coast line made by makecoast"""

import numpy as np
import matplotlib.pyplot as plt

coastfile = "demo10km_coast.dat"

X, Y = np.loadtxt(coastfile, unpack=True)

I = ~np.isnan(X)
print "X: min, max = ", X[I].min(), X[I].max()
print "Y: min, max = ", Y[I].min(), Y[I].max()


# Fill land
plt.fill(X, Y, color='green')
# Add a black coast line
plt.plot(X, Y, color='black')

plt.show()

