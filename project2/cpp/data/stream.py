import sys
import os
import matplotlib.pyplot as plt
import numpy as np

#filename = 're1000_central/final-snapshot'
filename = sys.argv[1]
data = np.loadtxt(filename, skiprows=1)
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]
v = data[:, 3]


xi = np.unique(x)
yi = np.unique(y)
#nx = len(xi)
#ny = len(yi)
nx = 1000
ny = 1000

xlin = np.linspace(x.min(),x.max(), nx)
ylin = np.linspace(y.min(),y.max(), ny)
X, Y = np.meshgrid(xlin, ylin)

from scipy.interpolate import griddata

U = griddata((x, y), u, (X, Y), method='cubic')
V = griddata((x, y), v, (X, Y), method='cubic')

MAG = np.sqrt((U*U) + (V*V))

plt.figure(figsize=(6,6))

plt.streamplot(X, Y, U, V, color='black', density=1.6)

#plt.title('RE = 1000 | Central')
plt.title(filename)
plt.xlabel('X')
plt.ylabel('Y')

plt.show()

