
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm

import codecs

#with codecs.open('u.mat', encoding='utf-8-sig') as f:
#    u = [[float(x) for x in line.split()] for line in f]
u = np.transpose(np.loadtxt('u.mat'))
#print(u)
x = np.loadtxt('original_mesh_x.mat')
y = np.loadtxt('original_mesh_y.mat')
tri = np.loadtxt('original_mesh_tri.mat').astype(int)
tri = tri - 1
#print(tri)

#triang = mtri.Triangulation(x.ravel(), y.ravel(), tri)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
surf = ax.plot_trisurf(x,y,tri,u,cmap=cm.coolwarm)
fig.colorbar(surf)
plt.show()

