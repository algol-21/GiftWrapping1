import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def read_file(file_name):
    res = []
    f = open(file_name, 'r')
    for line in f:
        lst = map(float, line.split())
        res.append(lst)
    f.close()
    return np.array(res)

points = read_file('data/points.txt')
faces  = read_file('data/faces.txt')

# -----
"""
points = np.array([[-0.692, 0.000, 0.427],
          [0.000, 0.427, -0.692],
          [0.000, 0.427, 0.692],
          [0.692, 0.000, -0.427],
          [-0.427, -0.692, 0.000],
          [-0.427, 0.692, 0.000],
          [0.000, -0.427, 0.692],
          [0.427, 0.692, 0.000],
          [0.000, -0.427, -0.692],
          [0.692, 0.000, 0.427],
          [0.427, -0.692, 0.000],
          [-0.692, 0.000, -0.427]])



faces = np.array([[9, 2, 6],
                 [1, 5, 11],
                 [11, 1, 8],
                 [0, 11, 4],
                 [3, 7, 1],
                 [3, 1, 8],
                 [9, 3, 7],
                 [0, 2, 6],
                 [4, 6, 10],
                 [1, 7, 5],
                 [7, 2, 5],
                 [8, 10, 3],
                 [4, 11, 8],
                 [9, 2, 7],
                 [10, 6, 9],
                 [0, 11, 5],
                 [0, 2, 5],
                 [8, 10, 4],
                 [3, 9, 10],
                 [6, 4, 0]])

print points
print faces
"""
# -----

poly_verts = np.array([points[map(int, faces[counter])] for counter in xrange(len(faces))])
poly = Poly3DCollection(poly_verts, facecolor=[0., 0., 1.], edgecolor=[0., 0., 0.], alpha=0.3)

fig = plt.figure()

ax = Axes3D(fig)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.scatter(*zip(*points))
ax.add_collection3d(poly)

plt.show()