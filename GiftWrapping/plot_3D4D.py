import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def onKeyEvent(event):
    angle_step = .1    
    
    if getattr(event, 'key') == 'q' or getattr(event, 'key') == 'z':
		if getattr(event, 'key') == 'q':
			phi =  angle_step
		else:
			phi = -angle_step
        
		rotation_matrix = np.array([[np.cos(phi), 0., 0., -np.sin(phi)], 
									[		  0., 1., 0., 			0.],
									[		  0., 0., 1., 			0.],
									[np.sin(phi), 0., 0., np.cos(phi)]])

    elif getattr(event, 'key') == 'w' or getattr(event, 'key') == 'x':
		if getattr(event, 'key') == 'w':
			phi =  angle_step
		else:
			phi = -angle_step
		
		rotation_matrix = np.array([[1., 		  0., 0., 			0.],
									[0., np.cos(phi), 0., -np.sin(phi)], 
									[0., 		  0., 1., 			0.],
									[0., np.sin(phi), 0., np.cos(phi)]])
		
    elif getattr(event, 'key') == 'e' or getattr(event, 'key') == 'c':
		if getattr(event, 'key') == 'e':
			phi =  angle_step
		else:
			phi = -angle_step
		
		rotation_matrix = np.array([[1., 0., 		  0., 			0.],
									[0., 1., 		  0., 			0.],
									[0., 0., np.cos(phi), -np.sin(phi)], 
									[0., 0., np.sin(phi), np.cos(phi)]])

    else:
		rotation_matrix = np.array([[1., 0., 0., 0.],
									[0., 1., 0., 0.],
									[0., 0., 1., 0.], 
									[0., 0., 0., 1.]])

    for i in xrange(len(points)):
	    points[i] = rotation_matrix.dot(points[i])  

    projected_points = np.array(map(lambda p: (p / (p[3]-3.))[0:3], points))

    ax.clear() 
    ax.set_xlim(-1., 1.)
    ax.set_ylim(-1., 1.)
    ax.set_zlim(-1., 1.)

	# ~~~~~
    """
    for i in xrange(len(projected_points)):
        ax.text(projected_points[i,0], projected_points[i,1], projected_points[i,2], i, None)
	"""
    # ~~~~~
       
    poly_verts = np.array([projected_points[map(int, faces[counter])] for counter in xrange(len(faces))])
    poly = Poly3DCollection(poly_verts, facecolor=[0., 1., 0.], edgecolor=[0., 0., 0.], alpha=0.1)

    ax.scatter(*zip(*projected_points))
    ax.add_collection3d(poly)
    fig.canvas.draw()


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

"""
# Icosahedron
points = np.array([[-0.692,  0.000,  0.427],
				   [ 0.000,  0.427, -0.692],
				   [ 0.000,  0.427,  0.692],
				   [ 0.692,  0.000, -0.427],
				   [-0.427, -0.692,  0.000],
				   [-0.427,  0.692,  0.000],
				   [ 0.000, -0.427,  0.692],
				   [ 0.427,  0.692,  0.000],
				   [ 0.000, -0.427, -0.692],
				   [ 0.692,  0.000,  0.427],
				   [ 0.427, -0.692,  0.000],
				   [-0.692,  0.000, -0.427]])

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

"""

"""
# Tesseract
points = np.array([[-1., -1., -1., -1.], # 0
				   [-1., -1., -1.,  1.], # 1
				   [-1., -1.,  1., -1.], # 2
				   [-1., -1.,  1.,  1.], # 3
                   [-1.,  1., -1., -1.], # 4
				   [-1.,  1., -1.,  1.], # 5
                   [-1.,  1.,  1., -1.], # 6
				   [-1.,  1.,  1.,  1.], # 7
                   [ 1., -1., -1., -1.], # 8
				   [ 1., -1., -1.,  1.], # 9
                   [ 1., -1.,  1., -1.], # 10
				   [ 1., -1.,  1.,  1.], # 11
                   [ 1.,  1., -1., -1.], # 12
				   [ 1.,  1., -1.,  1.], # 13
				   [ 1.,  1.,  1., -1.], # 14
				   [ 1.,  1.,  1.,  1.]])# 15


faces = np.array([[10,  8, 12, 14], # Inner cube
                  [14, 12,  4,  6],
                  [ 6,  4,  0,  2],
                  [ 2,  0,  8, 10],
                  [ 2, 10, 14,  6],
                  [ 8,  0,  4, 12],
                 
                  [11,  9, 13, 15], # Outer cube
                  [15, 13,  5,  7],
                  [ 7,  5,  1,  3],
                  [ 3,  1,  9, 11],
                  [ 3, 11, 15,  7],
                  [ 9,  1,  5, 13],
                  
                  [ 8,  9, 13, 12], # Bottom partitions
                  [12, 13,  5,  4],
                  [ 4,  5,  1,  0],
                  [ 0,  1,  9,  8],
                  
                  [10, 11, 15, 14], # Upper partitions
                  [14, 15,  7,  6],
                  [ 6,  7,  3,  2],
                  [ 2,  3, 11, 10],
                 
                  [11, 10,  8,  9], # Side partitions
                  [15, 14, 12, 13],
                  [ 7,  6,  4,  5],
                  [ 3,  2,  0,  1]])
"""

is_4D = (len(points[0]) == 4)

if is_4D:
    # $$$$$
    faces = np.array([[face[index] for index in xrange(4) if index != index_skip] for face in faces for index_skip in xrange(4)])
    # $$$$$

    projected_points = np.array(map(lambda p: (p / (p[3]-3.))[0:3], points))
    poly_verts = np.array([projected_points[map(int, faces[counter])] for counter in xrange(len(faces))])
    poly = Poly3DCollection(poly_verts, facecolor=[0., 1., 0.], edgecolor=[0., 0., 0.], alpha=0.1)
else:
	poly_verts = np.array([points[map(int, faces[counter])] for counter in xrange(len(faces))])
	poly = Poly3DCollection(poly_verts, facecolor=[0., 0., 1.], edgecolor=[0., 0., 0.], alpha=0.3)

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim(-1., 1.)
ax.set_ylim(-1., 1.)
ax.set_zlim(-1., 1.)

if is_4D:
	ax.scatter(*zip(*projected_points))
else:
	ax.scatter(*zip(*points))

ax.add_collection3d(poly)

# ~~~~~
"""
if is_4D:
	for i in xrange(len(projected_points)):
		ax.text(projected_points[i,0], projected_points[i,1], projected_points[i,2], i, None)	
else:	
	for i in xrange(len(points)):
		ax.text(points[i,0], points[i,1], points[i,2], i, None)
"""
# ~~~~~

if is_4D:
	key_press_event_id = fig.canvas.mpl_connect('key_press_event', onKeyEvent)
	key_release_event_id = fig.canvas.mpl_connect('key_release_event', onKeyEvent)

plt.show()

if is_4D:
	fig.canvas.mpl_disconnect(key_press_event_id)
	fig.canvas.mpl_disconnect(key_release_event_id)
