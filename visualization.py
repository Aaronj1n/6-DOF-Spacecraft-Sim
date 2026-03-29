from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import math_functions
import functools

def animate(quat_data):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection= '3d')

    def initialize_box_sc():
        vertices = np.zeros((8,3))
        vertices[0,:] = np.array([0.5,0.5,0.5])
        vertices[1,:] = np.array([0.5,-0.5,0.5])
        vertices[2,:] = np.array([-0.5,-0.5,0.5])
        vertices[3,:] = np.array([-0.5,0.5,0.5])
        vertices[4,:] = np.array([0.5,0.5,-0.5])
        vertices[5,:] = np.array([0.5,-0.5,-0.5])
        vertices[6,:] = np.array([-0.5,-0.5,-0.5])
        vertices[7,:] = np.array([-0.5,0.5,-0.5])
        return vertices
    
    def draw(k): 
        R = math_functions.q_to_dcm(quat_data[k, 1:])
        ax.cla()
        vertices_b = initialize_box_sc()
        vertices_i = np.zeros((8,3))
        for k in range(8):
            vertices_i[k,:] = (R.T @ vertices_b[k,:].reshape((-1,1))).flatten() 
            ax.scatter(vertices_i[k,0], vertices_i[k,1], vertices_i[k,2], c='red')
        
        
    

    # rotate box vertices by R and redraw
    ani = animation.FuncAnimation(fig, draw, frames=len(quat_data))
    ani.save('sim.gif', writer='pillow')