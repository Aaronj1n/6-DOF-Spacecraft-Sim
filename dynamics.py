import numpy as np

def kinematic_diffeq_quaternion(quat, w_b):
    #outputs the derivative of the quaternion
    #takes in current quaternion and the angular velocity of the body fixed frame w.r.t 
    #...the inertial frame, written in body fixed coordinates
    beta_0 = quat[0]
    beta_1 = quat[1]
    beta_2 = quat[2]
    beta_3 = quat[3]
    B_of_quaternion = np.array([-beta_1, -beta_2, -beta_3],
                               [beta_0, -beta_3, beta_2],
                               [beta_3, beta_0, -beta_1],
                               [-beta_2, beta_1, beta_0])
    quat_dot = (1/2)*np.matmul(B_of_quaternion, w_b)
    return quat_dot

def eulers_eq_of_rotation(I, w_b, M):
    I_inv = np.linalg.inv(I)
    w_b_dot = np.cross(np.matmul(-I_inv, w_b), np.matmul(I,w_b)) + np.matmul(I_inv, M)
    return w_b_dot

    