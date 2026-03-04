import numpy as np
from spacecraft import spacecraft, reaction_wheel_system_basic

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

def disturbances(I, #principal inertia tensor of spacecraft 
                 DCM_I2B, #DCM that transforms the inertial reference frame to the body fixed reference frame
                 sc_dipole_moment, #comes from spacecraft class
                 mu, #comes from circular orbit class
                 ECI_position, #comes from circular orbit class
                 B_0, #magnetic field at the equator
                 R_eq #Earth's equatorial radius
                 ):
    Rc = DCM_I2B @ ECI_position               #vector from center of Earth to spacecraft center, given in body coordinates
    gravity_grad_torque =((3*mu)/(np.linalg.norm(Rc))**5)*np.cross(Rc, (I @ Rc)) 
    #next solve for magnetic field torque
    latitude = np.arcsin(ECI_position([2])/np.linalg.norm(ECI_position))
    B_NED = B_0 * ((R_eq/np.linalg.norm(ECI_position))**3)*np.array([[np.cos(latitude)], [0], [2*np.sin(latitude)]])
    #calculate the DCM from NED coordinates to Earth Centered Inertial (ECI)
    n_z_i = -(ECI_position/np.linalg.norm(ECI_position))
    n_y_i = np.cross(n_z_i, np.array([0], [0], [1])) / np.linalg.norm(np.cross(n_z_i, np.array([0], [0], [1])))
    n_x_i = np.cross(n_y_i, n_z_i)
    dcm_i2ned = np.array([[n_z_i.T], [n_y_i.T], [n_x_i.T]])
    dcm_ned2i = dcm_i2ned.T
    dcm_ned2b = DCM_I2B @ dcm_ned2i
    B_B = dcm_ned2b @ B_NED
    mag_torque = np.cross(sc_dipole_moment, B_B) 
    #total disturbance torque is the sum of the gravity gradient torque and magnetic field torque 
    total_torque = gravity_grad_torque  + mag_torque
    return total_torque
    

def eulers_eq_of_rotation_RW( #equation of motion using S&J 4.140
                             I_RW, #get from the reaction_wheel class
                             angular_velocity, #body frame w.r.t the inertial, in body coordinates     
                             hs, #hs vector
                             u, #reaction wheel torque vector
                             l  #external torque vector
                             ):
    angular_acceleration = np.linalg.inv(I_RW) @ ((np.cross(-angular_velocity, (I_RW @ angular_velocity))) - np.cross(angular_velocity, hs) - u + l) 
    return angular_acceleration            