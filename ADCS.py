import numpy as np
from math_functions import create_skew

def TRIAD_AD( #Attitude Determination via the TRIAD Method
              mag_i, #Earth's magnetic field in inertial coordinates, calculated by astrodynamics
              mag_b,  #Earth's magnetic field in body coordinates, measured by the on-board magnetometer
              sun_i,        #direction of the Sun in inertial coordinates, calculated by astrodynamics
              sun_b         #direction of the Sun in body coordinates, calculated by the on-board sun sensors
              ): 
    #create triad in the body frame
    t1_b = sun_b/np.linalg.norm(sun_b)
    t2_b = np.cross(sun_b, mag_b)/np.linalg.norm(np.cross(sun_b, mag_b))
    t3_b = np.cross(t1_b, t2_b)
    #create triad in the inertial frame
    t1_i = sun_i/np.linalg.norm(sun_i)
    t2_i = np.cross(sun_i, mag_i)/np.linalg.norm(np.cross(sun_i, mag_i))
    t3_i = np.cross(t1_i, t2_i)
    DCM_BT = np.hstack((t1_b, t2_b, t3_b))
    DCM_IT = np.hstack((t1_i, t2_i, t3_i))
    DCM_TI = DCM_IT.T
    DCM_BI = DCM_BT @ DCM_TI
    return DCM_BI #a 3x3 numpy matrix describing the rotation from the inertial frame to the spacecraft body fixed frame

def simulate_magnetometer(B_i, #magnetic field in inertial coordinates, get from astrodynamics module
                          DCM_I2B, 
                          noise_std = 40E-9 #standard deviation of sensor, in Teslas
                          ):
    #noise is of Guassian Distribution with mean of 0 and a standard deviation equal to resolution of the magnetometer 
    random_noise_x = np.random.normal(0, noise_std)
    random_noise_y = np.random.normal(0, noise_std)
    random_noise_z = np.random.normal(0, noise_std)
    reading = DCM_I2B @ B_i + np.array([[random_noise_x], [random_noise_y], [random_noise_z]])
    return reading

def sensor_noise_DCM(accuracy): #converting rotational accuracy into a noise matrix. input rotational accuracy in degrees. 
                                #accuracy is assumed to be the same for all three body fixed axes
    accuracy = np.deg2rad(accuracy)
    x_error_angle = np.random.normal(0, accuracy)
    y_error_angle = np.random.normal(0, accuracy)
    z_error_angle = np.random.normal(0, accuracy)
    rotation_angle = np.sqrt(np.square(x_error_angle) + np.square(y_error_angle)+np.square(z_error_angle))
    unit_rotation_axis = np.array([[x_error_angle], [y_error_angle], [z_error_angle]]) / rotation_angle
    e = create_skew(unit_rotation_axis)
    #use Rodrigues' rotation formula to create the rotation matrix after obtaining the rotation angle and the unit rotation axis
    Error_Rotation_Matrix = np.eye(3) + (np.sin(rotation_angle) * e) + ((1-np.cos(rotation_angle)) * e @ e)
    return Error_Rotation_Matrix


def simulate_sunsensor(sun_i, #sun direction in inertial coordinates, get from astrodynamics module
                       DCM_I2B,
                       sun_sensor_accuracy #accuracy in degrees
                       ):
    sun_b = sensor_noise_DCM(sun_sensor_accuracy) @ DCM_I2B @ sun_i
    return sun_b

def create_IMU_bias(BR):
    BR = np.deg2rad(BR)
    return np.array([[np.random.normal(0,BR)],[np.random.normal(0,BR)],[np.random.normal(0,BR)]])
    

def simulate_IMU(bias, #deg/s
                 ARW, #deg/s^.5
                 t_step, true_angular_velocity):
    ARW = np.deg2rad(ARW) 
    gyro_noise_std = ARW/np.sqrt(t_step)
    gyro_noise = np.array([[np.random.normal(0,gyro_noise_std)],[np.random.normal(0,gyro_noise_std)],[np.random.normal(0,gyro_noise_std)]])
    measured_angular_velocity = true_angular_velocity + bias + gyro_noise
    return measured_angular_velocity

def simulate_imperfect_RW(u, torque_error):
    rw1_error = np.random.normal(0, torque_error)
    rw2_error = np.random.normal(0, torque_error)
    rw3_error = np.random.normal(0, torque_error)
    error = np.array([[rw1_error], [rw2_error], [rw3_error]])
    imperfect_u = u + error
    return imperfect_u


def PD_Control_RW(Kp, Kd, DCM_estimate, DCM_nominal, ang_vel_estimate, ang_vel_nominal, RW, I_s): #let this function output the required reaction wheel torques 'u'
    #algorithm obtained from section 6.2 of Gavin Martin ADCS Project
    Kp = Kp * np.eye(3)
    Kd = Kd * np.eye(3)
    error_DCM = DCM_estimate @ DCM_nominal.T
    error_vector = 0.5 * np.array([[error_DCM[2, 1]-error_DCM[1,2]],[error_DCM[0,2]-error_DCM[2,0]], [error_DCM[1,0]-error_DCM[0,1]]])
    error_ang_vel = ang_vel_estimate - ang_vel_nominal
    error_vector_dot = -np.cross(ang_vel_estimate, error_vector) + error_ang_vel
    alpha = -Kd*error_vector_dot - Kp*error_vector #alpha represents the angular acceleration that the controller would like to enact
    J_RW = RW.principal_moments()
    I_RW = RW.calculate_I_RW(I_s, J_RW)
    hs_vector = RW.calculate_hs_vector(J_RW, ang_vel_estimate)
    u = RW.calculate_u(alpha, I_RW, ang_vel_estimate, hs_vector)
    return u
    



