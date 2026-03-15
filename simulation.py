import numpy as np
import dynamics 
import math_functions
from astrodynamics import circular_orbit
import spacecraft
import ADCS

#create the circular orbit
earth_leo_orbit = circular_orbit(altitude = 1500e3, inclination_angle= (45*(np.pi / 180)),
                                  planet_radius = 6378e3, planet_mu = 3.986e14 )
#create the spacecraft 
my_spacecraft = spacecraft.spacecraft(controller=PD, q=0, w=0,r=0,v=0)
#define some important simulation parameters:
t_step = 1 #seconds
t_end = 600 #simulation duration is 600 seconds
time_space = np.linspace(0,t_step,t_end)


#get DCM from inertial to nominal body attitude at a time t:
def get_nominal_body(t):
    # start by getting the nominal body frame at t = 0
    initial_position = earth_leo_orbit.ECI_3d_position(0)
    initial_velocity = earth_leo_orbit.ECI_3d_velocity(0)
    initial_norm_ang_velocity = np.cross(initial_position, initial_velocity)/np.linalg.norm(np.cross(initial_position, initial_velocity))

    nominal_initial_z = -initial_position/np.linalg.norm(initial_position)
    nominal_initial_y = -initial_norm_ang_velocity
    nominal_initial_x = np.cross(nominal_initial_y, nominal_initial_z)
    #create a DCM that describes the transformation from the inertial to the initial nominal body frame
    dcm_initial_nominal = np.array([[nominal_initial_x.T], [nominal_initial_y.T], [nominal_initial_z.T]])
    #to remain pointing at the planet, the body frame of the spacecraft
    #...is rotating with the same angular velocity as the orbit. we can
    #...construct a dcm that takes into account this rotation
    #... this rotation occurs about the same direction as the normalized angular velocity of 
    ###...the orbit. The angular velocity of the orbit is constant
    rotation_angle = (earth_leo_orbit.w) * t
    dcm_orbit_rotation = np.array([[np.cos(rotation_angle),0,np.sin(rotation_angle)],
                            [0,1,0],
                            [-np.sin(rotation_angle),0,np.cos(rotation_angle)]])
    #to get the dcm that transforms the inertial frame to the nominal body frame at a certain
    #...time 't', we must combine dcm_initial_nominal with dcm_orbit_rotation
    dcm_inertial2body_nominal = np.matmul(dcm_orbit_rotation, dcm_initial_nominal)
    return dcm_inertial2body_nominal




def quaternion_evolution(current_quaternion, t_step, 
                         angular_velocity #body wrt inertial, given in body coordinates
                        ): #a function that uses the kinematic differential equation to find the quaternion at the next time step
    quat_dot = dynamics.kinematic_diffeq_quaternion(current_quaternion, angular_velocity)
    future_quaternion = quat_dot * t_step + current_quaternion
    return future_quaternion

def angular_velocity_evolution(current_angular_velocity, t_step, 
                                rw_sys, #insert reaction wheel class
                                I_s, #principal spacecraft inertia
                                L #external torques get from disturbances function within dynamics
                                ): 
    J_RW = rw_sys.principal_moments()
    I_RW = rw_sys.calculate_I_RW(I_s, J_RW)
    hs = rw_sys.calculate_hs_vector(J_RW, rw_sys.spin_speed, current_angular_velocity)
    angular_accel = dynamics.eulers_eq_of_rotation_RW(I_RW, current_angular_velocity, hs, rw_sys.u, L )
    next_angular_velocity = current_angular_velocity + angular_accel * t_step
    #we need to update the spin speeds of our reaction wheels:
    wheel_accel = rw_sys.calculate_wheel_acceleration(angular_accel, J_RW)
    rw_sys.spin_speed = rw_sys.spin_speed + wheel_accel * t_step
    return next_angular_velocity

###INITIAL CONDITIONS###
#step 1: put the spacecraft in a random orientation at time = 0s 
axis = np.array([np.random.normal(0, 10), np.random.normal(0,10), np.random.normal(0,10)])
axis_n = axis / np.linalg.norm(axis)
random_angle = np.random.normal(0, 100) * (np.pi/180) #random angle where 100 degrees is the standard deviation
random_rotation_quaternion = np.array([np.cos(random_angle/2), axis_n[0]*np.sin(random_angle/2), axis_n[1]*np.sin(random_angle/2), axis_n[2]*np.sin(random_angle/2)])
initial_position_quat = math_functions.quaternion_multiply(random_rotation_quaternion, get_nominal_body(0))

#step 2: put the spacecraft in a random angular velocity at time = 0s
initial_angular_velocity = np.array([[np.random.normal(0, (15 * np.pi/180))],[np.random.normal(0, (15 * np.pi/180))],[np.random.normal(0, (15 * np.pi/180))]])


for t in time_space: #no state estimation
    nominal_body_DCM = get_nominal_body(t)
    nominal_ang_vel = np.array([[0],[-earth_leo_orbit.w], [0]]) #requirement to keep the body z axis pointing at Earth
    if t == 0: 
        true_current_position_DCM = math_functions.q_to_dcm(initial_position_quat)
        true_current_angular_velocity = initial_angular_velocity
        imu_bias = ADCS.create_IMU_bias(5.56E-4) #input bias repeatability as deg/s
    #step 3: estimate positions and velocity
    sun_dir_inertial = earth_leo_orbit.sun_i
    sun_sensor_reading = ADCS.simulate_sunsensor(sun_dir_inertial, true_current_position_DCM, .2)
    mag_field_inertial = earth_leo_orbit.calculate_magfield(earth_leo_orbit.ECI_3d_position(t))
    magnetometer_reading = ADCS.simulate_magnetometer(mag_field_inertial, true_current_position_DCM)
    measured_current_position_DCM = ADCS.TRIAD_AD(mag_field_inertial, magnetometer_reading, sun_dir_inertial, sun_sensor_reading)
    measured_angular_velocity = ADCS.simulate_IMU(imu_bias, 6.33E-3,t_step, true_current_angular_velocity)
    #step 4: create control signal (a.k.a reaction wheel torques 'u') 
    u = ADCS.PD_Control_RW(Kp=.2, Kd = .01, DCM_estimate=measured_current_position_DCM, DCM_nominal=nominal_body_DCM, 
                           ang_vel_estimate=measured_angular_velocity, ang_vel_nominal=nominal_ang_vel, principal_inertia=my_spacecraft.I )
    
        


