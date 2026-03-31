import numpy as np
import matplotlib.pyplot as plt
import dynamics 
import math_functions
from astrodynamics import circular_orbit
from spacecraft import spacecraft, reaction_wheel_system_basic
import ADCS
import math
from functools import partial
from visualization import animate



#create the circular orbit
earth_leo_orbit = circular_orbit(altitude = 1500e3, inclination_angle= (45*(np.pi / 180)),
                                  planet_radius = 6378e3, planet_mu = 3.986e14, B_0=3E-5 )
#create the spacecraft 
my_spacecraft = spacecraft(controller='PD', q=0, w=0,r=0,v=0)
rw= reaction_wheel_system_basic(wheel_diameter=.049, wheel_height=.0175, wheel_mass= .197, max_torque = 8E-3)
#define some important simulation parameters:
t_step = .01 #seconds
t_end = 600 #simulation duration is 600 seconds
time_space = np.arange(0,t_end,t_step)
print(len(time_space))
#set up data collection:
N = len(time_space)
quat_data = np.ones((N,4))
rw_u_data = np.ones((N,3))
rw_spin_speed_data = np.ones((N,3))
vel_data = np.ones((N,3))
error_data = np.ones((N,1))

#clear previous test_log
open('test_log.txt', 'w').close()


#get DCM from inertial to nominal body attitude at a time t:
def get_nominal_body(t):
    # start by getting the nominal body frame at t = 0
    initial_position = earth_leo_orbit.ECI_3d_position(0)
    initial_velocity = earth_leo_orbit.ECI_3d_velocity(0)
    initial_norm_ang_velocity = math_functions.cross(initial_position, initial_velocity)/np.linalg.norm(math_functions.cross(initial_position, initial_velocity))

    nominal_initial_z = -initial_position/np.linalg.norm(initial_position)
    nominal_initial_y = -initial_norm_ang_velocity
    nominal_initial_x = math_functions.cross(nominal_initial_y, nominal_initial_z)
    #create a DCM that describes the transformation from the inertial to the initial nominal body frame
    dcm_initial_nominal = np.vstack([nominal_initial_x.T, nominal_initial_y.T, nominal_initial_z.T])
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
initial_position_quat = math_functions.quaternion_multiply(random_rotation_quaternion, math_functions.dcm_to_q(get_nominal_body(0)))

#step 2: put the spacecraft in a random angular velocity at time = 0s
initial_angular_velocity = np.array([[np.random.normal(0, (15 * np.pi/180))],[np.random.normal(0, (15 * np.pi/180))],[np.random.normal(0, (15 * np.pi/180))]])

print('orbital w:', -earth_leo_orbit.w)
for r in range(N): #no state estimation
    t = time_space[r]
    nominal_body_DCM = get_nominal_body(t)
    nominal_ang_vel = np.array([[0],[-earth_leo_orbit.w], [0]]) #requirement to keep the body z axis pointing at Earth
    if t == 0: 
        current_quat = initial_position_quat
        true_current_position_DCM = math_functions.q_to_dcm(current_quat)
        true_current_angular_velocity = initial_angular_velocity
        imu_bias = ADCS.create_IMU_bias(5.56E-4) #input bias repeatability as deg/s
    #if r > 0 and np.linalg.norm(true_current_angular_velocity) > 100:  # threshold for "too big"
        # print(f"Blowup at step {r}, t={r*t_step:.3f}s")
        # print(f"  sun_b norm:              {np.linalg.norm(sun_sensor_reading)}")
        # print(f"  mag_b norm:              {np.linalg.norm(mag_b)}")
        # print(f"  cross(sun,mag) norm:     {np.linalg.norm(cross(sun_b, mag_b))}")
        # print(f"  Rc norm:                 {np.linalg.norm(Rc)}")

    
    #step 3: estimate positions and velocity
    sun_dir_inertial = earth_leo_orbit.sun_i
    sun_sensor_reading = ADCS.simulate_sunsensor(sun_dir_inertial, true_current_position_DCM, .2)
    if math.isnan(sun_sensor_reading[1,0]):
        print('breaking')
        break
    mag_field_inertial = earth_leo_orbit.calculate_magfield(earth_leo_orbit.ECI_3d_position(t))
    magnetometer_reading = ADCS.simulate_magnetometer(mag_field_inertial, true_current_position_DCM)
    measured_current_position_DCM = ADCS.TRIAD_AD(mag_field_inertial, magnetometer_reading, sun_dir_inertial, sun_sensor_reading)
    measured_angular_velocity = ADCS.simulate_IMU(imu_bias, 6.33E-3,t_step, true_current_angular_velocity)
    #step 4: create control signal (a.k.a reaction wheel torques 'u') 
    u, error_vector_magnitude = ADCS.PD_Control_RW(Kp=10, Kd = 10, DCM_estimate=measured_current_position_DCM, DCM_nominal=nominal_body_DCM, 
                           ang_vel_estimate=measured_angular_velocity, ang_vel_nominal=nominal_ang_vel, RW=rw, I_s=my_spacecraft.I )
    #imperfect_u = ADCS.simulate_imperfect_RW(u, torque_error= .04E-3)
    imperfect_u = ADCS.simulate_imperfect_RW(u, torque_error= 0)
    #step 5: create disturbances
    disturbances = dynamics.disturbances(I=my_spacecraft.I, DCM_I2B=true_current_position_DCM, sc_dipole_moment= my_spacecraft.b, 
                              mu=earth_leo_orbit.planet_mu, ECI_position=earth_leo_orbit.ECI_3d_position(t), 
                              B_i=earth_leo_orbit.calculate_magfield(earth_leo_orbit.ECI_3d_position(t)))
    #step 6: obtain the true angular acceleration due to environment + rw actuators
    
    rw.hs = rw.calculate_hs_vector(J_RW=rw.principal_moments(),angular_velocity_body=true_current_angular_velocity)
    I_RW = rw.calculate_I_RW(I_s=my_spacecraft.I, J_RW= rw.principal_moments())
    angular_acceleration = dynamics.eulers_eq_of_rotation_RW(I_RW=I_RW, angular_velocity= true_current_angular_velocity, hs = rw.hs, 
                                                             u=imperfect_u, l=disturbances)
    
    #step 7: update the reaction wheel speed
    rw.u = imperfect_u
    wheel_accel = rw.calculate_wheel_acceleration(angular_acceleration, rw.principal_moments())
    rw.spin_speed = rw.spin_speed + wheel_accel * t_step


    #step 7: propagate  the angular velocity using the angular acceleration
    true_future_angular_velocity = true_current_angular_velocity + angular_acceleration * t_step

    #step 8: propagate the attitude using the angular velocity

    quat_dot = dynamics.kinematic_diffeq_quaternion(current_quat, true_current_angular_velocity)
    future_quat = (current_quat.reshape((-1,1)) + quat_dot * t_step).flatten()
    future_quat = future_quat/np.linalg.norm(future_quat)
    
    #step 9: update data
    quat_data[r,:] = current_quat
    rw_u_data[r,:] = imperfect_u.T
    rw_spin_speed_data[r,:] = (rw.spin_speed).T
    vel_data[r,:] = (true_current_angular_velocity).T
    error_data[r,:] = error_vector_magnitude
    
    #test log:
    with open('test_log.txt', 'a') as f:
        f.write(f't={t:.4f} \n')
        f.write(f'current quat={current_quat} \n')
       # f.write(f'true_current_position_DCM={true_current_position_DCM} \n')
        f.write(f'true current angular velocity={true_current_angular_velocity} \n')
        f.write(f'quat_dot={quat_dot} \n')
        f.write(f'wheel u={rw.u} \n')
        f.write(f'true current angular acceleration = {angular_acceleration} \n')
        f.write(f'Error mag: {error_vector_magnitude}\n\n')

    # print('future quat:', future_quat)
    # print('true future angular velocity:', true_future_angular_velocity)
    

    #update variables for the next iteration of for loop
    current_quat = future_quat
    true_current_position_DCM = math_functions.q_to_dcm(current_quat)
    true_current_angular_velocity = true_future_angular_velocity
   


#give all the data arrays a time column for readability
stackable_time = time_space.reshape(-1,1)
# quat_data = np.hstack([stackable_time, quat_data])
# rw_u_data = np.hstack([stackable_time, rw_u_data])
# rw_spin_speed_data = np.hstack([stackable_time, rw_spin_speed_data])
# vel_data = np.hstack([stackable_time, vel_data])


# animate(quat_data)

#plot error:
plt.plot(stackable_time, error_data, color='red')
plt.xlim(0,t_end)
plt.grid(True, alpha=.3)
plt.xlabel('Time (s)')
plt.ylabel('Error (radians)')
plt.title('Error vs Time')


#plot quaternion

fig_quat, axes_quat = plt.subplots(4, 1, figsize=(10, 8), sharex=True)

axes_quat[0].plot(stackable_time, quat_data[:,0], 'y')
axes_quat[0].set_ylabel('q0')

axes_quat[1].plot(stackable_time, quat_data[:,1], 'r')
axes_quat[1].set_ylabel('q1')

axes_quat[2].plot(stackable_time, quat_data[:,2], 'b')
axes_quat[2].set_ylabel('q2')

axes_quat[3].plot(stackable_time, quat_data[:,3], 'g')
axes_quat[3].set_ylabel('q3')
axes_quat[3].set_xlabel('Time (s)')

plt.tight_layout()   # prevents labels overlapping


#plot angular velocity
fig_omega, axes_omega = plt.subplots(3,1)
axes_omega[0].plot(stackable_time,vel_data[:,0])
axes_omega[0].set_ylabel('omega x (rad/s)')

axes_omega[1].plot(stackable_time, vel_data[:,1])
axes_omega[1].set_ylabel('omega y (rad/s)')

axes_omega[2].plot(stackable_time, vel_data[:,2])
axes_omega[2].set_ylabel('omega z (rad/s)')
axes_omega[2].set_xlabel('Time (s)')



#plot reaction wheel torques

fig_rw, axes_rw = plt.subplots(3,1)
axes_rw[0].plot(stackable_time,rw_u_data[:,0])
axes_rw[0].set_ylabel('rw #1 torque (add units)')

axes_rw[1].plot(stackable_time, rw_u_data[:,1])
axes_rw[1].set_ylabel('rw #2 torque')

axes_rw[2].plot(stackable_time, rw_u_data[:,2])
axes_rw[2].set_ylabel('rw #3 torque')
axes_rw[2].set_xlabel('Time (s)')

plt.tight_layout()
plt.show()


