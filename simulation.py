import numpy as np
import dynamics 
import frame_transforms
from astrodynamics import circular_orbit


#create the circular orbit
earth_leo_orbit = circular_orbit(altitude = 1500e3, inclination_angle= (45*(np.pi / 180)),
                                  planet_radius = 6378e3, planet_mu = 3.986e14 )


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
    dcm_initial_nominal = np.array([nominal_initial_x, nominal_initial_y, nominal_initial_z])
    #to remain pointing at the planet, the body frame of the spacecraft
    #...is rotating with the same angular velocity as the orbit. we can
    #...construct a dcm that takes into account this rotation
    #... this rotation occurs about the same direction as the normalized angular velocity of 
    ###...the orbit. The angular velocity of the orbit is constant
    rotation_angle = (earth_leo_orbit.w) * t
    dcm_orbit_rotation = np.array([np.cos(rotation_angle),0,np.sin(rotation_angle)],
                            [0,1,0],
                            [-np.sin(rotation_angle),0,np.cos(rotation_angle)])
    #to get the dcm that transforms the inertial frame to the nominal body frame at a certain
    #...time 't', we must combine dcm_initial_nominal with dcm_orbit_rotation
    dcm_inertial2body_nominal = np.matmul(dcm_orbit_rotation, dcm_initial_nominal)
    return dcm_inertial2body_nominal

     

 