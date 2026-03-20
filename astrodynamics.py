import numpy as np
from math_functions import cross
class circular_orbit:
    def __init__(self, altitude, inclination_angle, planet_radius, planet_mu, 
                 B_0, #magnetic field at equator
                 sun_i = np.array([[1], [0], [0]])#sun direction in ECI coordinates, at March 21st
                 ):
        self.altitude = altitude
        self.inclination_angle = inclination_angle #inclination angle is about ECI Y-Axis
        self.planet_radius = planet_radius
        self.planet_mu = planet_mu
        self.type = type
        self.orbit_radius = (self.altitude + self.planet_radius) 
        self.w = np.sqrt(self.planet_mu)/(self.orbit_radius**(3/2))
        self.B_0 = B_0
        self.sun_i = sun_i
        ##def perifocal_2d_motion(t):
            #returns p and q coordinates of the spacecraft in the perifocal plane
           ## angle = (self.w)*t
           ## p = self.orbit_radius * np.cos(angle)
          ##  q = self.orbit_radius * np.sin(angle)
            #third direction (p cross q) is zero due to definition of the perifocal plane
        ##    return np.array([p,q,0])
      ##  def perifocal_to_ECI(inclination_angle):
            #returns the DCM that takes the spacecraft position from perifocal frame to Earth
            #...Centered Inertial frame (ECI).
      ##      dcm = np.array([np.cos(inclination_angle), 0, np.cos(inclination_angle)],
    #                       [0, 1, 0],
     #                      [-np.cos(inclination_angle), 0, np.cos(inclination_angle)])
       ##     return dcm
    def ECI_3d_position(self,t):
            #perifocal_coords = perifocal_2d_motion(t)
           # dcm = perifocal_to_ECI(self.inclination_angle)
            #return np.matmul(dcm, perifocal_coords)
            x = np.cos(self.inclination_angle)* self.orbit_radius*np.cos((self.w)*t)
            y = self.orbit_radius*np.sin((self.w)*t)
            z = -np.cos(self.inclination_angle) * self.orbit_radius*np.cos((self.w)*t)
            return np.array([[x],[y],[z]])
    def ECI_3d_velocity(self,t):
            xdot = -np.cos(self.inclination_angle)*self.orbit_radius*self.w*np.sin((self.w)*t)
            ydot = self.orbit_radius * self.w * np.cos(self.w * t)
            zdot = np.cos(self.inclination_angle)*self.orbit_radius*self.w*np.sin((self.w)*t)
            return np.array([[xdot], [ydot], [zdot]])
    def calculate_magfield(self, ECI_position): #returns the magnetic field in inertial frame
      latitude = np.arcsin(ECI_position[2,0]/np.linalg.norm(ECI_position))
      B_NED = self.B_0 * ((self.planet_radius/np.linalg.norm(ECI_position))**3)*np.array([[np.cos(latitude)], [0], [2*np.sin(latitude)]])
      #calculate the DCM from North-East-Down (NED) coordinates to Earth Centered Inertial (ECI)
      n_z_i = -(ECI_position/np.linalg.norm(ECI_position))
      n_y_i = cross(n_z_i, np.array([[0],[0],[1]])) / np.linalg.norm(cross(n_z_i, np.array([[0], [0], [1]])))
      n_x_i = cross(n_y_i, n_z_i)
      dcm_i2ned = np.vstack([n_z_i.T, n_y_i.T, n_x_i.T])
      dcm_ned2i = dcm_i2ned.T
      B_i = dcm_ned2i @ B_NED
      return B_i
    