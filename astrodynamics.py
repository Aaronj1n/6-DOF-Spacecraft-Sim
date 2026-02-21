import numpy as np

class orbit:
    def __init__(self, altitude, inclination_angle, planet, type = "circular" ):
        self.altitude = altitude
        self.inclination_angle = inclination_angle #inclination angle is about ECI Y-Axis
        self.planet = planet
        self.type = type
    if self.type == "circular":

        self.orbit_radius = (self.altitude + self.planet.radius) 
        self.w = np.sqrt(self.planet.mu)/(r**(3/2))
        def perifocal_2d_motion(t):
            #returns p and q coordinates of the spacecraft in the perifocal plane
            angle = (self.w)*t
            p = self.orbit_radius * np.cos(angle)
            q = self.orbit_radius * np.sin(angle)
            #third direction (p cross q) is zero due to definition of the perifocal plane
            return np.array([p,q,0])
        def perifocal_to_ECI(inclination_angle):
            #returns the DCM that takes the spacecraft position from perifocal frame to Earth
            #...Centered Inertial frame (ECI).
            dcm = np.array([np.cos(inclination_angle), 0, np.cos(inclination_angle)],
                           [0, 1, 0],
                           [-np.cos(inclination_angle), 0, np.cos(inclination_angle)])
            return dcm
        def ECI_3d_motion(t):
            perifocal_coords = perifocal_2d_motion(t)
            dcm = perifocal_to_ECI(self.inclination_angle)
            return np.matmul(dcm, perifocal_coords)
        



class planet:
    def __init__(self, radius, mu):
        self.radius = radius
        self.mu = mu

