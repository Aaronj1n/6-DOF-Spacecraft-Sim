import numpy as np 
from math_functions import dot, cross

class spacecraft:
    def __init__(self, controller, q, w, r, v, 
                 b = np.array([[3.6E-2], [0], [0]]), #units of A x m^2
                 I = np.diag([.13, .13, .08]) #units of kg x m^2
                 ):
        self.I = I
        self.controller = controller
        self.q = q
        self.w = w
        self.r = r
        self.v = v 
        self.b = b


class reaction_wheel_system_basic:    #basic meaning there are three RW's, each aligned with a 
                                      #...body principal axis
    def __init__(self, wheel_diameter, wheel_height, wheel_mass, max_torque, max_ang_momentum, max_spin_speed, spin_speed = np.array([[0],[0],[0]]),
                  u= np.array([[0],[0],[0]]), hs = np.array([[0],[0],[0]])
                  ):
        self.wheel_diameter = wheel_diameter
        self.wheel_radius = wheel_diameter / 2
        self.wheel_height = wheel_height
        self.wheel_mass = wheel_mass
        self.max_torque = max_torque
        self.max_ang_momentum = max_ang_momentum
        self.max_spin_speed = max_spin_speed
        self.spin_speed = spin_speed
        self.u = u
        self.hs= hs

    def principal_moments(self): #calculates moment of inertia of reaction wheel about 
                                 #...reaction wheel coordinate system 
        tensor = np.zeros((3,3))
        tensor[0,0] = ((1/12) *  self.wheel_mass * (self.wheel_height)**2) + ((1/4) *  self.wheel_mass * (self.wheel_radius)**2)
        tensor[1,1] = ((1/12) *  self.wheel_mass * (self.wheel_height)**2) + ((1/4) *  self.wheel_mass * (self.wheel_radius)**2)
        tensor [2,2] = (1/2) * self.wheel_mass * (self.wheel_radius)**2 
        return tensor #transverse, gimbal, spin axes 
    def calculate_I_RW(self, #calculates the spacecraft + RW moment of inertia tensor [I_RW] described in S&J 4.140
                       I_s, #spacecraft moment of inertia about body fixed frame
                       J_RW #principal inertia tensor for RW about reaction wheel frame
                       ): 
        #get the reaction wheel principal moments for transverse and gimbal axes:
        J_t = J_RW[0,0]
        J_g = J_RW[1,1]
        #grab body frame inertia tensor minus spin axis for reaction wheel #1:
        g_t_1 = np.array([[1],[0],[0]])
        g_t_1_tranpose = np.transpose(g_t_1)
        g_g_1 = np.array([[0], [0], [-1]])
        g_g_1_tranpose = np.transpose(g_g_1)
        bf_inertia_rw1 = (J_t * g_t_1 @ g_t_1_tranpose) + (J_g * g_g_1 @ g_g_1_tranpose)
        #grab body frame inertia tensor minus spin axis for reaction wheel #2:
        g_t_2 = np.array([[0],[1],[0]])
        g_t_2_tranpose = np.transpose(g_t_2)
        g_g_2 = np.array([[0],[0],[1]])
        g_g_2_tranpose = np.transpose(g_g_2)
        bf_inertia_rw2 = (J_t * g_t_2 @ g_t_2_tranpose) + (J_g * g_g_2 @ g_g_2_tranpose)
        #grab body frame inertia tensor minus spin axis for reaction wheel #3
        g_t_3 = np.array([[1],[0],[0]])
        g_t_3_tranpose = np.transpose(g_t_3)
        g_g_3 = np.array([[0],[1],[0]])
        g_g_3_tranpose = np.transpose(g_g_3)
        bf_inertia_rw3 = (J_t * g_t_3 @ g_t_3_tranpose) + (J_g * g_g_3 @ g_g_3_tranpose)
        # add all body frame inertia matrices of RW's to the body frame inertia matrix of the entire spacecraft
        I_RW = I_s + bf_inertia_rw1 + bf_inertia_rw2 + bf_inertia_rw3
        return I_RW
    def calculate_wheel_acceleration(self, #function to return a 3x1 vector of the reaction wheel spin accelerations
                                      angular_acceleration_body, #angular acceleration of the spacecraft in body coordinates, 3x1 vector
                                      J_RW #the principal inertia tensor for each reaction wheel in reaction wheel frame coordinates
                                       ):
        J_s = J_RW[2,2]
        g_s_tensor = np.array([[0,1,0], [1,0,0], [0,0,1]])
        wheel_accel_vector = (self.u/(J_s)) - (g_s_tensor @ angular_acceleration_body)
        return wheel_accel_vector
    def calculate_hs_vector(self,  #returns a 3x1 column vector
                           J_RW, 
                           angular_velocity_body #angular velocity of the spacecraft in body coordinates, 3x1 vector
                           ):
        J_s = J_RW[2,2] # the same for each reaction wheel
        s1 = np.array([[0],[1],[0]])  #RW #1 spin direction, written in body coordinates 
        s2 = np.array([[1],[0],[0]])
        s3= np.array([[0],[0],[1]])
        w_s1 = dot(angular_velocity_body, s1)   #component of spacecraft angular velocity in the RW #1 spin direction
        w_s2 = dot(angular_velocity_body, s2)  #component of spacecraft angular velocity in the RW #2 spin direction
        w_s3 = dot(angular_velocity_body, s3) #component of spacecraft angular velocity in the RW #3 spin direction
        h_vector = J_s * (np.array([[w_s1],[w_s2], [w_s3]]) + self.spin_speed)
        return h_vector
    
    def calculate_u(self, desired_ang_accel, I_RW, estimated_angular_vel, hs_vector, J_RW):
        desired_u =  cross(-estimated_angular_vel, (I_RW @ estimated_angular_vel)) - cross(estimated_angular_vel, hs_vector) - (I_RW @ desired_ang_accel)
        empty_list = []
        #check that max angular momentum not exceeded for wheels
        J_s = J_RW[2,2]
        for k in range(3):
            spin_speed = self.spin_speed[k,0]
            # L = J_s * spin_speed
            # if np.abs(L) >= self.max_ang_momentum:
            if spin_speed >= self.max_spin_speed:
                u_c = 0
            else:
                u_c = desired_u[k]
                if np.abs(u_c.item()) <= self.max_torque:
                    u_c = u_c.item()
                else:
                    if u_c.item() > 0:
                        u_c = self.max_torque
                    elif u_c.item() < 0:
                        u_c = -self.max_torque 
                # if u_c.item() < self.max_torque and u_c.item() >= 0:
                #     u_c = u_c.item()
                # elif u_c.item() > -self.max_torque and u_c.item() < 0:
                #     u_c = u_c.item()
                # elif u_c.item() >= self.max_torque:
                #     u_c = self.max_torque
                # elif u_c.item() <=  -self.max_torque:
                #     u_c = -self.max_torque
            empty_list.append(u_c)
        u = np.array(empty_list).reshape(-1,1)
        return u 
            
            


