import numpy as np 

def q_to_dcm(quaternion):
#please give quaternion as a 4x1 column np array
    beta_0 = quaternion[0] #the scalar part 
    beta_1 = quaternion[1] #i vector
    beta_2 = quaternion[2] #j vector
    beta_3 = quaternion[3] #k vector
    dcm = np.array[[(np.square(beta_0)+np.square(beta_1)-np.square(beta_2)-np.square(beta_3)), 
                     2*(beta_1*beta_2 + beta_0*beta_3), 2*(beta_1*beta_3-beta_0*beta_2)], 
                    [2*(beta_1*beta_2-beta_0*beta_3), np.square(beta_0)-np.square(beta_1)+np.square(beta_2)-np.square(beta_3),
                     2*(beta_2*beta_3+beta_0*beta_1)], [2*(beta_1*beta_3+beta_0*beta_2), 2*(beta_2*beta_3-beta_0*beta_1),
                                                        np.square(beta_0)-np.square(beta_1)-np.square(beta_2)+np.square(beta_3)]]
    return dcm
def dcm_to_q(dcm):
    #uses Shepperd's Selection Algorithm (https://motoq.github.io/doc/tnotes/dcmq.pdf)
    [beta_0, beta_1, beta_2, beta_3] = [0,0,0,0]
    trace = np.trace(dcm)
    c11 = dcm[0,0]
    c22 = dcm[1,1]
    c33 = dcm[2,2]
    c12 = dcm[0,1]
    c13 = dcm[0,2]
    c21 = dcm[1,0]
    c23 = dcm[1,2]
    c31 = dcm[2,0]
    c32 = dcm[2,1]
    if trace>c11 and trace>c22 and trace>c33:
        #eq 15
        beta_0 = np.sqrt((1+c11+c22+c33)/4)
        beta_1 = (c23 - c32)/(4*beta_0)
        beta_2 = (c31-c13)/(4*beta_0)
        beta_3 = (c12-c21)/(4*beta_0)
    elif c11>c22 and c11>c33:
        #eq 16
        beta_1 = np.sqrt((1+c11-c22-c33)/4)
        beta_0 = (c23-c32)/(4*beta_1)
        beta_2 = (c12+c21)/(4*beta_1)
        beta_3 = (c31+c13)/(4*beta_1)
    elif c22>c33:
        #eq 17 
        beta_2 = np.sqrt((1-c11+c22-c33)/4)
        beta_0 = (c31-c13)/(4*beta_2)
        beta_1 = (c12+c21)/(4*beta_2)
        beta_3 = (c23+c32)/(4*beta_2)
    else:
        #eq 18
        beta_3 = np.sqrt((1-c11-c22+c33)/4)
        beta_0 = (c12-c21)/(4*beta_3)
        beta_1 = (c31+c13)/(4*beta_3)
        beta_2 = (c23+c32)/(4*beta_3)

    return np.array([beta_0, beta_1, beta_2, beta_3]) #quaternion

def create_skew(vector): #vector is 3x1 column np array
    x = vector[0]
    y = vector[1]
    z = vector[2]
    skew = np.array([[0, -z, y], [z,0,-x], [-y, x, 0]])
    return skew

def quaternion_multiply(Q0,Q1):
    """
    Multiplies two quaternions.
 
    Input
    :param Q0: A 4 element array containing the first quaternion (q01,q11,q21,q31) 
    :param Q1: A 4 element array containing the second quaternion (q02,q12,q22,q32) 
 
    Output
    :return: A 4 element array containing the final quaternion (q03,q13,q23,q33) 
 
    """
    #credit to https://automaticaddison.com/how-to-multiply-two-quaternions-together-using-python/

    # Extract the values from Q0
    w0 = Q0[0]
    x0 = Q0[1]
    y0 = Q0[2]
    z0 = Q0[3]
     
    # Extract the values from Q1
    w1 = Q1[0]
    x1 = Q1[1]
    y1 = Q1[2]
    z1 = Q1[3]
     
    # Computer the product of the two quaternions, term by term
    Q0Q1_w = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1
    Q0Q1_x = w0 * x1 + x0 * w1 + y0 * z1 - z0 * y1
    Q0Q1_y = w0 * y1 - x0 * z1 + y0 * w1 + z0 * x1
    Q0Q1_z = w0 * z1 + x0 * y1 - y0 * x1 + z0 * w1
     
    # Create a 4 element array containing the final quaternion
    final_quaternion = np.array([Q0Q1_w, Q0Q1_x, Q0Q1_y, Q0Q1_z])
     
    # Return a 4 element array containing the final quaternion (q02,q12,q22,q32) 
    return final_quaternion


#please give dcm as a 4x4 np array

