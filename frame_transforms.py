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


#please give dcm as a 4x4 np array

