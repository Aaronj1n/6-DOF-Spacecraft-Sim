import numpy as np 

class spacecraft:
    def __init__(self, controller, q, w, r, v,
                 I = np.diag([.13, .13, .08])):
        self.I = I
        self.controller = controller
        self.q = q
        self.w = w
        self.r = r
        self.v = v 




