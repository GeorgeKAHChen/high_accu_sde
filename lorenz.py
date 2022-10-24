import numpy as np
import mpmath as mp


dim = 3
para_dim = 3
rand_para_dim = 0
rand_para_len = 0

group_size = 1
init_vals = [[1, 1, 4]]
para_vals = [[10, 28, 8/3]]
rand_para_vals = [[]]

def f(curr_x, curr_para, curr_t, delta_t):
    x, y, z = curr_x[0], curr_x[1], curr_x[2]
    s, r, b = curr_para[0], curr_para[1], curr_para[2]

    x_new = - s * x + s * y
    y_new = - x * z + r * x - y
    z_new = x * y - b * z

    return [x_new, y_new, z_new]

def sf():
    return 

def Jf(curr_x, curr_para, curr_t, delta_t):
    x, y, z = curr_x[0], curr_x[1], curr_x[2]
    s, r, b = curr_para[0], curr_para[1], curr_para[2]

    return [mp.mpf(float(1.0))-s * delta_t,   s * delta_t, mp.mpf(float(0.0)),
            (r-z) * delta_t, mp.mpf(float(1.0)) - delta_t, -x * delta_t,
            y * delta_t,     x * delta_t, mp.mpf(float(1.0))-b * delta_t]

def psf():
    return 
