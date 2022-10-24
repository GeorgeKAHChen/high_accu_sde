import numpy as np
from copy import deepcopy
import os

import mpmath as mp

import config
from layers import rk4
from layers import le


mp.mp.dps = config.accu

o_ob = config.o_ob
o_le = config.o_le
o_ps = config.o_ps

f = config.f
sf = config.sf
Jf = config.Jf
psf = config.psf

dim = config.dim
para_dim = config.para_dim
rand_para_dim = config.rand_para_dim
rand_para_len = config.rand_para_len

group_size = config.group_size

str_write_time = 1e6

def accu_arr(arr):
    for i in range(0, len(arr)):
        tmp_arr = []
        for j in range(0, len(arr[i])):
            tmp_arr.append(mp.mpf(float(arr[i][j])))
        arr[i] = deepcopy(tmp_arr)

    return arr


def main():
    init_vals = accu_arr(config.init_vals) 
    para_vals = accu_arr(config.para_vals) 
    rand_para_vals = accu_arr(config.rand_para_vals)

    delta_t = mp.mpf(float(config.delta_t))
    delta_t_ob = mp.mpf(float(config.delta_t_ob))

    t_max = mp.mpf(float(config.t_max))
    t_ob = mp.mpf(float(config.t_ob)) * t_max
    t_le = mp.mpf(float(config.t_le)) * t_max
    t_ps = mp.mpf(float(config.t_ps)) * t_max
    
    print(init_vals, para_vals, rand_para_vals) 

    output_le = []
    output_ps = []

    

    if not os.path.exists("output"):
        os.system("mkdir output")

    print(group_size)

    for kase in range(0, group_size):
        """
            Initialization and definition

        """
        curr_x = init_vals[kase]
        curr_para = para_vals[kase]
        curr_rand_para = rand_para_vals[kase]

        curr_t = mp.mpf(float(0))
        curr_t_ob = mp.mpf(float(0))

        ttl = 0

        val_le = []
        eye = []
        str_orbits = ""
        ps = []
        last_val = []
        last_last_val = []
        save_file_name = ""

        """
            Pre-treatment
        
        """
        if o_ob:
            save_file_name = "output/" + str(kase) + ".info"
            if os.path.exists(save_file_name):
                os.system("rm -rf " + save_file_name)
            
            file = open(save_file_name, "w")
            tmp_str = ""
            tmp_str += "curr_x" + "\n"
            tmp_str += str(curr_x) + "\n"
            tmp_str += "curr_para" + "\n"
            tmp_str += str(curr_para) + "\n"
            tmp_str += "curr_rand_para" + "\n"
            tmp_str += str(curr_rand_para) + "\n"
            tmp_str += "curr_rand_para" + "\n"
            tmp_str += str(curr_rand_para) + "\n"
            tmp_str += "config.accu, o_ob, o_le, o_ps, delta_t, t_ob, t_le, t_ps, t_max" + "\n"
            tmp_str += str([config.accu, o_ob, o_le, o_ps, delta_t, t_ob, t_le, t_ps, t_max]) + "\n"
            tmp_str += "dim, para_dim, rand_para_dim, rand_para_len" + "\n"
            tmp_str += str([dim, para_dim, rand_para_dim, rand_para_len]) + "\n"
            file.write(tmp_str)
            file.close()

            save_file_name = "output/" + str(kase) + ".dat"
            if os.path.exists(save_file_name):
                os.system("rm -rf " + save_file_name)
            file = open(save_file_name, "w")
            file.write("")
            file.close()

        if o_le:
            val_le = [mp.mpf(float(0.0)) for n in range(dim)]
            eye = [mp.mpf(float(0.0)) for n in range(dim * dim)]
            tmp_val = 0
            while 1:
                if tmp_val > dim * dim:
                    break
                eye[tmp_val] = mp.mpf(float(1.0))
                tmp_val += (dim + 1)

        """
            Main algorithm
        
        """
        while 1:
            ttl += 1 
            if ttl % 1e5 == 0:
                print(kase, group_size, float(curr_t), t_max, float(curr_x[0]))

            if curr_t > t_max:
                break
            curr_t += delta_t
            curr_t_ob += delta_t

            curr_x = rk4.rk4(curr_x, curr_para, curr_t, f, delta_t)

            if o_ob and curr_t > t_ob and curr_t_ob > delta_t_ob:
                curr_t_ob = 0
                str_orbits += str(curr_t) + " "
                for i in range(0, len(curr_x)):
                    str_orbits += str(curr_x[i]) + " "
                str_orbits += "\n"
                if ttl % str_write_time == 0:
                    #print(str_orbits)
                    file = open(save_file_name, "a")
                    file.write(str_orbits)
                    file.close()
                    str_orbits = ""

            if o_le and curr_t > t_le:
                val_le, eye = le.le(curr_x, curr_para, curr_t, Jf, delta_t, eye, val_le, t_le, dim)

        """
            After treatment
        
        """
        if o_ob: 
            file = open(save_file_name, "a")
            file.write(str_orbits)
            file.close()
            str_orbits = ""

        if o_le:
            save_file_name = "output/" + str(kase) + ".dat"
            file = open(save_file_name, "a")
            tmp_str = ""
            for i in range(0, len(val_le)):
                tmp_str += str(float(val_le[i])) + " "
            file.write(tmp_str)
            file.close()
            print(tmp_str)
    return 

if __name__ == '__main__':
    main()