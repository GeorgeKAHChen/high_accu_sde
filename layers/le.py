from copy import deepcopy
import math
import mpmath as mp


def le(curr_x, curr_para, curr_t, Jf, delta_t, eye, val_le, t_le, dim):
    jacobian = Jf(curr_x, curr_para, curr_t, delta_t)
    
    """matrix times"""
    mat_result = [mp.mpf(float(0.0)) for n in range(dim * dim)]
    for i in range(0, dim):
        for j in range(0, dim):
            mat_result[i * dim + j] = mp.mpf(float(0.0));
            for k in range(0, dim):
                mat_result[i * dim + j] += jacobian[i * dim + k] * eye[j + k * dim];
    eye = deepcopy(mat_result)
    #print(eye)
    """gram_schmidt"""
    for kase in range(0, dim):
        for i in range(0, kase):
            inner_beta = mp.mpf(float(0.0))
            inner_ab = mp.mpf(float(0.0))
            for j in range(0, dim):
                inner_beta += eye[i+j*dim] * eye[i+j*dim]
                inner_ab += eye[i+j*dim] * mat_result[kase+j*dim]
            for j in range(0, dim):
                eye[kase + j*dim] -= (inner_ab/inner_beta) * eye[i+j*dim]
    #print(eye)
    """Normalization"""
    para_norm = [0 for n in range(dim)]
    for i in range(0, dim):
        for j in range(0, dim):
            para_norm[i] += math.pow(eye[i + dim*j], 2)
        para_norm[i] = math.sqrt(para_norm[i])          
        for j in range(0, dim):
            eye[i + dim*j] /= para_norm[i]
    #print(eye)
    """Calculate current LE"""
    new_spec = [0 for n in range(dim)]
    for i in range(0, dim):
        for j in range(0, dim):
            new_spec[i] += eye[i + j * dim] * mat_result[i + j * dim]
    #print(new_spec)
    """Calculate final LE"""
    time_minus = curr_t - t_le
    for i in range(0, len(val_le)):
        val_le[i] = (time_minus * val_le[i] + math.log(new_spec[i])) / (time_minus + delta_t)
    #print(val_le)
    #input()
    return val_le, eye



