def rk4(curr_x, curr_para, curr_t, f, delta_t):
    x1 = curr_x
    t1 = curr_t
    k1 = f(x1, curr_para, t1, delta_t)
    
    x2 = []
    for i in range(0, len(curr_x)):
        x2.append(curr_x[i] + k1[i] * delta_t / 2)
    t2 = curr_t + delta_t / 2
    k2 = f(x2, curr_para, t2, delta_t)

    x3 = []
    for i in range(0, len(curr_x)):
        x3.append(curr_x[i] + k2[i] * delta_t / 2)
    t3 = curr_t + delta_t / 2
    k3 = f(x3, curr_para, t3, delta_t)

    x4 = []
    for i in range(0, len(curr_x)):
        x4.append(curr_x[i] + k3[i] * delta_t)
    t4 = curr_t + delta_t
    k4 = f(x4, curr_para, t4, delta_t)

    x_final = [] 
    for i in range(0, len(curr_x)):
        x_final.append(curr_x[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * delta_t / 6)

    return x_final
