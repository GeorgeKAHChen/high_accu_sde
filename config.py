#------------------------------------------------
#
#      System parameters
#
#
#------------------------------------------------

from model import fdrs as dyn
accu = 45
delta_t = 1e-5
delta_t_ob = 1e-3
t_ob = 0.5
t_le = 0.9
t_ps = 0.8
t_max = 300

o_ob = 1
o_le = 0
o_ps = 0

#------------------------------------------------
#
#      Model Initial
#
#
#------------------------------------------------
import numpy as np

f = dyn.f
sf = dyn.sf
Jf = dyn.Jf
psf = dyn.psf

dim = dyn.dim
para_dim = dyn.para_dim
rand_para_dim = dyn.rand_para_dim
rand_para_len = dyn.rand_para_len

init_vals = dyn.init_vals
para_vals = dyn.para_vals
rand_para_vals = dyn.rand_para_vals

group_size = dyn.group_size
#------------------------------------------------
