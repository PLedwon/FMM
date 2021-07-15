#!/usr/bin/python

# ==============
# import section
# ==============

#path_to_fmm_src = '/users/stud/kiel/diss_works/beamsplitter/fmm_source/'

import sys

#for os interaction:
import os
# for math:
import numpy as np
from cmath import *



#for matlab output:
import scipy.io as sio
#for plotting
#from matplotlib import pyplot as plt

# ======================
# import section for fmm
# ======================

sys.path.append('/users/stud/ledwon/Documents/Bachelorarbeit/Magnus2D/fmm_source/')
import fmm1d_py as fmm1d



# =================
# helping functions
# =================

def zero_func(x):
  return np.zeros(x.shape)

def one_function(x):
  return np.ones(x.shape)

def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))

# =========
# init data
# =========

wavelength = np.arange(5.5,6.5,1)
#numpts = np.array([16,22,32,46,64,90,128,182,256,362,512])#,724,1024])
numpts = 4
glass_thickness = np.array([1])
substrate_thickness = np.array([3])
setup_period = np.array([(5)])
glass_width = 2.5
incAngle = np.arange(30,60,30)
eps_glass = 1.
eps_air = 1.


# =========
# filenames
# =========

save_filename = 'glass_grating_air.mat'

# =========

omega = 2*pi/wavelength

TEy = np.zeros([len(omega),len(incAngle)])
REy = np.zeros([len(omega),len(incAngle)])
THy = np.zeros([len(omega),len(incAngle)])
RHy = np.zeros([len(omega),len(incAngle)])

TEy_ZO = np.zeros([len(omega),len(incAngle)])
REy_ZO = np.zeros([len(omega),len(incAngle)])
THy_ZO = np.zeros([len(omega),len(incAngle)])
RHy_ZO = np.zeros([len(omega),len(incAngle)])

for incAng_idx in range(0,len(incAngle)):
    for omega_idx in range(0,len(omega)):

        kappa = omega[omega_idx]*sin(pi*incAngle[incAng_idx]/180)#*sqrt(2.25)*(1.0+1e-15j)
        kztolfac = 1e-10

        layer_thickness = np.array([5, glass_thickness, substrate_thickness, 5])

        def eps1_func(x):
            return eps_air*one_function(x)
        def eps2_func(x):
            return eps_air*one_function(x) + (eps_glass- eps_air)*box_func(np.mod(x,setup_period),0.0,glass_width)
        def eps3_func(x):
            return eps_glass*one_function(x)


        def mu_func(x):
            return one_function(x)

        def input_pw(x):
            return one_function(x)*np.exp( 1.0j * kappa *x)

        matmu1 = fmm1d.Material_struct(mu_func,mu_func,mu_func)
        mateps1 = fmm1d.Material_struct(eps1_func,eps1_func,eps1_func)
        mateps2 = fmm1d.Material_struct(eps2_func,eps2_func,eps2_func)
        mateps3 = fmm1d.Material_struct(eps3_func,eps3_func,eps3_func)


        eps_list = [mateps1,mateps2,mateps3,mateps1]
        mu_list = [matmu1,matmu1,matmu1,matmu1]

        lstack = fmm1d.Layer_stack(numpts,kappa,omega[omega_idx],setup_period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})
        for i in range(0,len(layer_thickness)):
            lstack.add_single_layer(layer_thickness[i],eps_list[i],mu_list[i])

        (p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy) = lstack.get_power_data()
        difforder_data = lstack.get_power_data_by_orders()


        TEy_tmp = pEndscatEy/p0incEy
        REy_tmp = - p0scatEy/p0incEy
        THy_tmp = pEndscatHy/p0incHy
        RHy_tmp = - p0scatHy/p0incHy

        do_idx = difforder_data.order_number.index(0)

        TEy_ZO_tmp = difforder_data.pEndscatEy[do_idx]/p0incEy
        REy_ZO_tmp = - difforder_data.p0scatEy[do_idx]/p0incEy
        THy_ZO_tmp = difforder_data.pEndscatHy[do_idx]/p0incHy
        RHy_ZO_tmp = - difforder_data.p0scatHy[do_idx]/p0incHy

        TEy[omega_idx,incAng_idx] = TEy_tmp
        REy[omega_idx,incAng_idx] = REy_tmp
        THy[omega_idx,incAng_idx] = THy_tmp
        RHy[omega_idx,incAng_idx] = RHy_tmp

        TEy_ZO[omega_idx,incAng_idx] = TEy_ZO_tmp
        REy_ZO[omega_idx,incAng_idx] = REy_ZO_tmp
        THy_ZO[omega_idx,incAng_idx] = THy_ZO_tmp
        RHy_ZO[omega_idx,incAng_idx] = RHy_ZO_tmp

#Feldverteilung berechnen
xvec = np.arange(0,5.01,0.01)
zvec = np.arange(0,14.01,0.01)
(ex,ey,ez,hx,hy,hz)=lstack.get_modes_at(xvec,zvec)

#save to savefile
sio.matlab.savemat(save_filename,
              {"TEy": TEy, "REy": REy, "THy": THy, "RHy": RHy,
                "TEy_ZO": TEy_ZO, "REy_ZO": REy_ZO, "THy_ZO": THy_ZO, "RHy_ZO": RHy_ZO,
                "wavelength": wavelength, "omega": omega, "incAngle": incAngle#, "ex":ex, "ey":ey, "ez":ez, "xvec": xvec, "zvec":zvec
                })
