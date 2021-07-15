#!/usr/bin/python

# ==============
# import section
# ==============

import sys

#for os interaction:
import os
# for math:
import numpy as np
from cmath import *
#for parallelization:
import multiprocessing
#for timeing:
from timeit import default_timer as timer
#for matlab output:
import scipy.io as sio
#for parsing (commandline input) arguments
import argparse
#singal handling, e.g. keyboard interrupts
import signal
from matplotlib import pyplot as plt


# ======================
# import section for fmm
# ======================

sys.path.append('/users/stud/ledwon/Documents/Bachelorarbeit/Magnus2D/fmm_source/')
import magnus_fmm1d_py as fmm1d

####### helping functions ######################################################
def zero_func(x):
  return np.zeros(x.shape)
def one_function(x):
  return np.ones(x.shape)
def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))

####### parameters #############################################################
wavelength = 5.5
#wavelength = np.arange(2,10.025,0.025)
numpts = 16
#numpts = np.array([16,22,32,46,64,90,128])#,182,256])
intervallVec = np.array([[0.0,5.],[5.,6.],[6., 8.],[8.,13.]])
setup_period = 4.0
incAngle = 0
#incAngle = np.arange(0,90,1.0)
glassblock_width = 0.5*setup_period
magnus_order = 0
period = setup_period
kztolfac = 1e-10
omega = 2.*pi/wavelength
kappa = omega*sin(pi*incAngle/180)#*sqrt(2.25)*(1.0+1e-15j)

####### functions ##############################################################
def input_pw(x):
      return one_function(x)*np.exp( 1.0j * kappa *x)
eps_glass = 2.25
eps_air = 1.
def eps_func(x,z):
        return one_function(x)
def eps_func_glass_struct(x,z):
        return one_function(x)#+ box_func(x,0,glassblock_width)*(eps_glass-eps_air)
def eps_func_glass(x,z):
        return one_function(x)*eps_glass
def mu_func(x,z):
        return one_function(x)
matmu = fmm1d.Material_struct(mu_func,mu_func,mu_func)
mateps = fmm1d.Material_struct(eps_func,eps_func,eps_func)
mateps_glass_struct = fmm1d.Material_struct(eps_func_glass_struct,eps_func_glass_struct,eps_func_glass_struct)
mateps_glass = fmm1d.Material_struct(eps_func_glass,eps_func_glass,eps_func_glass)
eps_list = [mateps,mateps_glass_struct,mateps_glass,mateps]

lstack = fmm1d.Layer_stack(numpts,kappa,omega,period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})
for i in range(0,len(intervallVec)):
  lstack.add_single_layer(eps_list[i],matmu,intervallVec[i],magnus_order)

#air
xvec = np.arange(0,6.1,0.1)
zvec = np.arange(0,5.1,0.1)
ex_air,ey_air,ez_air,hx_air,hy_air,hz_air = lstack.get_field_constant_single_layer(0,xvec,zvec)
ex,ey,ez,hx,hy,hz = lstack.get_all_fields(xvec)

#glassblock
xvec_gb = np.arange(0,6.1,0.1)
zvec_gb = np.arange(5.0,6.1,0.1)
ex_glassblock,ey_glassblock,ez_glassblock,hx_glassblock,hy_glassblock,hz_glassblock = lstack.get_field_constant_single_layer(1,xvec_gb,zvec_gb)

#glasslayer
xvec_g = np.arange(0,6.1,0.1)
zvec_g = np.arange(6.0,8.1,0.1)
ex_glass,ey_glass,ez_glass,hx_glass,hy_glass,hz_glass = lstack.get_field_constant_single_layer(2,xvec_g,zvec_g)

#air2
xvec_2 = np.arange(0,6.1,0.1)
zvec_2 = np.arange(8.,13.1,0.1)
ex_air2,ey_air2,ez_air2,hx_air2,hy_air2,hz_air2 = lstack.get_field_constant_single_layer(3,xvec_2,zvec_2)

ey = np.zeros((len(zvec)+len(zvec_gb)+len(zvec_g)+len(zvec_2),(len(xvec))),complex)
ey[0:len(zvec),:] = ey_air.transpose()
ey[len(zvec):len(zvec)+len(zvec_gb),:] = ey_glassblock.transpose()
ey[len(zvec)+len(zvec_gb):len(zvec)+len(zvec_gb)+len(zvec_g),:] = ey_glass.transpose()
ey[len(zvec)+len(zvec_gb)+len(zvec_g):len(zvec)+len(zvec_gb)+len(zvec_g)+len(zvec_2),:] = ey_air2.transpose()

plt.figure('Magnus Feldverteilung')
plt.imshow(-ey.real)
plt.colorbar()
plt.show()
