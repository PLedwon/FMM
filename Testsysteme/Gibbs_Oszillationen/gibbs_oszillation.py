
import sys
import os
import numpy as np
from cmath import *
import scipy.linalg
#for matlab output:
import scipy.io as sio
sys.path.append('/users/stud/ledwon/Documents/Bachelorarbeit/Magnus2D/fmm_source/')
import magnus_fmm1d_py as fmm1d

####### helping functions ######################################################
def zero_func(x):
  return np.zeros(x.shape)
def one_function(x):
  return np.ones(x.shape)
def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))
def heaviside_circle(x,z):
    heavisidevec=np.zeros(np.size(x))
    for i in range(0,len(x)):
        if 1.*1.-np.power(x[i],2)-np.power(z,2) < 0:
            heavisidevec[i]=1
        else:
            heavisidevec[i]=0
    return heavisidevec


####### parameters #############################################################
wavelength = 5.5
numpts = np.array([16,22,32,46,64,90,182,256,362,512])
N = numpts
intervallVec = np.array([[0.,5.],[5.,6.],[6.,11.]])
setup_period = 4.
incAngle = 0
omega = 2*pi/wavelength
kappa = omega*sin(pi*incAngle/180)
kztolfac = 1e-10
magnus_order = 0

####### functions ##############################################################
def input_pw(x):
      return one_function(x)*np.exp( 1.0j * kappa *x)
      #return (one_function(x) + np.exp(2*np.pi*1j*x*(5.0/setup_period)))*np.exp( 1.0j * kappa *x)

eps_glass = 2.25
eps_air = 1.
def eps_func(x,z):
        return one_function(x)
def eps_func_glass(x,z):
        return one_function(x) + box_func(x,1,3)*(eps_glass-eps_air)
def mu_func(x,z):
        return one_function(x)

matmu = fmm1d.Material_struct(mu_func,mu_func,mu_func)
mateps = fmm1d.Material_struct(eps_func,eps_func,eps_func)
mateps_glass = fmm1d.Material_struct(eps_func_glass,eps_func_glass,eps_func_glass)

eps_list = [mateps,mateps_glass,mateps]
step = 1./512/3
xvec = np.arange(0,setup_period+step,step)

exMat = np.zeros((len(xvec),len(numpts)),complex)

for n_idx in range(0,len(numpts)):

    lstack = fmm1d.Layer_stack(numpts[n_idx],kappa,omega,setup_period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})
    for i in range(0,len(intervallVec)):
        lstack.add_single_layer(eps_list[i],matmu,intervallVec[i],magnus_order)

    ex,ey,ez,hx,hy,hz = lstack.get_all_fields(xvec)
    exMat[:,n_idx] = ex[:,1]

save_filename = 'gibbs_oszillation.mat'
sio.matlab.savemat(save_filename,
              {"wavelength": wavelength, "omega": omega, "incAngle": incAngle, "exMat": exMat, "numpts": numpts, "xvec": xvec, "step": step
                })
