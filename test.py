
import sys
import os
import numpy as np
from cmath import *
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
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
circleradius = 1.5
def heaviside_circle(x,z):
    heavisidevec=np.zeros(np.size(x))
    for i in range(0,len(x)):
        if np.power(x[i],2)+np.power(z,2) <= np.power(circleradius,2):
            heavisidevec[i]=1
        else:
            heavisidevec[i]=0
    return heavisidevec

####### generate staircase #############################################################
def smooth_step_function_staircase(N,w,h,r1,r2):
  if N < 1:
    raise ValueError('N must be at least 1')

  x = np.zeros(N+1)
  for i in range(0,len(x)):
    x[i] = i*1.0*w/N

  y = np.zeros(len(x))

  #if r1*2 > w:
    #raise ValueError('curvature radius r1 must be smaller or equal to half the width.')
  #if r2*2 > w:
    #raise ValueError('curvature radius r2 must be smaller or equal to half the width.')

  if r1 + r2 > w:
    raise ValueError('sum of curvature radii must be smaller or equal to the width')

  if r1+r2 <= h:
    # determine midpoints line and tangential line intersection point
    length_of_midpt_line = np.sqrt(w**2 + (h - r1 - r2)**2)
    alpha = np.arcsin((r1+r2)/length_of_midpt_line)
    beta = np.pi/2.0 - alpha
    gamma = np.arcsin(w/length_of_midpt_line)
    delta = np.pi - beta - gamma
  else:
    length_of_midpt_line = np.sqrt(w**2 + (h - r1 - r2)**2)
    alpha = np.arcsin((r1+r2)/length_of_midpt_line)
    beta = np.pi/2.0 - alpha
    gamma = np.arccos(w/length_of_midpt_line)
    delta = np.pi/2.0 - beta - gamma

  x1 = r1 *np.sin(delta)
  x2 = w - r2 * np.sin(delta)

  y1 = h - r1 + np.sqrt(r1**2 - x1**2)
  y2 = r2 - np.sqrt(r2**2 - (x2 - w)**2)

  for i in range(0,len(x)):
    if x[i] <= 0.:
      y[i] = h
    elif x[i] < x1:
      y[i] = h - r1 + np.sqrt(r1**2 - x[i]**2)
    elif x[i] < x2:
      y[i] = (y1 - y2)/(x1-x2) *(x[i] - x2) + y2
    elif x[i] < w:
      y[i] = r2 - np.sqrt(r2**2 - (x[i] - w)**2)
    else:
      y[i] = 0.0

  heights = np.zeros(N)
  for i in range(0,len(y)-1):
    heights[i] = (y[i] + y[i+1])*0.5


  xvals = np.zeros(N)
  for i in range(0,len(heights)):
    if heights[i] > y1:
      xvals[i] = np.sqrt(r1**2 - (heights[i] - h +r1)**2)
    elif heights[i] > y2:
      xvals[i] = (x1-x2)/(y1-y2) *(heights[i] - y2) + x2
    else:
      xvals[i] = w - np.sqrt(r2**2 - (heights[i] - r2)**2)

  thicknesses = np.zeros(N)
  for i in range(0,len(y)-1):
    thicknesses[i] = np.abs(y[i] - y[i+1])

  pltx = []
  plty = []

  pltx.append(0.0)
  plty.append(y[0])

  for i in range(0,len(xvals)):
    pltx.append(xvals[i])
    pltx.append(xvals[i])
    plty.append(y[i])
    plty.append(y[i+1])

  pltx.append(w)
  plty.append(y[-1])

  pltx = np.array(pltx)
  plty = np.array(plty)

  return xvals,heights,x,y,thicknesses,pltx,plty

####### parameters #############################################################
wavelength = 5.5
numpts = 32
N = numpts

setup_period = 6.
centerpoint = setup_period*0.5
incAngle = 0
omega = 2*pi/wavelength
kappa = omega*sin(pi*incAngle/180)
kztolfac = 1e-10
magnus_order = 1
Ns = 10

xvals,heights,x,y,thicknesses,pltx,plty = smooth_step_function_staircase(Ns,circleradius,circleradius,circleradius,0)


####### functions ##############################################################

intervallVec = np.zeros((Ns,2),float)
eps_metal = 2 + 10j
eps_air = 1.
air_thickness = 5.

intervallVec[0,0] = air_thickness
intervallVec[0,1] = air_thickness + thicknesses[0]

for i in range(1,Ns):
    intervallVec[i,0] = intervallVec[i-1,1]
    intervallVec[i,1] = intervallVec[i,0] + thicknesses[i]

#intervallVec[Ns-1,0] = intervallVec[Ns-2,1]
#intervallVec[Ns-1,1] = intervallVec[Ns-1,0] + xvals[Ns-1]-xvals[Ns-2]


def input_pw(x):
      return one_function(x)*np.exp( 1.0j * kappa *x)

def eps_func(x,z):
    return one_function(x)

def mu_func(x,z):
    return one_function(x)

matmu = fmm1d.Material_struct(mu_func,mu_func,mu_func)
mateps = fmm1d.Material_struct(eps_func,eps_func,eps_func)
lstack = fmm1d.Layer_stack(numpts,kappa,omega,setup_period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})


#single layer of air
intervall_air = np.array([0.0, air_thickness ])
lstack.add_single_layer(mateps,matmu,intervall_air,0)



# semicircle
for i in range(0,len(intervallVec)):
    def eps_func_semicircle(x,z):
        return one_function(x) + heaviside_circle(x-centerpoint,z-air_thickness-circleradius)*(eps_metal-eps_air) #box_func(x,centerpoint-xvals[i],centerpoint+xvals[i])*(eps_metal-eps_air)#

    mateps_semicircle = fmm1d.Material_struct(eps_func_semicircle,eps_func_semicircle,eps_func_semicircle)
    lstack.add_single_layer(mateps_semicircle,matmu,intervallVec[i],magnus_order)

intervall_air = np.array([air_thickness+np.sum(thicknesses), 2*air_thickness+circleradius])
lstack.add_single_layer(mateps,matmu,intervall_air,0)

#xvec = np.arange(0,setup_period,0.02)

#ex,ey,ez,hx,hy,hz = lstack.get_all_fields(xvec)
p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy = lstack.get_power_data()

print p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy



TEy = pEndscatEy/p0incEy
REy = p0scatEy/p0incEy
THy = pEndscatHy/p0incHy
RHy = p0scatHy/p0incHy

print 'M:',TEy, REy,THy,RHy

xvec = np.arange(0,6.1,0.1)
zvec = np.arange(0,5.1,0.1)

#ex,ey,ez,hx,hy,hz = lstack.get_field_constant_single_layer(0,xvec,zvec)
#plt.figure(1)
#plt.imshow(hy.real)
#plt.colorbar()
#plt.show()

#TMat=np.exp(1j*np.tensordot(xvec,(lstack.layer_list[5].kxvec + lstack.kappa),0))
#HyModes = np.dot(TMat,lstack.layer_list[5].hycoeffs[N:2*N,:])
#EyModes = np.dot(TMat,lstack.layer_list[5].eycoeffs[N:2*N,:])
#plt.figure('Ey Modes')
#plt.plot(xvec,EyModes)
#plt.figure('Hy Modes')
#plt.plot(xvec,HyModes)
#plt.show()



#plt.figure('M Ex')
#plt.title('Ex real')
#plt.plot(xvec, ex.real) #xvec, ey.real,xvec, ez.real, xvec, hx.real,xvec, hy.real, xvec, hz.real)
#plt.figure('M Hx')
#plt.title('Hx real')
#plt.plot(xvec, hx.real)



#plt.figure(2)
#plt.plot(xvec, ez.real,label='Ez M')
#plt.legend()

#plt.figure('EV M Ey')
#plt.plot(lstack.layer_list[0].kzvalsEy.real[0:N],lstack.layer_list[0].kzvalsEy.imag[0:N],'.')

k=2

plt.figure('EV M Ey')
plt.plot(lstack.layer_list[k].kzvalsEy.real[0:N],lstack.layer_list[k].kzvalsEy.imag[0:N],'.b',lstack.layer_list[k].kzvalsEy.real[N:2*N],lstack.layer_list[k].kzvalsEy.imag[N:2*N],'.r')
plt.plot(np.array([0,0]),np.array([np.min(lstack.layer_list[k].kzvalsEy.imag),np.max(lstack.layer_list[k].kzvalsEy.imag)]),'k-')
plt.plot(np.array([np.min(lstack.layer_list[k].kzvalsEy.real),np.max(lstack.layer_list[k].kzvalsEy.real)]),np.array([0,0]),'k-')

plt.figure('EV M Hy')
plt.plot(lstack.layer_list[k].kzvalsHy.real[0:N],lstack.layer_list[k].kzvalsHy.imag[0:N],'.b',lstack.layer_list[k].kzvalsHy.real[N:2*N],lstack.layer_list[k].kzvalsHy.imag[N:2*N],'.r')
plt.plot(np.array([0,0]),np.array([np.min(lstack.layer_list[k].kzvalsHy.imag),np.max(lstack.layer_list[k].kzvalsHy.imag)]),'k-')
plt.plot(np.array([np.min(lstack.layer_list[k].kzvalsHy.real),np.max(lstack.layer_list[k].kzvalsHy.real)]),np.array([0,0]),'k-')
plt.show()

save_path = '/users/stud/ledwon/Documents/Bachelorarbeit/Magnus2D/Testsysteme/Spektren/'
save_filename = 'evM1'
sio.matlab.savemat(save_path +  save_filename,
            {"EV_TE_M1": lstack.layer_list[k].kzvalsEy,"EV_TM_M1": lstack.layer_list[k].kzvalsHy})
