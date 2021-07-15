
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

def init_worker():
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        print 'signal: %s' % signal_num
    signal.signal(signal.SIGINT, sig_int)

# =================
# parsing input
# =================

parser = argparse.ArgumentParser()
parser.add_argument('-j','--numproc'    ,type=int,help='number of cores used',action='store',metavar='N'        ,dest='num_proc',required=True)
parser.add_argument('-c','--calcnum'    ,type=int,help='run calculation number N',action='store',metavar='N'    ,dest='calc_idx',required=True)
parser.add_argument('-o','--savepath'   ,help='path to save the output file(s)',action='store'                  ,dest='save_path',required=True)
parser.add_argument('-l','--logpath'    ,help='path to save the log file (deleted after calculation)',action='store',dest='log_path',required=True)
parser.add_argument('-t','--logtime'    ,type=float,help='time interval for logoutput (sec)',action='store'     ,dest='log_time',default=120.)
parser.add_argument('-s','--sourcepath' ,help='path to the fmm1d_py packet',action='store'                      ,dest='source_path',required=True)

inputargs = parser.parse_args()


# ======================
# import section for fmm
# ======================
sys.path.append(inputargs.source_path)
import fmm1d_py as fmm1d


def zero_func(x):
  return np.zeros(x.shape)
def one_function(x):
  return np.ones(x.shape)
def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))



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

def get_box_and_layer_thickness(N,w,h,r1,r2):
    xvals,heights,x,y,thicknesses,pltx,plty = smooth_step_function_staircase(N,w,h,r1,r2)

    box_thickness = y[0:len(y)]
    box_thickness = box_thickness[np.arange(len(box_thickness)-2,-1,-1)]
    box_thickness[len(box_thickness)-1] = r1

    layer_thickness = np.zeros(len(xvals))
    for i in range(len(xvals)-1,0,-1):
        layer_thickness[i] =  xvals[i] - xvals[i-1]
    layer_thickness[0] = xvals[0]

    return box_thickness, layer_thickness
################################################################################

#wavelength = np.arange(2,10,2)
wavelength = np.arange(2,10.025,0.025)
numpts = np.array([16,22,32,46,64,90,128])
#numpts = np.array([16])
setup_period = np.array([6.0])
centerpoint = setup_period/2
#incAngle = np.array([0])
incAngle = np.arange(0,90,1.0)
omega = 2.*pi/wavelength
kztolfac = 1e-10
eps_air = 1.
eps_glass = 2.25
circleradius = 1.5


#number of layers
NsVec = np.array([1,2,4,8,12,16])
#NsVec = np.array([1,2])

num_of_runs = len(NsVec)
if inputargs.calc_idx >= num_of_runs:
  raise IndexError('There is no' + str(inputargs.calc_idx) + '-th calculation to be run out of' + str(num_of_runs) + ' calculations!')

Ns = NsVec[inputargs.calc_idx]

log_filname = 'metallsubstrate'+str(inputargs.calc_idx)+'.log'
save_filname = 'metallsubstrate'+str(inputargs.calc_idx)+'.mat'

if not os.path.exists(inputargs.log_path):
  os.makedirs(inputargs.log_path)

box_thickness, layer_thickness = get_box_and_layer_thickness(Ns,circleradius,circleradius,circleradius,0)

time_needed = np.zeros([len(omega),len(incAngle),len(numpts)])

TEy = np.zeros([len(omega),len(incAngle),len(numpts)])
REy = np.zeros([len(omega),len(incAngle),len(numpts)])
THy = np.zeros([len(omega),len(incAngle),len(numpts)])
RHy = np.zeros([len(omega),len(incAngle),len(numpts)])

def run_func(omega_idx,incAng_idx):
  TEy_tmp = np.zeros([len(numpts)])
  REy_tmp = np.zeros([len(numpts)])
  THy_tmp = np.zeros([len(numpts)])
  RHy_tmp = np.zeros([len(numpts)])

  time_needed_tmp = np.zeros([len(numpts)])


  for j in range(0,len(numpts)):
    numpts_tmp = numpts[j]
    kappa = omega[omega_idx]*sin(pi*incAngle[incAng_idx]/180)#*sqrt(2.25)*(1.0+1e-15j)
    period = setup_period
    kztolfac = 1e-10
    eps_metal = eps_ag_vals[omega_idx]

    def input_pw(x):
        return one_function(x)*np.exp( 1.0j * kappa *x)

    def eps_func(x):
        return one_function(x)
    def mu_func(x):
        return one_function(x)
    #glass
    def eps_func2(x):
        return one_function(x)*eps_glass


    matmu = fmm1d.Material_struct(mu_func,mu_func,mu_func)
    mateps = fmm1d.Material_struct(eps_func,eps_func,eps_func)
    mateps2 = fmm1d.Material_struct(eps_func2,eps_func2,eps_func2)


    starttime1 = timer()
    lstack = fmm1d.Layer_stack(numpts_tmp,kappa,omega[omega_idx],period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})

    #single layer of air
    lstack.add_single_layer(5,mateps,matmu)

    # add semicircle
    for i in range(0,len(layer_thickness)):
        def eps_func_semicircle(x):
                return one_function(x) + box_func(x,centerpoint-box_thickness[i],centerpoint+box_thickness[i])*(eps_glass-eps_air)
        mateps_semicircle = fmm1d.Material_struct(eps_func_semicircle,eps_func_semicircle,eps_func_semicircle)
        lstack.add_single_layer(layer_thickness[i],mateps_semicircle,matmu)

    #single layer of glass
    lstack.add_single_layer(2,mateps2,matmu)
    #single layer of air
    lstack.add_single_layer(5,mateps,matmu)

    (p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy) = lstack.get_power_data()

    stoptime = timer()

    time_needed_tmp[j] = stoptime-starttime1

    TEy_tmp[j] = pEndscatEy/p0incEy
    REy_tmp[j] = - p0scatEy/p0incEy
    THy_tmp[j] = pEndscatHy/p0incHy
    RHy_tmp[j] = - p0scatHy/p0incHy
#    print 'box:',box_thickness,'layer:', layer_thickness
  return (TEy_tmp,REy_tmp,THy_tmp,RHy_tmp,time_needed_tmp)

if __name__== '__main__':
  pool = multiprocessing.Pool(inputargs.num_proc,init_worker)
  total_calc_num = len(incAngle)*len(omega)
  try:
    tstart = timer()
    time_idx = 0
    tot_calc_idx = 0
    for incAng_idx in range(0,len(incAngle)):
      res = []
      for omega_idx in range(0,len(omega)):
          #print run_func(omega_idx,incAng_idx)
          res.append(pool.apply_async(run_func,(omega_idx,incAng_idx,)))
      res_idx = 0

      for omega_idx in range(0,len(omega)):
        (TEy_tmp,REy_tmp,THy_tmp,RHy_tmp,time_needed_tmp) = res[res_idx].get()
        res_idx = res_idx + 1
        tot_calc_idx = tot_calc_idx + 1
        TEy[omega_idx,incAng_idx,:] = TEy_tmp
        REy[omega_idx,incAng_idx,:] = REy_tmp
        THy[omega_idx,incAng_idx,:] = THy_tmp
        RHy[omega_idx,incAng_idx,:] = RHy_tmp
        time_needed[omega_idx,incAng_idx,:] = time_needed_tmp
        tstop_curr = timer()
        if int((tstop_curr - tstart ) / inputargs.log_time) > time_idx:
            time_idx = int((tstop_curr - tstart ) / inputargs.log_time)
            #save to logfile
            flog = open(inputargs.log_path + log_filname,'a')
            flog.write("Status: %(indx)i / %(num_of_calcs)i calcs run, i.e. %(percentage)2.3f %%; needed time: %(time_need).1f min, ETA in: %(time_eta).1f min \n" %{"indx": tot_calc_idx, "num_of_calcs": total_calc_num,"percentage": tot_calc_idx*100./float(total_calc_num), "time_need": (tstop_curr - tstart)/60.,"time_eta": (tstop_curr - tstart)/60. * ((float(total_calc_num)/tot_calc_idx) - 1.0) })
            flog.close()


    #save to savefile
    sio.matlab.savemat(inputargs.save_path +  save_filname,
                {"TEy": TEy, "REy": REy, "THy": THy, "RHy": RHy,
                   "wavelength": wavelength, "omega": omega, "incAngle": incAngle,
                   "numpts": numpts, "setup_period": setup_period, "time_needed": time_needed})
    #delete logfile
    if os.path.isfile(inputargs.log_path + log_filname):
        os.system('rm -f ' + inputargs.log_path + log_filname)

  except KeyboardInterrupt:
    print "Caught KeyboardInterrupt, terminating workers"
    pool.terminate()
    pool.join()
