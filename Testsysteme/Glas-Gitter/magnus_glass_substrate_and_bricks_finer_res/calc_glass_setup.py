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

def init_worker():
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        print 'signal: %s' % signal_num
    signal.signal(signal.SIGINT, sig_int)

#from matplotlib import pyplot as plt

# =================
# parsing input
# =================

parser = argparse.ArgumentParser()
parser.add_argument('-j','--numproc'    ,type=int,help='number of cores used',action='store',metavar='N'        ,dest='num_proc',required=True)
#parser.add_argument('-c','--calcnum'    ,type=int,help='run calculation number N',action='store',metavar='N'    ,dest='clac_idx',required=True)
parser.add_argument('-o','--savepath'   ,help='path to save the output file(s)',action='store'                  ,dest='save_path',required=True)
parser.add_argument('-l','--logpath'    ,help='path to save the log file (deleted after calculation)',action='store',dest='log_path',required=True)
parser.add_argument('-t','--logtime'    ,type=float,help='time interval for logoutput (sec)',action='store'     ,dest='log_time',default=120.)
parser.add_argument('-s','--sourcepath' ,help='path to the fmm1d_py packet',action='store'                      ,dest='source_path',required=True)

inputargs = parser.parse_args()

# ======================
# import section for fmm
# ======================

sys.path.append(inputargs.source_path)
import magnus_fmm1d_py as fmm1d

####### helping functions ######################################################
def zero_func(x):
  return np.zeros(x.shape)
def one_function(x):
  return np.ones(x.shape)
def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))

####### parameters #############################################################
wavelength = np.arange(2,10.008,0.008)
#wavelength = np.arange(2,10.025,0.025)
numpts = np.array([16,22,32,46,64,90,128,182,256])
#numpts = np.array([16,22,32,46,64,90,128,182,256])
intervallVec = np.array([[0.0,5.],[5.,6.],[6., 8.],[8.,13.]])
setup_period = 4
incAngle = np.arange(0,90,0.5)
#incAngle = np.arange(0,90,1.0)
glassblock_width = 0.5*setup_period
magnus_order = 3

####### filenames ##############################################################
log_filname = 'magnus_glasssystem.log'
save_filname = 'magnus_glassystem.mat'

if not os.path.exists(inputargs.log_path):
  os.makedirs(inputargs.log_path)

# =========

omega = 2.*pi/wavelength


time_needed = np.zeros([len(omega),len(incAngle),len(numpts)])

TEy = np.zeros([len(omega),len(incAngle),len(numpts)])
REy = np.zeros([len(omega),len(incAngle),len(numpts)])
THy = np.zeros([len(omega),len(incAngle),len(numpts)])
RHy = np.zeros([len(omega),len(incAngle),len(numpts)])

####### functions ##############################################################

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

    def input_pw(x):
          return one_function(x)*np.exp( 1.0j * kappa *x)

    eps_glass = 2.25
    eps_air = 1.
    def eps_func(x,z):
            return one_function(x)
    def eps_func_glass_struct(x,z):
            return one_function(x) #+ box_func(x,1,3)*(eps_glass-eps_air)
    def eps_func_glass(x,z):
            return one_function(x)*eps_glass
    def mu_func(x,z):
            return one_function(x)

    matmu = fmm1d.Material_struct(mu_func,mu_func,mu_func)
    mateps = fmm1d.Material_struct(eps_func,eps_func,eps_func)
    mateps_glass_struct = fmm1d.Material_struct(eps_func_glass_struct,eps_func_glass_struct,eps_func_glass_struct)
    mateps_glass = fmm1d.Material_struct(eps_func_glass,eps_func_glass,eps_func_glass)
    eps_list = [mateps,mateps_glass_struct,mateps_glass,mateps]

    starttime1 = timer()
    lstack = fmm1d.Layer_stack(numpts_tmp,kappa,omega[omega_idx],period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})
    for i in range(0,len(intervallVec)):
      lstack.add_single_layer(eps_list[i],matmu,intervallVec[i],magnus_order)

    (p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy) = lstack.get_power_data()

    stoptime = timer()

    time_needed_tmp[j] = stoptime-starttime1

    TEy_tmp[j] = pEndscatEy/p0incEy
    REy_tmp[j] = p0scatEy/p0incEy
    THy_tmp[j] = pEndscatHy/p0incHy
    RHy_tmp[j] = p0scatHy/p0incHy


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
            flog.write("Status: %(indx)i / %(num_of_calcs)i clacs run, i.e. %(percentage)2.3f %%; needed time: %(time_need).1f min, ETA in: %(time_eta).1f min \n" %{"indx": tot_calc_idx, "num_of_calcs": total_calc_num,"percentage": tot_calc_idx*100./float(total_calc_num), "time_need": (tstop_curr - tstart)/60.,"time_eta": (tstop_curr - tstart)/60. * ((float(total_calc_num)/tot_calc_idx) - 1.0) })
            flog.close()


    #save to savefile
    sio.matlab.savemat(inputargs.save_path +  save_filname,
                  {"MagnusTEy": TEy, "MagnusREy": REy, "MagnusTHy": THy, "MagnusRHy": RHy,
                    "wavelength": wavelength, "omega": omega, "incAngle": incAngle,
                    "numpts": numpts, "setup_period": setup_period, "Magnus_time_needed": time_needed})
    #delete logfile
    if os.path.isfile(inputargs.log_path + log_filname):
      os.system('rm -f ' + inputargs.log_path + log_filname)

  except KeyboardInterrupt:
    print "Caught KeyboardInterrupt, terminating workers"
    pool.terminate()
    pool.join()
