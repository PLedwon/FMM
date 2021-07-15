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
parser.add_argument('-c','--calcnum'    ,type=int,help='run calculation number N',action='store',metavar='N'    ,dest='clac_idx',required=True)
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



# =================
# helping functions
# =================

def zero_func(x):
  return np.zeros(x.shape)

def one_function(x):
  return np.ones(x.shape)

def box_func(x,a,b):
  return 0.5*(np.sign(x-a) - np.sign(x-b))

# ==================
# load material data
# ==================

olmon_gold = np.loadtxt('/home/kiel/runs/beamsplitter_setup_simple_AuTi_grating_non-dispersive_diamond/Olmon-ev_gold.txt',skiprows = 9)
ordal_titanium = np.loadtxt('/home/kiel/runs/beamsplitter_setup_simple_AuTi_grating_non-dispersive_diamond/Ordal_titanium.txt',skiprows = 9)


def eps_au(lmbd):
  return np.interp(lmbd,olmon_gold[:,0],pow(olmon_gold[:,1]+ 1.0j * olmon_gold[:,2],2).real) + 1.0j*np.interp(lmbd,olmon_gold[:,0],pow(olmon_gold[:,1]+ 1.0j * olmon_gold[:,2],2).imag)

def eps_ti(lmbd):
  return np.interp(lmbd,ordal_titanium[:,0],pow(ordal_titanium[:,1]+ 1.0j * ordal_titanium[:,2],2).real) + 1.0j*np.interp(lmbd,ordal_titanium[:,0],pow(ordal_titanium[:,1]+ 1.0j * ordal_titanium[:,2],2).imag)

# =========
# init data
# =========

wavelength = np.arange(2.5,25.25,0.25)
#numpts = np.array([16,22,32,46,64,90,128,182,256,362,512])#,724,1024])
numpts = np.array([16,32,64,128,256,362,512])#,724,1024])
gold_thickness = np.array([0.104])
titanium_thickness = np.array([0.015])
setup_period = np.array([(1.93+1.05)])
au_ti_widths_scanpar = np.arange(0.0,1.93+1.05,0.025)
# just add a single width combination to the beginning
au_ti_widths = np.append([[1.93,1.93],[1.93+1.05,1.93+1.05]],np.transpose(np.array([au_ti_widths_scanpar,au_ti_widths_scanpar])),axis=0)
incAngle = np.array([67.2])

#test data
#wavelength = np.arange(10.0,25.1,10.0)
#numpts = np.array([16,22,32,46,64,90,128,182,256,362,512,724,1024])
#gold_thickness = np.array([0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25])
#setup_period = np.array([0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0])
#coverage_factor = np.arange(0.0,1.025,0.25)
#incAngle = np.array([0.0,30.,45.0,60.,62.7,67.2])

num_of_runs = len(gold_thickness)*len(setup_period)
if inputargs.clac_idx >= num_of_runs:
  raise IndexError('There is no' + str(inputargs.clac_idx) + '-th calculation to be run out of' + str(num_of_runs) + ' calculations!')

current_gold_thickness_idx = int(inputargs.clac_idx)%len(gold_thickness)
current_setup_period_idx = int(inputargs.clac_idx)/len(gold_thickness)

current_gold_thickness = gold_thickness[current_gold_thickness_idx]
current_setup_period   = setup_period[current_setup_period_idx]

current_titanium_thickness = titanium_thickness[0]


# =========
# filenames
# =========

log_filname = 'beamsplitter_scan_data_period__%(periodpar)10.6E__thickness__%(thickpar)10.6E__.log' %{"periodpar": current_setup_period, "thickpar": current_gold_thickness}
save_filname = 'beamsplitter_scan_data_period__%(periodpar)10.6E__thickness__%(thickpar)10.6E__.mat' %{"periodpar": current_setup_period, "thickpar": current_gold_thickness}

if not os.path.exists(inputargs.log_path):
  os.makedirs(inputargs.log_path)

# =========

omega = 2*pi/wavelength

eps_au_vals = np.zeros([len(omega)],complex)
eps_ti_vals = np.zeros([len(omega)],complex)

for i in range(0,len(omega)):
  eps_au_vals[i] = eps_au(wavelength[i])
  eps_ti_vals[i] = eps_ti(wavelength[i])

time_needed = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])

TEy = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
REy = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
THy = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
RHy = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])

TEy_ZO = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
REy_ZO = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
THy_ZO = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])
RHy_ZO = np.zeros([len(omega),len(au_ti_widths),len(incAngle),len(numpts)])

def run_func(omega_idx,auti_idx,incAng_idx):
  TEy_tmp = np.zeros([len(numpts)])
  REy_tmp = np.zeros([len(numpts)])
  THy_tmp = np.zeros([len(numpts)])
  RHy_tmp = np.zeros([len(numpts)])

  TEy_ZO_tmp = np.zeros([len(numpts)])
  REy_ZO_tmp = np.zeros([len(numpts)])
  THy_ZO_tmp = np.zeros([len(numpts)])
  RHy_ZO_tmp = np.zeros([len(numpts)])

  time_needed_tmp = np.zeros([len(numpts)])

  for j in range(0,len(numpts)):
    numpts_tmp = numpts[j]
    kappa = omega[omega_idx]*sin(pi*incAngle[incAng_idx]/180)#*sqrt(2.25)*(1.0+1e-15j)
    period = current_setup_period
    kztolfac = 1e-10

    layer_thickness = np.array([500, current_gold_thickness, current_titanium_thickness, 500])
    eps1 = np.array([1.0,5.65,eps_au_vals[omega_idx],eps_ti_vals[omega_idx]])
    def eps1_func(x):
      return eps1[0]*one_function(x)
    def eps2_func(x):
      return eps1[0]*one_function(x) + (eps1[2]-eps1[0])*box_func(np.mod(x,period),0.0*period,au_ti_widths[auti_idx,0])
    def eps3_func(x):
      return eps1[0]*one_function(x) + (eps1[3]-eps1[0])*box_func(np.mod(x,period),0.0*period,au_ti_widths[auti_idx,1])
    def eps4_func(x):
      return eps1[1]*one_function(x)
    def mu_func(x):
      return one_function(x)

    def input_pw(x):
      return one_function(x)*np.exp( 1.0j * kappa *x)

    matmu1 = fmm1d.Material_struct(mu_func,mu_func,mu_func)
    mateps1 = fmm1d.Material_struct(eps1_func,eps1_func,eps1_func)
    mateps2 = fmm1d.Material_struct(eps2_func,eps2_func,eps2_func)
    mateps3 = fmm1d.Material_struct(eps3_func,eps3_func,eps3_func)
    mateps4 = fmm1d.Material_struct(eps4_func,eps4_func,eps4_func)

    eps_list = [mateps1,mateps2,mateps3,mateps4]
    mu_list = [matmu1,matmu1,matmu1,matmu1]

    starttime1 = timer()
    lstack = fmm1d.Layer_stack(numpts_tmp,kappa,omega[omega_idx],period,kztolfac,{"uinfuncEy": input_pw,"uinfuncHy": input_pw,"dinfuncEy": zero_func,"dinfuncHy": zero_func})
    for i in range(0,len(layer_thickness)):
      lstack.add_single_layer(layer_thickness[i],eps_list[i],mu_list[i])

    (p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy) = lstack.get_power_data()
    difforder_data = lstack.get_power_data_by_orders()
    stoptime = timer()

    time_needed_tmp[j] = stoptime-starttime1

    TEy_tmp[j] = pEndscatEy/p0incEy
    REy_tmp[j] = - p0scatEy/p0incEy
    THy_tmp[j] = pEndscatHy/p0incHy
    RHy_tmp[j] = - p0scatHy/p0incHy

    do_idx = difforder_data.order_number.index(0)

    TEy_ZO_tmp[j] = difforder_data.pEndscatEy[do_idx]/p0incEy
    REy_ZO_tmp[j] = - difforder_data.p0scatEy[do_idx]/p0incEy
    THy_ZO_tmp[j] = difforder_data.pEndscatHy[do_idx]/p0incHy
    RHy_ZO_tmp[j] = - difforder_data.p0scatHy[do_idx]/p0incHy

    #print str(numpts[j]) + ' ' + str(time_needed_tmp[j])

  return (TEy_tmp,REy_tmp,THy_tmp,RHy_tmp,TEy_ZO_tmp,REy_ZO_tmp,THy_ZO_tmp,RHy_ZO_tmp,time_needed_tmp)

if __name__== '__main__':
  pool = multiprocessing.Pool(inputargs.num_proc,init_worker)
  total_calc_num = len(incAngle)*len(au_ti_widths)*len(omega)
  try:
    tstart = timer()
    time_idx = 0
    tot_calc_idx = 0
    for incAng_idx in range(0,len(incAngle)):
      res = []
      for auti_idx in range(0,len(au_ti_widths)):
        for omega_idx in range(0,len(omega)):
          #print run_func(omega_idx,auti_idx,incAng_idx)
          res.append(pool.apply_async(run_func,(omega_idx,auti_idx,incAng_idx,)))
      res_idx = 0
      for auti_idx in range(0,len(au_ti_widths)):
        for omega_idx in range(0,len(omega)):
          (TEy_tmp,REy_tmp,THy_tmp,RHy_tmp,TEy_ZO_tmp,REy_ZO_tmp,THy_ZO_tmp,RHy_ZO_tmp,time_needed_tmp) = res[res_idx].get()
          res_idx = res_idx + 1
          tot_calc_idx = tot_calc_idx + 1

          TEy[omega_idx,auti_idx,incAng_idx,:] = TEy_tmp
          REy[omega_idx,auti_idx,incAng_idx,:] = REy_tmp
          THy[omega_idx,auti_idx,incAng_idx,:] = THy_tmp
          RHy[omega_idx,auti_idx,incAng_idx,:] = RHy_tmp

          TEy_ZO[omega_idx,auti_idx,incAng_idx,:] = TEy_ZO_tmp
          REy_ZO[omega_idx,auti_idx,incAng_idx,:] = REy_ZO_tmp
          THy_ZO[omega_idx,auti_idx,incAng_idx,:] = THy_ZO_tmp
          RHy_ZO[omega_idx,auti_idx,incAng_idx,:] = RHy_ZO_tmp

          time_needed[omega_idx,auti_idx,incAng_idx,:] = time_needed_tmp
          tstop_curr = timer()
          if int((tstop_curr - tstart ) / inputargs.log_time) > time_idx:
            time_idx = int((tstop_curr - tstart ) / inputargs.log_time)
            #save to logfile
            flog = open(inputargs.log_path + log_filname,'a')
            flog.write("Status: %(indx)i / %(num_of_calcs)i clacs run, i.e. %(percentage)2.3f %%; needed time: %(time_need).1f min, ETA in: %(time_eta).1f min \n" %{"indx": tot_calc_idx, "num_of_calcs": total_calc_num,"percentage": tot_calc_idx*100./float(total_calc_num), "time_need": (tstop_curr - tstart)/60.,"time_eta": (tstop_curr - tstart)/60. * ((float(total_calc_num)/tot_calc_idx) - 1.0) })
            flog.close()


    #save to savefile
    sio.matlab.savemat(inputargs.save_path +  save_filname,
                  {"TEy": TEy, "REy": REy, "THy": THy, "RHy": RHy,
                    "TEy_ZO": TEy_ZO, "REy_ZO": REy_ZO, "THy_ZO": THy_ZO, "RHy_ZO": RHy_ZO,
                    "wavelength": wavelength, "omega": omega, "incAngle": incAngle,
                    "numpts": numpts, "gold_thickness": gold_thickness, "setup_period": setup_period,
                    "au_ti_widths": au_ti_widths, "current_gold_thickness_idx": current_gold_thickness_idx,
                    "titanium_thickness": current_titanium_thickness,
                    "current_setup_period_idx": current_setup_period_idx, "time_needed": time_needed})
    #delete logfile
    if os.path.isfile(inputargs.log_path + log_filname):
      os.system('rm -f ' + inputargs.log_path + log_filname)

  except KeyboardInterrupt:
    print "Caught KeyboardInterrupt, terminating workers"
    pool.terminate()
    pool.join()
