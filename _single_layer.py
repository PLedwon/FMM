import numpy as np
from cmath import *

from timeit import default_timer as timer

class Single_layer():
  def __init__(self,setup_params,thickness,epsstruct,mustruct):
    self.numpts = setup_params['numpts']
    self.kappa = setup_params['kappa']
    self.omega = setup_params['omega']
    self.period = setup_params['period']
    self.kztolfactor = setup_params['kztolfactor']

    self.thickness = thickness
    self.epsstruct = epsstruct
    self.mustruct = mustruct

    delta = self.period/self.numpts;
    self.xvec = np.arange(0,self.numpts)*delta;
    self.kxvec = 2.0*np.pi*(np.arange(-self.numpts/2.,(self.numpts)/2.))/self.period

    self._process_layer()

  def _process_layer(self):
    self.kzvalsHy, self.hycoeffs = self._tesinglelayer()
    self.kzvalsEy, self.eycoeffs = self._tmsinglelayer()

  def _kzfilter(self,kz_old,kztol):
    kz_new = np.zeros(np.size(kz_old),complex);
    for i in range(0,np.size(kz_old)):
      if kz_old[i].real >= kztol:
        if kz_old[i].imag < -kztol:
          kz_new[i] = - kz_old[i]
        else:
          kz_new[i] = kz_old[i]
      elif kz_old[i].real <=  - kztol:
        if kz_old[i].imag <= kztol:
          kz_new[i] = - kz_old[i]
        else:
          kz_new[i] = kz_old[i]
      else:
        if kz_old[i].imag > kztol:
          kz_new[i] = kz_old[i]
        elif kz_old[i].imag < - kztol:
          kz_new[i] = -kz_old[i]
        else:
          kz_new[i] = -kz_old[i] + 1.0j *(1.0 + np.finfo(np.complex_).eps)*kztol;
    return kz_new

  def _tesinglelayer(self):
    epsXX = self.epsstruct.xx(self.xvec)
    epsZZ = self.epsstruct.zz(self.xvec)
    muYY = self.mustruct.yy(self.xvec)

    pdmat = 1j*np.diag(self.kxvec+self.kappa)
    epsXXft = np.fft.fft(epsXX)/self.numpts
    muYYft = np.fft.fft(muYY)/self.numpts
    invepsZZft = np.fft.fft(1.0/epsZZ)/self.numpts

    epsXXmat = np.zeros([self.numpts,self.numpts],complex)
    muYYmat = np.zeros([self.numpts,self.numpts],complex)
    invepsZZmat = np.zeros([self.numpts,self.numpts],complex)
    for i in range(0,self.numpts):
      for j in range(0,self.numpts):
        k = np.mod(i-j,self.numpts)
        epsXXmat[i,j] = epsXXft[k]
        muYYmat[i,j] = muYYft[k]
        invepsZZmat[i,j] = invepsZZft[k]

    deigvals, vhycoeffs = np.linalg.eig(np.dot(epsXXmat,(muYYmat*self.omega*self.omega + np.dot(np.dot(pdmat,invepsZZmat),pdmat))))

    deigvals = np.sqrt(deigvals)
    didx = np.argsort(deigvals)

    kzvals = deigvals[didx]
    kztol = np.max(np.abs(kzvals))*self.kztolfactor
    kzvals = self._kzfilter(kzvals,kztol)

    hycoeffs = vhycoeffs[:,didx]

    return (kzvals,hycoeffs)

  def _tmsinglelayer(self):
    epsYY = self.epsstruct.yy(self.xvec)
    muXX = self.mustruct.xx(self.xvec)
    muZZ = self.mustruct.zz(self.xvec)

    pdmat = 1j*np.diag(self.kxvec+self.kappa)
    epsYYft = np.fft.fft(epsYY)/self.numpts
    muXXft = np.fft.fft(muXX)/self.numpts
    invmuZZft = np.fft.fft(1.0/muZZ)/self.numpts

    epsYYmat = np.zeros([self.numpts,self.numpts],complex)
    muXXmat = np.zeros([self.numpts,self.numpts],complex)
    invmuZZmat = np.zeros([self.numpts,self.numpts],complex)
    for i in range(0,self.numpts):
      for j in range(0,self.numpts):
        k = np.mod(i-j,self.numpts)
        epsYYmat[i,j] = epsYYft[k]
        muXXmat[i,j] = muXXft[k]
        invmuZZmat[i,j] = invmuZZft[k]

    deigvals, veycoeffs = np.linalg.eig(np.dot(muXXmat,(epsYYmat*self.omega*self.omega + np.dot(np.dot(pdmat,invmuZZmat),pdmat))))

    deigvals = np.sqrt(deigvals)
    didx = np.argsort(deigvals)

    kzvals = deigvals[didx]
    kztol = np.max(np.abs(kzvals))*self.kztolfactor
    kzvals = self._kzfilter(kzvals,kztol)

    eycoeffs = veycoeffs[:,didx]

    return (kzvals,eycoeffs)

  def get_layer_modes(self,x):

    Ex = np.zeros([np.size(x),self.numpts],complex)
    Hx = np.zeros([np.size(x),self.numpts],complex)
    Ez = np.zeros([np.size(x),self.numpts],complex)
    Hz = np.zeros([np.size(x),self.numpts],complex)

    dftmatW = np.zeros([np.size(x),self.numpts],complex)
    pdXdftmatW = np.zeros([np.size(x),self.numpts],complex)
    for i in range(0,self.numpts):
      for j in range(0,np.size(x)):
        dftmatW[j,i] = np.exp( 1.0j*(self.kxvec[i]+self.kappa)*x[j])
        pdXdftmatW[j,i] = 1.0j*(self.kxvec[i]+self.kappa)*dftmatW[j,i]

    Ey = np.dot(dftmatW,self.eycoeffs)
    Hy = np.dot(dftmatW,self.hycoeffs)

    for i in range(0,self.numpts):
      Ex[:,i] = self.kzvalsHy[i]/self.omega/self.epsstruct.xx(x)*(np.dot(dftmatW,self.hycoeffs[:,i]))
      Hx[:,i] = -self.kzvalsEy[i]/self.omega/self.mustruct.xx(x)*(np.dot(dftmatW,self.eycoeffs[:,i]))

      Ez[:,i] = 1.0j/self.omega/self.epsstruct.zz(x)*(np.dot(pdXdftmatW,self.hycoeffs[:,i]))
      Hz[:,i] = - 1.0j/self.omega/self.mustruct.zz(x)*(np.dot(pdXdftmatW,self.eycoeffs[:,i]))

    return Ex,Ey,Ez,Hx,Hy,Hz
