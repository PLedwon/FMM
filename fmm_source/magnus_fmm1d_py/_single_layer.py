import numpy as np
import scipy.linalg
from cmath import *
from timeit import default_timer as timer
#import matplotlib.pyplot as plt

class Single_layer():
  def __init__(self,setup_params,epsstruct,mustruct,intervall,magnus_order):
    self.numpts = setup_params['numpts']
    self.kappa = setup_params['kappa']
    self.omega = setup_params['omega']
    self.period = setup_params['period']
    self.kztolfactor = setup_params['kztolfactor']

    self.intervall = intervall
    self.thickness = self.intervall[1]-self.intervall[0]
    self.epsstruct = epsstruct
    self.mustruct = mustruct
    self.n = magnus_order

    delta = self.period/(1.0*self.numpts);
    self.xvec = np.arange(0,self.numpts)*delta;
    self.kxvec = 2.0*np.pi*(np.arange(-self.numpts/2.,(self.numpts)/2.))/self.period



    self._process_layer()

  def _process_layer(self):
      self.kzvalsHy, self.hycoeffs = self.Magnus_te()
      self.kzvalsEy, self.eycoeffs = self.Magnus_tm()


  def get_layer_modes_EzHz(self,xvec):

    EzUp = np.zeros([np.size(xvec),self.numpts],complex)
    HzUp = np.zeros([np.size(xvec),self.numpts],complex)
    EzDown = np.zeros([np.size(xvec),self.numpts],complex)
    HzDown = np.zeros([np.size(xvec),self.numpts],complex)
    N = self.numpts

    pdXdftmatW = np.zeros([np.size(xvec),self.numpts],complex)
    for i in range(0,self.numpts):
      for j in range(0,np.size(xvec)):
        pdXdftmatW[j,i] = 1.0j*(self.kxvec[i]+self.kappa)*np.exp( 1.0j*(self.kxvec[i]+self.kappa)*xvec[j])


    for i in range(0,self.numpts):

      EzUp[:,i] = 1.0j/self.omega/self.epsstruct.zz(xvec,self.thickness)*(np.dot(pdXdftmatW,self.hycoeffs[N:2*N,i]))
      HzUp[:,i] = - 1.0j/self.omega/self.mustruct.zz(xvec,self.thickness)*(np.dot(pdXdftmatW,self.eycoeffs[N:2*N,i]))

      EzDown[:,i] = 1.0j/self.omega/self.epsstruct.zz(xvec,self.thickness)*(np.dot(pdXdftmatW,self.hycoeffs[N:2*N,N+i]))
      HzDown[:,i] = - 1.0j/self.omega/self.mustruct.zz(xvec,self.thickness)*(np.dot(pdXdftmatW,self.eycoeffs[N:2*N,N+i]))

    return EzUp,EzDown,HzUp,HzDown


##################################################HIGHLY EXPERIMENTAL AREA BELOW!!!!!##########################################################################################################

  def _L_te_Interpolation(self):

    if self.n == 0:
        zm = np.array([self.intervall[0]])
    else:
        zm = self.intervall[0] + self.thickness*np.arange(0,self.n+1)/self.n

    if self.n == 0:
        B = np.matrix([1])
        L_list = [[]]
    elif self.n == 1:
        B = np.matrix([np.ones(np.size(zm)),zm]).transpose()
        L_list=[[],[]]
    elif self.n == 2:
        B = np.matrix([np.ones(np.size(zm)),zm,0.5*np.power(zm,2)]).transpose()
        L_list=[[],[],[]]
    elif self.n == 3:
        B = np.matrix([np.ones(np.size(zm)),zm,0.5*np.power(zm,2),1.0/6.0*np.power(zm,3)]).transpose()
        L_list=[[],[],[],[]]



    pdmat = 1j*np.diag(self.kxvec+self.kappa)


    for zi in range(0,len(zm)): #Erzeugung der Liste der L_te(z_i)
        epsXX = self.epsstruct.xx(self.xvec,zm[zi])
        epsZZ = self.epsstruct.zz(self.xvec,zm[zi])
        muYY = self.mustruct.yy(self.xvec,zm[zi])

        epsXXft = np.fft.fft(epsXX)/self.numpts
        invepsZZft = np.fft.fft(1.0/epsZZ)/self.numpts
        muYYft = np.fft.fft(muYY)/self.numpts

        epsXXmat = np.zeros([self.numpts,self.numpts],complex)
        muYYmat = np.zeros([self.numpts,self.numpts],complex)
        invepsZZmat = np.zeros([self.numpts,self.numpts],complex)

        for i in range(0,self.numpts):
          for j in range(0,self.numpts):
            k = np.mod(i-j,self.numpts)
            epsXXmat[i,j] = epsXXft[k]
            invepsZZmat[i,j] = invepsZZft[k]
            muYYmat[i,j] = muYYft[k]

        L01 = muYYmat + 1.0/self.omega/self.omega*np.dot(np.dot(pdmat,invepsZZmat),pdmat)
        L10 = epsXXmat
        L_list[zi]=[L01, L10]

    L_te_ord = np.zeros((len(zm),2*self.numpts,2*self.numpts),complex)
    L_ord = np.zeros((2,len(zm),self.numpts,self.numpts),complex)
    # temp field for storing L01 and L10 for different k-entries
    lpq2 = np.zeros(len(zm),complex)

    for p in range(0,self.numpts): #Berechnung der Interpolationskoeffizienten
      for q in range(0,self.numpts):
            for i in range(0,2):
                for zi in range(0,len(zm)):
                    lpq2[zi]=L_list[zi][i][p,q]
                lpq1=np.linalg.solve(B,lpq2)
                for o in range(0,len(zm)):
                    L_ord[i,o,p,q] = lpq1[o]

    L_te_ord[:,0:self.numpts,self.numpts:2*self.numpts] = L_ord[0,:,:,:]*self.omega
    L_te_ord[:,self.numpts:2*self.numpts,0:self.numpts] = L_ord[1,:,:,:]*self.omega
    self.L_te_ord=L_te_ord




    return self.L_te_ord

  def _L_tm_Interpolation(self):

    if self.n == 0:
        zm = np.array([self.intervall[0]])
    else:
        zm = self.intervall[0] + self.thickness*np.arange(0,self.n+1)/self.n

    if self.n == 0:
        B = np.matrix([1])
        L_list = [[]]
    elif self.n == 1:
        B = np.matrix([np.ones(np.size(zm)),zm]).transpose()
        L_list=[[],[]]
    elif self.n == 2:
        B = np.matrix([np.ones(np.size(zm)),zm,0.5*np.power(zm,2)]).transpose()
        L_list=[[],[],[]]
    elif self.n == 3:
        B = np.matrix([np.ones(np.size(zm)),zm,0.5*np.power(zm,2),1.0/6.0*np.power(zm,3)]).transpose()
        L_list=[[],[],[],[]]



    pdmat = 1j*np.diag(self.kxvec+self.kappa)

    for zi in range(0,len(zm)): #Erzeugung der Liste der L_te(z_i)
        muXX = self.mustruct.xx(self.xvec,zm[zi])
        muZZ = self.mustruct.zz(self.xvec,zm[zi])
        epsYY = self.epsstruct.yy(self.xvec,zm[zi])

        muXXft = np.fft.fft(muXX)/self.numpts
        invmuZZft = np.fft.fft(1.0/muZZ)/self.numpts
        epsYYft = np.fft.fft(epsYY)/self.numpts

        muXXmat = np.zeros([self.numpts,self.numpts],complex)
        invmuZZmat = np.zeros([self.numpts,self.numpts],complex)
        epsYYmat = np.zeros([self.numpts,self.numpts],complex)

        for i in range(0,self.numpts):
          for j in range(0,self.numpts):
            k = np.mod(i-j,self.numpts)
            muXXmat[i,j] = muXXft[k]
            invmuZZmat[i,j] = invmuZZft[k]
            epsYYmat[i,j] = epsYYft[k]

        L01 = -1.0*epsYYmat -1.0/self.omega/self.omega*np.dot(np.dot(pdmat,invmuZZmat),pdmat)
        L10 = -1.0*muXXmat
        L_list[zi]=[L01, L10]

    L_tm_ord = np.zeros((len(zm),2*self.numpts,2*self.numpts),complex)
    L_ord = np.zeros((2,len(zm),self.numpts,self.numpts),complex)
    lpq2 = np.zeros(len(zm),complex)

    for p in range(0,self.numpts): #Berechnung der Interpolationskoeffizienten
      for q in range(0,self.numpts):
            for i in range(0,2):
                for zi in range(0,len(zm)):
                  lpq2[zi]=L_list[zi][i][p,q]
                lpq1=np.linalg.solve(B,lpq2)
                for o in range(0,len(zm)):
                    L_ord[i,o,p,q] = lpq1[o]

    L_tm_ord[:,0:self.numpts,self.numpts:2*self.numpts] = L_ord[0,:,:,:]*self.omega
    L_tm_ord[:,self.numpts:2*self.numpts,0:self.numpts] = L_ord[1,:,:,:]*self.omega
    self.L_tm_ord=L_tm_ord



    return L_tm_ord


  def _Magnus_Expansion(self,LMat):

    def com(A,B):
        return np.dot(A,B)-np.dot(B,A)

    self.OmegaW = np.zeros((4,2*self.numpts,2*self.numpts),complex)

    if self.n ==0:
            self.OmegaW[0,:,:] = LMat[0,:,:]

    elif self.n ==1:
            self.OmegaW[0,:,:] = LMat[0,:,:]
            self.OmegaW[1,:,:] = 0.5*LMat[1,:,:]

    elif self.n == 2:
        self.OmegaW[0,:,:] = LMat[0,:,:]
        self.OmegaW[1,:,:] = 0.5*LMat[1,:,:]
        self.OmegaW[2,:,:] = 1.0/6.0 * (LMat[2,:,:]-com(LMat[0,:,:],LMat[1,:,:]))

    elif self.n == 3:
        self.OmegaW[0,:,:] = LMat[0,:,:]
        self.OmegaW[1,:,:] = 0.5*LMat[1,:,:]
        self.OmegaW[2,:,:] = 1.0/6.0 * (LMat[2,:,:]-com(LMat[0,:,:],LMat[1,:,:]))
        self.OmegaW[3,:,:] = 1.0/24.0* (LMat[3,:,:]- com(LMat[0,:,:],LMat[2,:,:])-1.0/3.0*com(LMat[0,:,:],com(LMat[0,:,:],LMat[1,:,:])))

    else:
        raise IndexError('Magnus order has to be an integer between 0 and 3')
    return self.OmegaW

  def _Magnus_kzvals_te(self,OW):

        Omega=OW[0,:,:]+OW[1,:,:]*np.power(self.thickness,1)+OW[2,:,:]*np.power(self.thickness,2)+OW[3,:,:]*np.power(self.thickness,3) #Propagator-Matrix

        deigvals, vycoeffs = scipy.linalg.eig(Omega)
        print np.amax(np.abs(deigvals)), np.amin(np.abs(deigvals))
        kzvals, indlist  = self._kzclassification(deigvals,vycoeffs)
        ycoeffs = vycoeffs[:,indlist]

        return (kzvals,ycoeffs)

  def _Magnus_kzvals_tm(self,OW):

        Omega=OW[0,:,:]+OW[1,:,:]*np.power(self.thickness,1)+OW[2,:,:]*np.power(self.thickness,2)+OW[3,:,:]*np.power(self.thickness,3) #Propagator-Matrix

        deigvals, vycoeffs = scipy.linalg.eig(Omega)
        print  np.amax(np.abs(deigvals)), np.amin(np.abs(deigvals))
        kzvals, indlist = self._kzclassification(deigvals,vycoeffs)
        ycoeffs = vycoeffs[:,indlist]

        return (kzvals,ycoeffs)

  def Magnus_te(self):

        self._L_te_Interpolation()
        LMat=self.L_te_ord
        OW=self._Magnus_Expansion(LMat)
        kzvalste, hycoeffste=self._Magnus_kzvals_te(OW)

        return(kzvalste,hycoeffste)

  def Magnus_tm(self):

        self._L_tm_Interpolation()
        LMat=self.L_tm_ord
        OW=self._Magnus_Expansion(LMat)
        kzvalstm, eycoeffstm=self._Magnus_kzvals_tm(OW)

        return(kzvalstm,eycoeffstm)

  def _kzclassification(self,kz_old,vcoeffs):

      kztol = np.max(np.abs(kz_old))*self.kztolfactor

      kzu= []
      kzd=[]

      indu = []
      indd = []
      indzero = []
      is_left = []

      for i in range(0,np.size(kz_old)):
          if kz_old[i].real >= kztol:
              if kz_old[i].imag < -kztol:
                  kzd.append(kz_old[i])
                  indd.append(i)
              else:
                  kzu.append(kz_old[i])
                  indu.append(i)
          elif kz_old[i].real <=  - kztol:
              if kz_old[i].imag <= kztol:
                  kzd.append(kz_old[i])
                  indd.append(i)
              else:
                  kzu.append(kz_old[i])
                  indu.append(i)
          else:
              if kz_old[i].imag > kztol:
                  kzu.append(kz_old[i])
                  indu.append(i)
              elif kz_old[i].imag < - kztol:
                  kzd.append(kz_old[i])
                  indd.append(i)
              else:
                  indzero.append(i)
                  tmp_idx = np.argmax(np.abs(vcoeffs[:,i]))
                  if ((tmp_idx < self.numpts/2) or (tmp_idx > self.numpts and tmp_idx < 3*self.numpts/2)):
                    is_left.append(True)
                  else:
                    is_left.append(False)

#                  if kz_old[i].real < 0.0:
#                      kzd.append(kz_old[i] - (1.0 + np.finfo(np.complex_).eps)*kztol)
#                      indd.append(i)
#                  else:
#                      kzu.append(kz_old[i] + (1.0 + np.finfo(np.complex_).eps)*kztol)
#                      indu.append(i)
      leftcounter = 0
      rightcounter = 0
      for i in range(0,len(indzero)):
        if is_left[i]:
            if leftcounter%2==0:
                indu.append(indzero[i])
                kzu.append(kz_old[indzero[i]] + (1.0 + np.finfo(np.complex_).eps)*kztol)
            else:
                kzd.append(kz_old[indzero[i]] - (1.0 + np.finfo(np.complex_).eps)*kztol)
                indd.append(indzero[i])
            leftcounter += 1
        else:
            if rightcounter%2==0:
                indu.append(indzero[i])
                kzu.append(kz_old[indzero[i]] + (1.0 + np.finfo(np.complex_).eps)*kztol)
            else:
                kzd.append(kz_old[indzero[i]] - (1.0 + np.finfo(np.complex_).eps)*kztol)
                indd.append(indzero[i])
            rightcounter += 1

      self.Nup=len(kzu)
      self.Ndown=len(kzd)

      kzlist = []
      kzlist.extend(kzu)
      kzlist.extend(kzd)

      indlist = []
      indlist.extend(indu)
      indlist.extend(indd)

      kz_new = np.array(kzlist)
      if len(indzero)%2==1:
          raise IndexError('number of zero-EV is uneven, up-down-classification is unclear')
      if self.Nup != self.Ndown:
          raise IndexError('Number of up-modes != number of down-modes')

      return kz_new, indlist
