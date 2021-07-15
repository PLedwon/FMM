import numpy as np
from cmath import *
from _single_layer import *

class Layer_stack():
  def __init__(self,numpts,kappa,omega,period,kztolfactor,input_amplitudes):
    self.num_of_layers = 0
    self.numpts = numpts
    self.kappa = kappa
    self.omega = omega
    self.period = period
    self.kztolfactor = kztolfactor

    self.uinfuncEy = input_amplitudes['uinfuncEy']
    self.uinfuncHy = input_amplitudes['uinfuncHy']
    self.dinfuncEy = input_amplitudes['dinfuncEy']
    self.dinfuncHy = input_amplitudes['dinfuncHy']

    self.layer_list = []

    self._smat_calculated = False
    self._all_ampls_calculated = False

    self._SmatEy = None
    self._SmatHy = None

    self._uamplsEy = []
    self._uamplsHy = []
    self._damplsEy = []
    self._damplsHy = []

  def add_single_layer(self,thickness,epsstruct,mustruct):

    self.layer_list.append(Single_layer({'numpts': self.numpts,'kappa': self.kappa,'omega': self.omega,'period': self.period,'kztolfactor': self.kztolfactor},thickness,epsstruct,mustruct))

    self.num_of_layers += 1

    self._set_calculation_status_false()


  def change_input_amplitudes(self,input_amplitudes):
    self.uinfuncEy = input_amplitudes['uinfuncEy']
    self.uinfuncHy = input_amplitudes['uinfuncHy']
    self.dinfuncEy = input_amplitudes['dinfuncEy']
    self.dinfuncHy = input_amplitudes['dinfuncHy']

    self._all_ampls_calculated = False
    self._reset_amplitudes()


  def _calculate_smatrix(self):
    if not(self._smat_calculated):
      self._SmatEy = np.eye(2*self.numpts)
      self._SmatHy = np.eye(2*self.numpts)

      if self.num_of_layers < 1:
        raise IndexError('number of layers must be at least one!')

      self._SmatEy = self._s_star_prop_prod(self._SmatEy, self.layer_list[0].kzvalsEy, self.layer_list[0].thickness)
      self._SmatHy = self._s_star_prop_prod(self._SmatHy, self.layer_list[0].kzvalsHy, self.layer_list[0].thickness)
      for i in range(1,self.num_of_layers):
        (transferI_Ey,transferI_Hy) = self._calculate_interface_matrix(self.layer_list[i-1],self.layer_list[i])
        self._SmatEy = self._s_star_t_prod(self._SmatEy,transferI_Ey)
        self._SmatHy = self._s_star_t_prod(self._SmatHy,transferI_Hy)
        self._SmatEy = self._s_star_prop_prod(self._SmatEy, self.layer_list[i].kzvalsEy, self.layer_list[i].thickness)
        self._SmatHy = self._s_star_prop_prod(self._SmatHy, self.layer_list[i].kzvalsHy, self.layer_list[i].thickness)

      self._smat_calculated = True

  def _calculate_all_ampls(self):
    self._uamplsEy[0], self._uamplsHy[0] = self._calculate_up_input_ampls(self.uinfuncEy,self.uinfuncHy,0.0)
    self._damplsEy[self.num_of_layers-1], self._damplsHy[self.num_of_layers-1] = self._calculate_down_input_ampls(self.dinfuncEy,self.dinfuncHy,0.0)

    upper_uampls0_Ey = np.exp( 1.0j * self.layer_list[0].kzvalsEy*self.layer_list[0].thickness)*self._uamplsEy[0]
    upper_uampls0_Hy = np.exp( 1.0j * self.layer_list[0].kzvalsHy*self.layer_list[0].thickness)*self._uamplsHy[0]

    N = self.numpts

    for j in range(self.num_of_layers-1,0,-1):
      if j < self.num_of_layers-1:
        (invtransferI_Ey, invtransferI_Hy) = self._calculate_interface_matrix(self.layer_list[j+1], self.layer_list[j])
        self._damplsEy[j] = np.dot(invtransferI_Ey[0:N,N:2*N],self._uamplsEy[j+1]) + np.dot(invtransferI_Ey[N:2*N,N:2*N],np.exp( 1.0j * self.layer_list[j+1].kzvalsEy*self.layer_list[j+1].thickness)*self._damplsEy[j+1])
        self._damplsHy[j] = np.dot(invtransferI_Hy[0:N,N:2*N],self._uamplsHy[j+1]) + np.dot(invtransferI_Hy[N:2*N,N:2*N],np.exp( 1.0j * self.layer_list[j+1].kzvalsHy*self.layer_list[j+1].thickness)*self._damplsHy[j+1])

      # calculate S-Matrix for specific setup
      StotEy = np.eye(2*N)
      StotHy = np.eye(2*N)
      for i in range(0,j-1):
        (transferI_Ey,transferI_Hy) = self._calculate_interface_matrix(self.layer_list[i],self.layer_list[i+1])
        StotEy = self._s_star_t_prod(StotEy,transferI_Ey)
        StotHy = self._s_star_t_prod(StotHy,transferI_Hy)
        StotEy = self._s_star_prop_prod(StotEy,self.layer_list[i+1].kzvalsEy,self.layer_list[i+1].thickness)
        StotHy = self._s_star_prop_prod(StotHy,self.layer_list[i+1].kzvalsHy,self.layer_list[i+1].thickness)
      (transferI_Ey,transferI_Hy) = self._calculate_interface_matrix(self.layer_list[j-1],self.layer_list[j])
      StotEy = self._s_star_t_prod(StotEy,transferI_Ey)
      StotHy = self._s_star_t_prod(StotHy,transferI_Hy)

      # calculate missing up and down ampls for the missing layer

      self._uamplsEy[j] = np.dot(StotEy[0:N,0:N],upper_uampls0_Ey) + np.dot(StotEy[0:N,N:2*N],np.exp( 1.0j * self.layer_list[j].kzvalsEy*self.layer_list[j].thickness)*self._damplsEy[j])
      self._uamplsHy[j] = np.dot(StotHy[0:N,0:N],upper_uampls0_Hy) + np.dot(StotHy[0:N,N:2*N],np.exp( 1.0j * self.layer_list[j].kzvalsHy*self.layer_list[j].thickness)*self._damplsHy[j])

      self._damplsEy[0] = np.dot(StotEy[N:2*N,0:N],upper_uampls0_Ey) + np.dot(StotEy[N:2*N,N:2*N],np.exp( 1.0j * self.layer_list[j].kzvalsEy*self.layer_list[j].thickness)*self._damplsEy[j])
      self._damplsHy[0] = np.dot(StotHy[N:2*N,0:N],upper_uampls0_Hy) + np.dot(StotHy[N:2*N,N:2*N],np.exp( 1.0j * self.layer_list[j].kzvalsHy*self.layer_list[j].thickness)*self._damplsHy[j])

    self._all_ampls_calculated = True

  def get_Ey_Smat(self):
    if not(self._smat_calculated):
      self.Smat = self._calculate_smatrix()
    return self._SmatEy

  def get_Hy_Smat(self):
    if not(self._smat_calculated):
      self.Smat = self._calculate_smatrix()
    return self._SmatHy

  def get_Smats(self):
    if not(self._smat_calculated):
      self.Smat = self._calculate_smatrix()
    return (self._SmatEy, self._SmatHy)

  def get_modes_at(self,xvec,zvec):
    if not(self._all_ampls_calculated):
      self._calculate_all_ampls()

    #init fields
    ex = np.zeros([np.size(xvec),np.size(zvec)],complex)
    ey = np.zeros([np.size(xvec),np.size(zvec)],complex)
    ez = np.zeros([np.size(xvec),np.size(zvec)],complex)
    hx = np.zeros([np.size(xvec),np.size(zvec)],complex)
    hy = np.zeros([np.size(xvec),np.size(zvec)],complex)
    hz = np.zeros([np.size(xvec),np.size(zvec)],complex)

    # lower halfplane

    (extmp,eytmp,eztmp,hxtmp,hytmp,hztmp) = self.layer_list[0].get_layer_modes(xvec)

    lowerzbdry = 0.0
    upperzbdry = self.layer_list[0].thickness

    for k in range(0,np.size(zvec)):
      if zvec[k] <= upperzbdry:
        #ey related fields  def get_power_data(self):

        current_Ey_uampl = np.exp( 1.0j * self.layer_list[0].kzvalsEy*(zvec[k]-lowerzbdry))*self._uamplsEy[0]
        current_Ey_dampl = np.exp( 1.0j * self.layer_list[0].kzvalsEy*(upperzbdry-zvec[k]))*self._damplsEy[0]

        hx[:,k] = np.dot(hxtmp,current_Ey_uampl) - np.dot(hxtmp,current_Ey_dampl)
        ey[:,k] = np.dot(eytmp,current_Ey_uampl) + np.dot(eytmp,current_Ey_dampl)
        hz[:,k] = np.dot(hztmp,current_Ey_uampl) + np.dot(hztmp,current_Ey_dampl)
        #hy related fields
        current_Hy_uampl = np.exp( 1.0j * self.layer_list[0].kzvalsHy*(zvec[k]-lowerzbdry))*self._uamplsHy[0]
        current_Hy_dampl = np.exp( 1.0j * self.layer_list[0].kzvalsHy*(upperzbdry-zvec[k]))*self._damplsHy[0]

        ex[:,k] = np.dot(extmp,current_Hy_uampl) - np.dot(extmp,current_Hy_dampl)
        hy[:,k] = np.dot(hytmp,current_Hy_uampl) + np.dot(hytmp,current_Hy_dampl)
        ez[:,k] = np.dot(eztmp,current_Hy_uampl) + np.dot(eztmp,current_Hy_dampl)

    # mid planes
    for j in range(1,self.num_of_layers-1):
      (extmp,eytmp,eztmp,hxtmp,hytmp,hztmp) = self.layer_list[j].get_layer_modes(xvec)

      lowerzbdry = lowerzbdry + self.layer_list[j-1].thickness
      upperzbdry = upperzbdry + self.layer_list[j].thickness

      for k in range(0,np.size(zvec)):
        if zvec[k] <= upperzbdry and zvec[k] > lowerzbdry:
          #ey related fields
          current_Ey_uampl = np.exp( 1.0j * self.layer_list[j].kzvalsEy*(zvec[k]-lowerzbdry))*self._uamplsEy[j]
          current_Ey_dampl = np.exp( 1.0j * self.layer_list[j].kzvalsEy*(upperzbdry-zvec[k]))*self._damplsEy[j]

          hx[:,k] = np.dot(hxtmp,current_Ey_uampl) - np.dot(hxtmp,current_Ey_dampl)
          ey[:,k] = np.dot(eytmp,current_Ey_uampl) + np.dot(eytmp,current_Ey_dampl)
          hz[:,k] = np.dot(hztmp,current_Ey_uampl) + np.dot(hztmp,current_Ey_dampl)
          #hy related fields
          current_Hy_uampl = np.exp( 1.0j * self.layer_list[j].kzvalsHy*(zvec[k]-lowerzbdry))*self._uamplsHy[j]
          current_Hy_dampl = np.exp( 1.0j * self.layer_list[j].kzvalsHy*(upperzbdry-zvec[k]))*self._damplsHy[j]

          ex[:,k] = np.dot(extmp,current_Hy_uampl) - np.dot(extmp,current_Hy_dampl)
          hy[:,k] = np.dot(hytmp,current_Hy_uampl) + np.dot(hytmp,current_Hy_dampl)
          ez[:,k] = np.dot(eztmp,current_Hy_uampl) + np.dot(eztmp,current_Hy_dampl)

    #upper half plane
    (extmp,eytmp,eztmp,hxtmp,hytmp,hztmp) = self.layer_list[self.num_of_layers-1].get_layer_modes(xvec)

    lowerzbdry = lowerzbdry + self.layer_list[self.num_of_layers-2].thickness
    upperzbdry = upperzbdry + self.layer_list[self.num_of_layers-1].thickness

    for k in range(0,np.size(zvec)):
      if zvec[k] > lowerzbdry:
        #ey related fields
        current_Ey_uampl = np.exp( 1.0j * self.layer_list[self.num_of_layers-1].kzvalsEy*(zvec[k]-lowerzbdry))*self._uamplsEy[self.num_of_layers-1]
        current_Ey_dampl = np.exp( 1.0j * self.layer_list[self.num_of_layers-1].kzvalsEy*(upperzbdry-zvec[k]))*self._damplsEy[self.num_of_layers-1]

        hx[:,k] = np.dot(hxtmp,current_Ey_uampl) - np.dot(hxtmp,current_Ey_dampl)
        ey[:,k] = np.dot(eytmp,current_Ey_uampl) + np.dot(eytmp,current_Ey_dampl)
        hz[:,k] = np.dot(hztmp,current_Ey_uampl) + np.dot(hztmp,current_Ey_dampl)
        #hy related fields
        current_Hy_uampl = np.exp( 1.0j * self.layer_list[self.num_of_layers-1].kzvalsHy*(zvec[k]-lowerzbdry))*self._uamplsHy[self.num_of_layers-1]
        current_Hy_dampl = np.exp( 1.0j * self.layer_list[self.num_of_layers-1].kzvalsHy*(upperzbdry-zvec[k]))*self._damplsHy[self.num_of_layers-1]

        ex[:,k] = np.dot(extmp,current_Hy_uampl) - np.dot(extmp,current_Hy_dampl)
        hy[:,k] = np.dot(hytmp,current_Hy_uampl) + np.dot(hytmp,current_Hy_dampl)
        ez[:,k] = np.dot(eztmp,current_Hy_uampl) + np.dot(eztmp,current_Hy_dampl)

    return (ex,ey,ez,hx,hy,hz)


  def get_modes_at_rotated_frame(self,xnvec,znvec,rot_angle):
    Xnmesh,Znmesh = np.meshgrid(xnvec,znvec);
    Xmesh = Xnmesh*cos(pi/180.0*rot_angle) - Znmesh*sin(pi/180.0*rot_angle)
    Zmesh = Xnmesh*sin(pi/180.0*rot_angle) + Znmesh*cos(pi/180.0*rot_angle)

    exn = np.zeros([np.size(xnvec),np.size(znvec)],complex)
    eyn = np.zeros([np.size(xnvec),np.size(znvec)],complex)
    ezn = np.zeros([np.size(xnvec),np.size(znvec)],complex)
    hxn = np.zeros([np.size(xnvec),np.size(znvec)],complex)
    hyn = np.zeros([np.size(xnvec),np.size(znvec)],complex)
    hzn = np.zeros([np.size(xnvec),np.size(znvec)],complex)

    for i in range(0,np.size(Xmesh,0)):
      for j in range(0,np.size(Xmesh,1)):
        (ex,ey,ez,hx,hy,hz)=self.get_modes_at(np.array(Xmesh[i][j]),np.array(Zmesh[i][j]))
        exn[j][i] = ex*cos(pi/180.0*rot_angle) + ez*sin(pi/180.0*rot_angle)
        eyn[j][i] = ey
        ezn[j][i] = -ex*sin(pi/180.0*rot_angle) + ez*cos(pi/180.0*rot_angle)

        hxn[j][i] = hx*cos(pi/180.0*rot_angle) + hz*sin(pi/180.0*rot_angle)
        hyn[j][i] = hy
        hzn[j][i] = -hx*sin(pi/180.0*rot_angle) + hz*cos(pi/180.0*rot_angle)

    return (exn,eyn,ezn,hxn,hyn,hzn)


  def get_power_data(self):
    if not(self._smat_calculated):
      self.Smat = self._calculate_smatrix()

    u1Ey, u1Hy = self._calculate_up_input_ampls(self.uinfuncEy,self.uinfuncHy,0.0)
    dNEy, dNHy = self._calculate_down_input_ampls(self.dinfuncEy,self.dinfuncHy,0.0)

    N = self.numpts

    #calculation of the amplitude vectors

    uNEy = np.dot(self._SmatEy[0:N,0:N],u1Ey) + np.dot(self._SmatEy[0:N,N:2*N],dNEy)
    d1Ey = np.dot(self._SmatEy[N:2*N,0:N],u1Ey) + np.dot(self._SmatEy[N:2*N,N:2*N],dNEy)

    uNHy = np.dot(self._SmatHy[0:N,0:N],u1Hy) + np.dot(self._SmatHy[0:N,N:2*N],dNHy)
    d1Hy = np.dot(self._SmatHy[N:2*N,0:N],u1Hy) + np.dot(self._SmatHy[N:2*N,N:2*N],dNHy)

    # calculation of the corresponding power transported
    # assumption of of non structered materials in 1st and last layer, i.e.
    # const in x-direction.

    #field expansion coefficients of the regions the poyntingvector shall
    #beevaluated

    zeropt = np.array([0.0])

    eycoeffs1 = self.layer_list[0].eycoeffs
    hxcoeffs1 = np.zeros([N,N],complex)
    for m in range(0,N):
      hxcoeffs1[:,m] = -self.layer_list[0].kzvalsEy[m]/self.omega/self.layer_list[0].mustruct.xx(zeropt)*eycoeffs1[:,m]

    eycoeffsN = self.layer_list[self.num_of_layers-1].eycoeffs
    hxcoeffsN = np.zeros([N,N],complex)
    for m in range(0,N):
      hxcoeffsN[:,m] = -self.layer_list[self.num_of_layers-1].kzvalsEy[m]/self.omega/self.layer_list[self.num_of_layers-1].mustruct.xx(zeropt)*eycoeffsN[:,m]

    hycoeffs1 = self.layer_list[0].hycoeffs
    excoeffs1 = np.zeros([N,N],complex)
    for m in range(0,N):
      excoeffs1[:,m] = self.layer_list[0].kzvalsHy[m]/self.omega/self.layer_list[0].epsstruct.xx(zeropt)*hycoeffs1[:,m]

    hycoeffsN = self.layer_list[self.num_of_layers-1].hycoeffs
    excoeffsN = np.zeros([N,N],complex)
    for m in range(0,N):
      excoeffsN[:,m] = self.layer_list[self.num_of_layers-1].kzvalsHy[m]/self.omega/self.layer_list[self.num_of_layers-1].epsstruct.xx(zeropt)*hycoeffsN[:,m]

    # from incedent fields
    p0incEy = -0.5 * (np.dot(np.dot(np.dot(u1Ey.conj().T,hxcoeffs1.conj().T),eycoeffs1),u1Ey)).real
    p0incHy = 0.5 * (np.dot(np.dot(np.dot(u1Hy.conj().T,excoeffs1.conj().T),hycoeffs1),u1Hy)).real
    pEndincEy = 0.5 * (np.dot(np.dot(np.dot(dNEy.conj().T,hxcoeffsN.conj().T),eycoeffsN),dNEy)).real
    pEndincHy = -0.5 * (np.dot(np.dot(np.dot(dNHy.conj().T,excoeffsN.conj().T),hycoeffsN),dNHy)).real

    # from scattered fields
    p0scatEy = 0.5 * (np.dot(np.dot(np.dot(d1Ey.conj().T,hxcoeffs1.conj().T),eycoeffs1),d1Ey)).real
    p0scatHy = -0.5 * (np.dot(np.dot(np.dot(d1Hy.conj().T,excoeffs1.conj().T),hycoeffs1),d1Hy)).real
    pEndscatEy = -0.5 * (np.dot(np.dot(np.dot(uNEy.conj().T,hxcoeffsN.conj().T),eycoeffsN),uNEy)).real
    pEndscatHy = 0.5 * (np.dot(np.dot(np.dot(uNHy.conj().T,excoeffsN.conj().T),hycoeffsN),uNHy)).real

    return (p0incEy,p0incHy,pEndincEy,pEndincHy,p0scatEy,p0scatHy,pEndscatEy,pEndscatHy)

  def get_power_data_by_orders(self):
    if not(self._smat_calculated):
      self.Smat = self._calculate_smatrix()

    u1Ey, u1Hy = self._calculate_up_input_ampls(self.uinfuncEy,self.uinfuncHy,0.0)
    dNEy, dNHy = self._calculate_down_input_ampls(self.dinfuncEy,self.dinfuncHy,0.0)

    N = self.numpts

    #calculation of the amplitude vectors

    uNEy = np.dot(self._SmatEy[0:N,0:N],u1Ey) + np.dot(self._SmatEy[0:N,N:2*N],dNEy)
    d1Ey = np.dot(self._SmatEy[N:2*N,0:N],u1Ey) + np.dot(self._SmatEy[N:2*N,N:2*N],dNEy)

    uNHy = np.dot(self._SmatHy[0:N,0:N],u1Hy) + np.dot(self._SmatHy[0:N,N:2*N],dNHy)
    d1Hy = np.dot(self._SmatHy[N:2*N,0:N],u1Hy) + np.dot(self._SmatHy[N:2*N,N:2*N],dNHy)

    # calculation of the corresponding power transported
    # assumption of of non structered materials in 1st and last layer, i.e.
    # const in x-direction.

    #field expansion coefficients of the regions the poyntingvector shall
    #beevaluated

    zeropt = np.array([0.0])

    eycoeffs1 = self.layer_list[0].eycoeffs
    hxcoeffs1 = np.zeros([N,N],complex)
    for m in range(0,N):
      hxcoeffs1[:,m] = -self.layer_list[0].kzvalsEy[m]/self.omega/self.layer_list[0].mustruct.xx(zeropt)*eycoeffs1[:,m]

    eycoeffsN = self.layer_list[self.num_of_layers-1].eycoeffs
    hxcoeffsN = np.zeros([N,N],complex)
    for m in range(0,N):
      hxcoeffsN[:,m] = -self.layer_list[self.num_of_layers-1].kzvalsEy[m]/self.omega/self.layer_list[self.num_of_layers-1].mustruct.xx(zeropt)*eycoeffsN[:,m]

    hycoeffs1 = self.layer_list[0].hycoeffs
    excoeffs1 = np.zeros([N,N],complex)
    for m in range(0,N):
      excoeffs1[:,m] = self.layer_list[0].kzvalsHy[m]/self.omega/self.layer_list[0].epsstruct.xx(zeropt)*hycoeffs1[:,m]

    hycoeffsN = self.layer_list[self.num_of_layers-1].hycoeffs
    excoeffsN = np.zeros([N,N],complex)
    for m in range(0,N):
      excoeffsN[:,m] = self.layer_list[self.num_of_layers-1].kzvalsHy[m]/self.omega/self.layer_list[self.num_of_layers-1].epsstruct.xx(zeropt)*hycoeffsN[:,m]

    difforder_data = Diffraction_order_data(N)

    for i in range(0,N):
      difforder_data.p0scatEy[i] = 0.5 * (np.dot(d1Ey.conj().T,hxcoeffs1[i,:].conj().T)*(np.dot(eycoeffs1[i,:],d1Ey))).real
      difforder_data.p0scatHy[i] = -0.5 * (np.dot(d1Hy.conj().T,excoeffs1[i,:].conj().T)*(np.dot(hycoeffs1[i,:],d1Hy))).real
      difforder_data.pEndscatEy[i] = -0.5 * (np.dot(uNEy.conj().T,hxcoeffsN[i,:].conj().T)*(np.dot(eycoeffsN[i,:],uNEy))).real
      difforder_data.pEndscatHy[i] = 0.5 * (np.dot(uNHy.conj().T,excoeffsN[i,:].conj().T)*(np.dot(hycoeffsN[i,:],uNHy))).real
      if abs((self.layer_list[0].kxvec[i]+self.kappa)/self.omega/np.sqrt(self.layer_list[0].epsstruct.xx(zeropt))) > 1:
        difforder_data.is_evanecent0[i] = True
      else:
        difforder_data.is_evanecent0[i] = False
      if abs((self.layer_list[self.num_of_layers-1].kxvec[i]+self.kappa)/self.omega/np.sqrt(self.layer_list[self.num_of_layers-1].epsstruct.xx(zeropt))) > 1:
        difforder_data.is_evanecentend[i] = True
      else:
        difforder_data.is_evanecentend[i] = False

      difforder_data.outangle0[i] = (180./pi*np.arcsin((self.layer_list[0].kxvec[i]+self.kappa)/self.omega/np.sqrt(self.layer_list[0].epsstruct.xx(zeropt)))).real
      difforder_data.outangleend[i] = (180./pi*np.arcsin((self.layer_list[self.num_of_layers-1].kxvec[i]+self.kappa)/self.omega/np.sqrt(self.layer_list[self.num_of_layers-1].epsstruct.xx(zeropt)))).real

    return difforder_data

  def remove_single_layer(self,idx = None):
    if idx == None:
      self.layer_list.pop()
      self.num_of_layers -= 1
    elif idx >= 0 and idx < self.num_of_layers:
      self.layer_list.pop(idx)
      self.num_of_layers -= 1
    else:
      raise IndexError('There is no layer of the given idx.')
    self._set_calculation_status_false()

  def update_single_layer(self,idx,thickness,epsstruct,mustruct):
    if idx >= 0 and idx < self.num_of_layers:
      self.layer_list[i].thickness = thickness
      self.layer_list[i].thickness = epsstruct
      self.layer_list[i].thickness = mustruct
    else:
      raise IndexError('There is no layer of the given idx.')
    self._set_calculation_status_false()

  def _set_calculation_status_false(self):
    self._smat_calculated = False
    self._all_ampls_calculated = False
    self._reset_amplitudes()
    self._reset_Smat()

  def _reset_amplitudes(self):
    self._uamplsEy = []
    self._uamplsHy = []
    self._damplsEy = []
    self._damplsHy = []
    for i in range(0,self.num_of_layers):
      self._uamplsEy.append(np.zeros(self.numpts,complex))
      self._uamplsHy.append(np.zeros(self.numpts,complex))
      self._damplsEy.append(np.zeros(self.numpts,complex))
      self._damplsHy.append(np.zeros(self.numpts,complex))

  def _reset_Smat(self):
    self._SmatEy = None
    self._SmatHy = None

  def _calculate_interface_matrix(self,layer1,layer2):
    xvec = layer1.xvec

    #First for Ey fields

    muXX1 = layer1.mustruct.xx(xvec)
    muXX2 = layer2.mustruct.xx(xvec)

    EyExpansionMat1 = np.zeros([layer1.numpts,layer1.numpts],complex)
    EyExpansionMat2 = np.zeros([layer2.numpts,layer2.numpts],complex)

    for i in range(0,layer1.numpts):
      EyExpansionMat1[:,i] = np.fft.ifft(np.fft.fftshift(layer1.eycoeffs[:,i]))*layer1.numpts
      EyExpansionMat2[:,i] = np.fft.ifft(np.fft.fftshift(layer2.eycoeffs[:,i]))*layer1.numpts

    L1 = np.zeros([2*layer1.numpts,2*layer1.numpts],complex)
    L2 = np.zeros([2*layer2.numpts,2*layer2.numpts],complex)

    L1[0:layer1.numpts,0:layer1.numpts] = EyExpansionMat1
    L2[0:layer2.numpts,0:layer2.numpts] = EyExpansionMat2

    L1[0:layer1.numpts,layer1.numpts:2*layer1.numpts] = EyExpansionMat1
    L2[0:layer2.numpts,layer2.numpts:2*layer2.numpts] = EyExpansionMat2

    for i in range(0,layer1.numpts):
      L1[layer1.numpts:2*layer1.numpts,i] = layer1.kzvalsEy[i]/layer1.omega/muXX1*EyExpansionMat1[:,i]
      L2[layer2.numpts:2*layer2.numpts,i] = layer2.kzvalsEy[i]/layer2.omega/muXX2*EyExpansionMat2[:,i]

      L1[layer1.numpts:2*layer1.numpts,layer1.numpts+i] = -layer1.kzvalsEy[i]/layer1.omega/muXX1*EyExpansionMat1[:,i]
      L2[layer2.numpts:2*layer2.numpts,layer2.numpts+i] = -layer2.kzvalsEy[i]/layer2.omega/muXX2*EyExpansionMat2[:,i]

      #L1[layer1.numpts:2*layer1.numpts,layer1.numpts+i] = -L1[layer1.numpts:2*layer1.numpts,i]
      #L2[layer2.numpts:2*layer2.numpts,layer2.numpts+i] = -L2[layer2.numpts:2*layer2.numpts,i]

    del EyExpansionMat1, EyExpansionMat2
    del muXX1,muXX2

    transferI_Ey = np.linalg.solve(L2,L1)

    #Second for Hy fields
    epsXX1 = layer1.epsstruct.xx(xvec)
    epsXX2 = layer2.epsstruct.xx(xvec)

    HyExpansionMat1 = np.zeros([layer1.numpts,layer1.numpts],complex)
    HyExpansionMat2 = np.zeros([layer2.numpts,layer2.numpts],complex)

    for i in range(0,layer1.numpts):
      HyExpansionMat1[:,i] = np.fft.ifft(np.fft.fftshift(layer1.hycoeffs[:,i]))*layer1.numpts
      HyExpansionMat2[:,i] = np.fft.ifft(np.fft.fftshift(layer2.hycoeffs[:,i]))*layer1.numpts

    L1 = np.zeros([2*layer1.numpts,2*layer1.numpts],complex)
    L2 = np.zeros([2*layer2.numpts,2*layer2.numpts],complex)

    L1[0:layer1.numpts,0:layer1.numpts] = HyExpansionMat1
    L2[0:layer2.numpts,0:layer2.numpts] = HyExpansionMat2

    L1[0:layer1.numpts,layer1.numpts:2*layer1.numpts] = HyExpansionMat1
    L2[0:layer2.numpts,layer2.numpts:2*layer2.numpts] = HyExpansionMat2

    for i in range(0,layer1.numpts):
      L1[layer1.numpts:2*layer1.numpts,i] = layer1.kzvalsHy[i]/layer1.omega/epsXX1*HyExpansionMat1[:,i]
      L2[layer2.numpts:2*layer2.numpts,i] = layer2.kzvalsHy[i]/layer2.omega/epsXX2*HyExpansionMat2[:,i]

      L1[layer1.numpts:2*layer1.numpts,layer1.numpts+i] = -layer1.kzvalsHy[i]/layer1.omega/epsXX1*HyExpansionMat1[:,i]
      L2[layer2.numpts:2*layer2.numpts,layer2.numpts+i] = -layer2.kzvalsHy[i]/layer2.omega/epsXX2*HyExpansionMat2[:,i]

      #L1[layer1.numpts:2*layer1.numpts,layer1.numpts+i] = -L1[layer1.numpts:2*layer1.numpts,i]
      #L2[layer2.numpts:2*layer2.numpts,layer2.numpts+i] = -L2[layer2.numpts:2*layer2.numpts,i]

    del HyExpansionMat1, HyExpansionMat2
    del epsXX1,epsXX2

    transferI_Hy = np.linalg.solve(L2,L1)

    del L1,L2

    return (transferI_Ey,transferI_Hy)

  def _s_star_prop_prod(self,S1,kzvals,z):
    Stot = np.zeros([np.size(S1,0),np.size(S1,1)],complex)
    N = np.size(S1,0)/2

    if N != np.size(kzvals):
      raise IndexError('Matrix dimension and number of propagation constants mismatch!')

    S2sub = np.diag(np.exp(1.0j * kzvals*z))

    Stot[0:N,0:N]     = np.dot(S2sub,S1[0:N,0:N])
    Stot[0:N,N:2*N]   = np.dot(np.dot(S2sub,S1[0:N,N:2*N]),S2sub)
    Stot[N:2*N,0:N]   = S1[N:2*N,0:N]
    Stot[N:2*N,N:2*N] = np.dot(S1[N:2*N,N:2*N],S2sub)

    return Stot

  def _s_star_s_prod(self,S1,S2):
    Stot = np.zeros([np.size(S1,0),np.size(S1,1)],complex)
    N = np.size(S1,0)/2

    S_tmp_ad_inv = np.eye(N) - np.dot(S1[0:N,N:2*N],S2[N:2*N,0:N])

    Stot[0:N,0:N]     = np.dot(S2[0:N,0:N],np.linalg.solve(S_tmp_ad_inv,S1[0:N,0:N]))
    Stot[0:N,N:2*N]   = np.dot(np.dot(S2[0:N,0:N],np.linalg.solve(S_tmp_ad_inv,S1[0:N,N:2*N])),S2[N:2*N,N:2*N]) + S2[0:N,N:2*N]
    Stot[N:2*N,0:N]   = S1[N:2*N,0:N] + np.dot(np.dot(S1[N:2*N,N:2*N],S2[N:2*N,0:N]),np.linalg.solve(S_tmp_ad_inv,S1[0:N,0:N]))
    Stot[N:2*N,N:2*N] = np.dot(np.dot(S1[N:2*N,N:2*N],np.dot(S2[N:2*N,0:N],np.linalg.solve(S_tmp_ad_inv,S1[0:N,N:2*N]))+np.eye(N)),S2[N:2*N,N:2*N])

    return Stot

  def _s_star_t_prod(self,S,Tin):
    Stot = np.zeros([np.size(S,0),np.size(S,1)],complex)
    N = np.size(S,0)/2

    Tmp_ad_inv_tr = (np.dot(Tin[N:2*N,0:N],S[0:N,N:2*N]) + Tin[N:2*N,N:2*N]).T

    Stot[0:N,N:2*N]   = (np.linalg.solve(Tmp_ad_inv_tr,(np.dot(Tin[0:N,0:N],S[0:N,N:2*N]) + Tin[0:N,N:2*N]).T)).T
    Stot[0:N,0:N]     = np.dot(Tin[0:N,0:N],S[0:N,0:N]) - np.dot(np.dot(Stot[0:N,N:2*N],Tin[N:2*N,0:N]),S[0:N,0:N])
    Stot[N:2*N,N:2*N] = (np.linalg.solve(Tmp_ad_inv_tr,(S[N:2*N,N:2*N]).T)).T
    Stot[N:2*N,0:N]   = S[N:2*N,0:N] - np.dot(np.dot(Stot[N:2*N,N:2*N],Tin[N:2*N,0:N]),S[0:N,0:N])

    return Stot

  def _calculate_down_input_ampls(self,infuncEy,infuncHy,z0):
    if self.num_of_layers < 1:
      raise IndexError('number of layers must be at least one!')

    xvec = self.layer_list[self.num_of_layers-1].xvec

    kzvalsEy = self.layer_list[self.num_of_layers-1].kzvalsEy
    kzvalsHy = self.layer_list[self.num_of_layers-1].kzvalsHy

    Eyinc = infuncEy(xvec)
    Hyinc = infuncHy(xvec)

    EyExpansionMat = np.zeros([self.numpts,self.numpts],complex)
    HyExpansionMat = np.zeros([self.numpts,self.numpts],complex)

    for i in range(0,self.numpts):
      EyExpansionMat[:,i] = np.exp(- 1.0j * kzvalsEy[i]*z0)*np.exp( 1.0j * self.kappa*xvec)*np.fft.ifft(np.fft.fftshift(self.layer_list[self.num_of_layers-1].eycoeffs[:,i]))*self.numpts
      HyExpansionMat[:,i] = np.exp(- 1.0j * kzvalsHy[i]*z0)*np.exp( 1.0j * self.kappa*xvec)*np.fft.ifft(np.fft.fftshift(self.layer_list[self.num_of_layers-1].hycoeffs[:,i]))*self.numpts

    ampl_vecEy = np.linalg.solve(EyExpansionMat,Eyinc)
    ampl_vecHy = np.linalg.solve(HyExpansionMat,Hyinc)

    return (ampl_vecEy,ampl_vecHy)

  def _calculate_up_input_ampls(self,infuncEy,infuncHy,z0):
    if self.num_of_layers < 1:
      raise IndexError('number of layers must be at least one!')

    xvec = self.layer_list[0].xvec

    kzvalsEy = self.layer_list[0].kzvalsEy
    kzvalsHy = self.layer_list[0].kzvalsHy

    Eyinc = infuncEy(xvec)
    Hyinc = infuncHy(xvec)

    EyExpansionMat = np.zeros([self.numpts,self.numpts],complex)
    HyExpansionMat = np.zeros([self.numpts,self.numpts],complex)

    for i in range(0,self.numpts):
      EyExpansionMat[:,i] = np.exp( 1.0j * kzvalsEy[i]*z0)*np.exp( 1.0j * self.kappa*xvec)*np.fft.ifft(np.fft.fftshift(self.layer_list[0].eycoeffs[:,i]))*self.numpts
      HyExpansionMat[:,i] = np.exp( 1.0j * kzvalsHy[i]*z0)*np.exp( 1.0j * self.kappa*xvec)*np.fft.ifft(np.fft.fftshift(self.layer_list[0].hycoeffs[:,i]))*self.numpts

    ampl_vecEy = np.linalg.solve(EyExpansionMat,Eyinc)
    ampl_vecHy = np.linalg.solve(HyExpansionMat,Hyinc)

    return (ampl_vecEy,ampl_vecHy)

class Diffraction_order_data():
  def __init__(self,numpts):
    self.order_number = range(-numpts/2,numpts/2)
    self.p0scatEy = np.zeros(numpts)
    self.p0scatHy = np.zeros(numpts)
    self.pEndscatEy = np.zeros(numpts)
    self.pEndscatHy = np.zeros(numpts)

    self.outangle0 = np.zeros(numpts)
    self.outangleend = np.zeros(numpts)

    self.is_evanecent0 = []
    self.is_evanecentend = []
    for i in range(0,numpts):
      self.is_evanecent0.append(True)
      self.is_evanecentend.append(True)

class Material_struct():
  def __init__(self,funcxx,funcyy,funczz):
    self.xx = funcxx
    self.yy = funcyy
    self.zz = funczz
