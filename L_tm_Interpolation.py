def _L_tm_Interpolation(self):

    zm = np.array([0.0, 1.0/3.0, 2.0/3.0, 1.0])
    B = np.matrix([np.ones(np.size(zm)),zm,1.0/2.0*np.power(zm,2),1.0/3.0*np.power(zm,3)])
    L_tm_list=[[],[],[],[]]
    for zi in range(0,len(zm)): #Erzeugung der Liste der L_te(z_i)
        muXX = self.mustruct.xx(self.xvec,zm[zi])
        muZZ = self.mustruct.zz(self.xvec,zm[zi])
        epsYY = self.epsstruct.yy(self.xvec,zm[zi])

        pdmat = 1j*np.diag(self.kxvec+self.kappa)
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

        L_tm = np.dot(muXXmat,(epsYYmat*self.omega*self.omega + np.dot(np.dot(pdmat,invmuZZmat),pdmat)))
        L_tm_list[zi]=L_tm

    L_tm_orders= np.zeros((len(zm),self.numpts,self.numpts),complex)
    lpq2 = np.zeros(len(zm),complex)

    for p in range(0,self.numpts): #Berechnung der Interpolationskoeffizienten
        for q in range(0,self.numpts):
            for i in range(0,len(zm)):
                lpq2[i]=L_te_list[i][p,q]
            lpq1=np.linalg.solve(B,lpq2)
            for o in range(0,len(zm)):
                L_tm_orders[o][p,q] = lpq1[o]

    return L_tm_orders
