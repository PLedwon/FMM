def _L_te_Interpolation(self):

    zm = np.array([0.0, 1.0/3.0, 2.0/3.0, 1.0])
    B = np.matrix([np.ones(np.size(zm)),zm,1.0/2.0*np.power(zm,2),1.0/3.0*np.power(zm,3)])
    L_te_list=[[],[],[],[]]
    for zi in range(0,len(zm)): #Erzeugung der Liste der L_te(z_i)
        epsXX = self.epsstruct.xx(self.xvec,zm[zi])
        epsZZ = self.epsstruct.zz(self.xvec,zm[zi])
        muYY = self.mustruct.yy(self.xvec,zm[zi])

        pdmat = 1j*np.diag(self.kxvec+self.kappa)
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

        L_te = np.dot(epsXXmat,(muYYmat*self.omega*self.omega + np.dot(np.dot(pdmat,invepsZZmat),pdmat)))
        L_te_list[zi]=L_te

    L_te_orders= np.zeros((len(zm),self.numpts,self.numpts),complex)
    lpq2 = np.zeros(len(zm),complex)

    for p in range(0,self.numpts): #Berechnung der Interpolationskoeffizienten
        for q in range(0,self.numpts):
            for i in range(0,len(zm)):
                lpq2[i]=L_te_list[i][p,q]
            lpq1=np.linalg.solve(B,lpq2)
            for o in range(0,len(zm)):
                L_te_orders[o][p,q] = lpq1[o]

    return L_te_orders
